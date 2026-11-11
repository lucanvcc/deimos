"""Bit-identical fast BLS kernel: search plug-in "deimos.fastbls:exact_select".

exact_select is the standard search (frozen_select) with astropy's BLS kernel
(astropy.timeseries.periodograms.bls.methods.bls_fast -> bls.c run_bls) swapped for the C kernel
below while the search runs. The kernel reproduces astropy 6.1's run_bls bit for bit: same
arithmetic (including the fused multiply-adds clang emits for the shipped wheel), same per-bin
summation order, same tie-breaking. It is faster because it
  1. computes y*ivar's operands and t - min_t once per call, not once per trial period;
  2. computes phase-bin indices in a separate vectorized pass, then scatters both bin sums with
     one 2-lane fma (fma(1.0, ivar, s) rounds exactly as s + ivar);
  3. prefilters the duration/phase scan without division or sqrt: a block of cells is skipped only
     when a cross-multiplied bound with a 1e-12 relative rounding margin proves that astropy would
     select none of them (y_out < y_in strictly, or snr strictly below the running best). Every
     cell that can be selected is evaluated with astropy's exact expressions.
On first load the compiled kernel is checked bitwise against astropy on synthetic data and refuses
to run on any difference. The C source is part of this file, so it is covered by code_hash().

Two variants, picked by platform.machine(): x86-64 uses an unfused multiply-add in binning and -ffp-contract=off
(matches the x86-64 astropy wheel, whose baseline has no FMA); everything else uses fma and -ffp-contract=on.
The variant is part of the C source, CFLAGS and the .so name/hash key.
"""
import ctypes, hashlib, os, pathlib, platform, subprocess, time
import numpy as np
from astropy.timeseries.periodograms.bls import methods as _bls_methods
from .core import clean_frame, extract_multi, diagnostics, serial

X86 = platform.machine() in ('x86_64', 'AMD64')
BUILD_DIR = pathlib.Path(os.environ.get('DEIMOS_BUILD_DIR') or '~/.cache/deimos').expanduser()

C_SOURCE = r'''
#include <math.h>
#include <float.h>
#include <stdlib.h>
#if defined(__aarch64__)
#include <arm_neon.h>
#endif

#define SCAN_BLOCK 64
#define INDEX_BLOCK 512

/* astropy bls.c compute_objective, verbatim. */
static void compute_objective(double y_in, double y_out, double ivar_in, double ivar_out, int obj_flag,
                              double* objective, double* log_likelihood, double* depth, double* depth_err, double* depth_snr) {
    if (obj_flag) {
        double arg = y_out - y_in;
        *log_likelihood = 0.5*ivar_in*arg*arg;
        *objective = *log_likelihood;
    } else {
        *depth = y_out - y_in;
        *depth_err = sqrt(1.0 / ivar_in + 1.0 / ivar_out);
        *depth_snr = *depth / *depth_err;
        *objective = *depth_snr;
    }
}

/* astropy bls.c wrap_into, verbatim (clang contracts it to one fnmsub, as in the wheel). */
static inline double wrap_into(double x, double period) {
    return x - period * floor(x / period);
}

/* Same signature and outputs as astropy bls.c run_bls. */
int fast_bls(int N, const double* t, const double* y, const double* ivar,
             int n_periods, const double* periods, int n_durations, const double* durations,
             int oversample, int obj_flag,
             double* best_objective, double* best_depth, double* best_depth_err, double* best_duration,
             double* best_phase, double* best_depth_snr, double* best_log_like) {
    double max_period = periods[0], min_period = periods[0];
    int k;
    for (k = 1; k < n_periods; ++k) {
        if (periods[k] < min_period) min_period = periods[k];
        if (periods[k] > max_period) max_period = periods[k];
    }
    if (min_period < DBL_EPSILON) return 1;
    double min_duration = durations[0], max_duration = durations[0];
    for (k = 1; k < n_durations; ++k) {
        if (durations[k] < min_duration) min_duration = durations[k];
        if (durations[k] > max_duration) max_duration = durations[k];
    }
    if ((max_duration > min_period) || (min_duration < DBL_EPSILON)) return 2;

    double bin_duration = min_duration / ((double)oversample);
    int max_n_bins = (int)(ceil(max_period / bin_duration)) + oversample;
    int blocksize = max_n_bins + 1;
    int npts = N > 0 ? N : 1;

    /* Bin n holds (sum y*ivar, sum ivar) at bins[2n], bins[2n+1]; yone holds (y, 1.0) per point. */
    double* bins = (double*)malloc(2*blocksize*sizeof(double));
    double* yone = (double*)malloc(2*npts*sizeof(double));
    double* tt = (double*)malloc(npts*sizeof(double));
    int* ind = (int*)malloc(INDEX_BLOCK*sizeof(int));
    int* durs = (int*)malloc((n_durations > 0 ? n_durations : 1)*sizeof(int));
    double* cy = (double*)malloc(blocksize*sizeof(double));
    double* ci = (double*)malloc(blocksize*sizeof(double));
    if (!bins || !yone || !tt || !ind || !durs || !cy || !ci) {
        free(bins); free(yone); free(tt); free(ind); free(durs); free(cy); free(ci);
        return -2;
    }

    double min_t = INFINITY;
    double sum_y = 0.0, sum_ivar = 0.0;
    int i;
    for (i = 0; i < N; ++i) {
        min_t = fmin(min_t, t[i]);
        sum_y += y[i] * ivar[i];
        sum_ivar += ivar[i];
    }
    for (i = 0; i < N; ++i) {
        yone[2*i] = y[i];
        yone[2*i+1] = 1.0;
        tt[i] = t[i] - min_t;
    }
    for (k = 0; k < n_durations; ++k) durs[k] = (int)(round(durations[k] / bin_duration));

    const double margin = 1e-12, shrink = 1.0 - 1e-9;
    int p;
    for (p = 0; p < n_periods; ++p) {
        double period = periods[p];
        int n_bins = (int)(ceil(period / bin_duration)) + oversample;
        int n, n0, j;
        for (n = 0; n < 2*(n_bins+1); ++n) bins[n] = 0.0;

        /* Bin the data. astropy's "mean_y += y*ivar" is contracted to one fma; fma(1.0, ivar, s) rounds
         * exactly as s + ivar. Indices are computed per L1-sized block, then scattered in point order,
         * so every bin is summed in the same order as astropy. */
        for (n0 = 0; n0 < N; n0 += INDEX_BLOCK) {
            int cnt = N - n0 < INDEX_BLOCK ? N - n0 : INDEX_BLOCK, c;
            const double* tb = tt + n0;
            for (c = 0; c < cnt; ++c)
                ind[c] = (int)(wrap_into(tb[c], period) / bin_duration) + 1;
#if defined(__aarch64__)
            const double* yb = yone + 2*n0;
            const double* ib = ivar + n0;
            for (c = 0; c < cnt; ++c) {
                double* b = bins + 2*ind[c];
                vst1q_f64(b, vfmaq_n_f64(vld1q_f64(b), vld1q_f64(yb + 2*c), ib[c]));
            }
#else
            for (c = 0; c < cnt; ++c) {
                int jj = ind[c];
                bins[2*jj] = @BIN@;
                bins[2*jj+1] += ivar[n0+c];
            }
#endif
        }

        for (n = 1, j = n_bins - oversample; n <= oversample; ++n, ++j) {
            bins[2*j] = bins[2*n];
            bins[2*j+1] = bins[2*n+1];
        }
        for (n = 1; n <= n_bins; ++n) {
            bins[2*n] += bins[2*n-2];
            bins[2*n+1] += bins[2*n-1];
        }
        for (n = 0; n <= n_bins; ++n) {
            cy[n] = bins[2*n];
            ci[n] = bins[2*n+1];
        }

        double objective, log_like, depth, depth_err, depth_snr;
        double best = -INFINITY;
        best_objective[p] = -INFINITY;
        for (k = 0; k < n_durations; ++k) {
            int dur = durs[k];
            int n_max = n_bins-dur;
            for (n0 = 0; n0 <= n_max; n0 += SCAN_BLOCK) {
                int n1 = n0 + SCAN_BLOCK - 1 < n_max ? n0 + SCAN_BLOCK - 1 : n_max;
                if (!obj_flag) {
                    /* Branch-free block prefilter against the best value at block start (best only grows).
                     * depth*Iin*Iout = a - b before rounding; m bounds all rounding with a 1e-12 relative margin.
                     * A cell is dropped only if y_out < y_in or its snr < best: astropy would not select it. */
                    double bb = best > 0 ? best * best * shrink : 0.0;
                    int any = 0;
                    for (n = n0; n <= n1; ++n) {
                        double Yin = cy[n+dur] - cy[n];
                        double Iin = ci[n+dur] - ci[n];
                        double Yout = sum_y - Yin;
                        double Iout = sum_ivar - Iin;
                        double a = Yout * Iin, b = Yin * Iout;
                        double D = a - b, m = margin * (fabs(a) + fabs(b));
                        double Dp = D + m, Q = Iin * Iout * (Iin + Iout);
                        any |= !(Iin < DBL_EPSILON) & !(Iout < DBL_EPSILON) & !(D < -m) & !(Dp * Dp < bb * Q);
                    }
                    if (!any) continue;
                }
                for (n = n0; n <= n1; ++n) {
                    /* astropy's loop body, verbatim. */
                    double y_in = cy[n+dur] - cy[n];
                    double ivar_in = ci[n+dur] - ci[n];
                    double y_out = sum_y - y_in;
                    double ivar_out = sum_ivar - ivar_in;
                    if ((ivar_in < DBL_EPSILON) || (ivar_out < DBL_EPSILON)) continue;
                    y_in /= ivar_in;
                    y_out /= ivar_out;
                    compute_objective(y_in, y_out, ivar_in, ivar_out, obj_flag, &objective, &log_like, &depth, &depth_err, &depth_snr);
                    if (y_out >= y_in && objective > best_objective[p]) {
                        best = best_objective[p] = objective;
                        compute_objective(y_in, y_out, ivar_in, ivar_out, (obj_flag == 0), &objective, &log_like, &depth, &depth_err, &depth_snr);
                        best_depth[p]     = depth;
                        best_depth_err[p] = depth_err;
                        best_depth_snr[p] = depth_snr;
                        best_log_like[p]  = log_like;
                        best_duration[p]  = dur * bin_duration;
                        best_phase[p]     = fmod(n*bin_duration + 0.5*best_duration[p] + min_t, period);
                    }
                }
            }
        }
    }

    free(bins); free(yone); free(tt); free(ind); free(durs); free(cy); free(ci);
    return 0;
}
'''.replace('@BIN@', 'y[n0+c] * ivar[n0+c] + bins[2*jj]' if X86 else 'fma(y[n0+c], ivar[n0+c], bins[2*jj])')

# fma variant: clang's "on" is the macOS wheel's contraction; gcc (Linux aarch64 wheel) spells it "fast".
_CONTRACT = 'off' if X86 else 'on' if platform.system() == 'Darwin' else 'fast'
CFLAGS = ['-O3', f'-ffp-contract={_CONTRACT}', '-shared', '-fPIC']
ASTROPY_BLS_FAST = _bls_methods.bls_fast
_LIB = None


def _compiler():
    return os.environ.get('CC', 'cc')


def build():
    """Compile the kernel once into BUILD_DIR, keyed by source, flags and compiler version."""
    version = subprocess.run([_compiler(), '--version'], capture_output=True, text=True, check=True).stdout
    key = hashlib.sha256('\0'.join([C_SOURCE, ' '.join(CFLAGS), version]).encode()).hexdigest()[:16]
    out = BUILD_DIR / f"fastbls{'_x86' if X86 else ''}-{key}.so"
    if not out.exists():
        out.parent.mkdir(parents=True, exist_ok=True)
        src = out.with_name(f'{out.stem}.{os.getpid()}.c'); tmp = out.with_name(f'{out.stem}.{os.getpid()}.so')
        src.write_text(C_SOURCE)
        try:
            subprocess.run([_compiler(), *CFLAGS, '-o', str(tmp), str(src)], check=True, capture_output=True, text=True)
            tmp.replace(out)   # atomic: concurrent workers may build the same key
        finally:
            src.unlink(missing_ok=True); tmp.unlink(missing_ok=True)
    return out


def _load():
    global _LIB
    if _LIB is None:
        lib = ctypes.CDLL(str(build()))
        dp = np.ctypeslib.ndpointer(np.float64, flags='C_CONTIGUOUS')
        lib.fast_bls.argtypes = [ctypes.c_int, dp, dp, dp, ctypes.c_int, dp, ctypes.c_int, dp, ctypes.c_int, ctypes.c_int] + [dp] * 7
        lib.fast_bls.restype = ctypes.c_int
        _LIB = lib
        self_check()
    return _LIB


def bls_fast(t, y, ivar, period, duration, oversample, use_likelihood):
    """Drop-in for astropy methods.bls_fast (same arguments, same 7 output arrays, same errors)."""
    lib = _LIB or _load()
    t, y, ivar, period, duration = (np.ascontiguousarray(a, dtype=np.float64) for a in (t, y, ivar, period, duration))
    outs = [np.empty_like(period) for _ in range(7)]
    flag = lib.fast_bls(len(t), t, y, ivar, len(period), period, len(duration), duration, int(oversample), int(bool(use_likelihood)), *outs)
    if flag < 0:
        raise MemoryError()
    if flag > 0:
        raise ValueError("Invalid inputs for period and/or duration")
    return tuple(outs)


def _check_data(seed):
    rng = np.random.default_rng(seed)
    t = np.concatenate([100 + 95 * q + np.arange(1200) * 0.0204 for q in range(3)])
    t = np.sort(np.delete(t, rng.choice(len(t), 300, replace=False)))
    dy = rng.uniform(5e-4, 2e-3, len(t))
    y = rng.normal(0, 1e-3, len(t)) - 3e-3 * (np.abs((t - 101.3 + 1.7) % 3.4 - 1.7) < 0.06)
    periods = np.geomspace(0.5, 30, 4000)
    durations = np.linspace(0.05, 0.3, 8)
    return t - t.min(), y - np.median(y), 1.0 / dy ** 2, periods, durations


def self_check():
    """Refuse to run unless the compiled kernel matches astropy bit for bit (both objectives, all outputs)."""
    for seed, use_likelihood in ((1, False), (2, True)):
        args = _check_data(seed)
        ours = bls_fast(*args, 10, use_likelihood)
        ref = ASTROPY_BLS_FAST(*args, 10, use_likelihood)
        ok = np.isfinite(ref[0])   # astropy leaves the other outputs uninitialized where no cell was selected
        if not np.array_equal(np.isfinite(ours[0]), ok):
            raise RuntimeError('fastbls kernel is not bit-identical to astropy bls_fast; refusing to run')
        for a, b in zip(ours, ref):
            if not np.array_equal(a[ok].view(np.int64), b[ok].view(np.int64)):
                raise RuntimeError('fastbls kernel is not bit-identical to astropy bls_fast; refusing to run')


class kernel:
    """Context manager: route astropy BoxLeastSquares.power through the fast kernel."""

    def __enter__(self):
        _load()
        _bls_methods.bls_fast = bls_fast
        return self

    def __exit__(self, *exc):
        _bls_methods.bls_fast = ASTROPY_BLS_FAST


def exact_select(raw, cfg):
    """frozen_select with the bit-identical fast BLS kernel."""
    with kernel():
        t = time.monotonic(); clean = clean_frame(raw, cfg['clean']); multi = extract_multi(clean, cfg['search'], max_peaks=5, residual=True); peaks = list(multi['peaks'])
        if multi['residual'] and multi['residual'].get('status') == 'ok': peaks.append(multi['residual'])
        for p in peaks: p['diagnostics'] = diagnostics(clean, p, cfg['validate'])
        return serial(dict(peaks=peaks, harmonic_groups=multi['harmonic_groups'], seconds=time.monotonic() - t, rows=len(clean), baseline_days=float(np.ptp(clean.time))))
