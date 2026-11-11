"""Numerical core: the frozen clean, BLS search, vetting and discovery/holdout primitives.

Function bodies are pinned by frozen-pins.json."""
from __future__ import annotations
import hashlib, json
import numpy as np
import pandas as pd
from astropy.timeseries import BoxLeastSquares
from scipy.ndimage import median_filter

def clean_frame(raw: pd.DataFrame, settings: dict) -> pd.DataFrame:
    """Retain finite positive-error samples; normalize and flatten each gap segment."""
    raw = raw.replace([np.inf, -np.inf], np.nan).dropna(subset=['time', 'flux', 'flux_err'])
    raw = raw[(raw.flux_err > 0) & (raw.flux > 0)].sort_values('time').drop_duplicates('time')
    segments = []
    for quarter, group in raw.groupby('quarter'):
        for _, segment in group.groupby((group.time.diff() > settings['gap_days']).cumsum()):
            if len(segment) < 5:
                continue
            segment = segment.copy()
            median = np.median(segment.flux)
            segment['normalized_flux'] = segment.flux / median
            cadence = np.median(np.diff(segment.time))
            window = max(3, int(settings['flatten_window_days'] / cadence) | 1)
            window = min(window, len(segment) if len(segment) % 2 else len(segment) - 1)
            trend = median_filter(segment.normalized_flux.to_numpy(), size=window, mode='nearest')
            segment['trend'] = trend
            segment['flux'] = segment.normalized_flux / trend
            segment['flux_err'] = segment.flux_err / median / trend
            residual = segment.flux - 1
            noise = max(1.4826 * np.median(np.abs(residual - np.median(residual))), float(np.median(segment.flux_err)))
            keep = (residual < settings['sigma_upper'] * noise) & (residual > -settings['sigma_lower'] * noise)
            segments.append(segment[keep])
    if not segments:
        raise ValueError('No usable segments')
    clean = pd.concat(segments).sort_values('time').reset_index(drop=True)
    if len(clean) < settings['min_cadences']:
        raise ValueError('Insufficient cleaned cadences')
    return clean


def transit_masks(time: np.ndarray, period: float, epoch: float, duration: float) -> tuple:
    """Return centered phase in days, integer event index and in-transit mask."""
    cycles = np.floor((time - epoch) / period + 0.5).astype(int)
    phase = time - epoch - cycles * period
    return phase, cycles, np.abs(phase) < duration / 2


def observed_transits(events: np.ndarray, inside: np.ndarray, min_points: int) -> int:
    """Count sampled events only, including sparse/gapped observing windows."""
    _, counts = np.unique(events[inside], return_counts=True)
    return int(np.sum(counts >= min_points))


def extract_bls(frame: pd.DataFrame, settings: dict, retain_details: bool = False) -> tuple[dict, pd.DataFrame]:
    """Search a frequency grid limiting accumulated phase drift to a duration fraction.

    Errors are floored by robust point-to-point scatter. SNR is a white-noise
    depth/error statistic, not a trials-corrected detection significance.
    """
    t, y, dy = (frame[col].to_numpy(dtype=float) for col in ('time', 'flux', 'flux_err'))
    baseline = np.ptp(t)
    if baseline <= settings['min_period_days'] * 2:
        raise ValueError('Insufficient baseline')
    scatter = 1.4826 * np.median(np.abs(np.diff(y) - np.median(np.diff(y)))) / np.sqrt(2)
    dy = np.maximum(dy, scatter)
    model = BoxLeastSquares(t, y, dy)
    duration = np.linspace(settings['min_duration_days'], settings['max_duration_days'], settings['duration_steps'])
    pmin, pmax = settings['min_period_days'], settings['max_period_days']
    # df/f <= phase_drift_fraction * shortest_duration / baseline.
    step = np.log1p(settings['phase_drift_fraction'] * min(duration) / baseline)
    count = int(np.ceil(np.log(pmax / pmin) / step)) + 1
    periods = np.geomspace(pmin, pmax, count)
    best, grids = None, []
    for start in range(0, count, settings['chunk_size']):
        r = model.power(periods[start:start + settings['chunk_size']], duration, objective='snr')
        valid = np.isfinite(r.power)
        if not valid.any():
            continue
        i = int(np.nanargmax(np.where(valid, r.power, np.nan)))
        if best is None or float(r.depth_snr[i]) > best['snr']:
            best = {'period_days': float(r.period[i]), 'epoch_bkjd': float(r.transit_time[i]),
                    'duration_days': float(r.duration[i]), 'depth': float(r.depth[i]), 'snr': float(r.depth_snr[i])}
        grid = {'period_days': r.period, 'power_snr': r.power}
        if retain_details:
            grid.update(epoch_bkjd=r.transit_time, duration_days=r.duration, depth=r.depth, snr=r.depth_snr)
        grids.append(pd.DataFrame(grid))
    if best is None:
        raise ValueError('No finite BLS result')
    phase, events, inside = transit_masks(t, best['period_days'], best['epoch_bkjd'], best['duration_days'])
    best['number_of_transits'] = observed_transits(events, inside, settings['min_points_per_transit'])
    best['initial_candidate'] = bool(all(np.isfinite(list(best.values()))) and best['depth'] > 0 and best['depth'] < 1
                                     and best['snr'] >= settings['min_snr'] and best['number_of_transits'] >= settings['min_transits']
                                     and settings['min_duration_days'] - 1e-10 <= best['duration_days'] <= settings['max_duration_days'] + 1e-10
                                     and pmin <= best['period_days'] <= pmax)
    best.update(status='ok', grid_points=count, noise_floor=float(scatter), time_system='BKJD = BJD_TDB - 2454833')
    return best, pd.concat(grids, ignore_index=True)


def extract_multi(frame: pd.DataFrame, settings: dict, max_peaks: int = 5, residual: bool = True) -> dict:
    """Retain distinct local maxima on the unchanged BLS grid plus one masked search."""
    from scipy.signal import find_peaks
    best, grid = extract_bls(frame, settings, retain_details=True)
    powers = grid.power_snr.to_numpy()
    indices = list(find_peaks(powers)[0]) + [0, len(grid)-1, int(np.argmax(powers))]
    indices = sorted(set(indices), key=lambda i: powers[i], reverse=True)
    retained, grouped = [], []
    for i in indices:
        row = grid.iloc[i]
        duplicate = next((j for j,p in enumerate(retained) if any(
            abs(row.period_days/p['period_days']/ratio-1) <= .01 for ratio in (.25,.5,1.,2.,4.))), None)
        if duplicate is not None:
            grouped.append({'period_days': float(row.period_days), 'parent_rank': duplicate+1})
            continue
        peak = {k: float(row[k]) for k in ('period_days','epoch_bkjd','duration_days','depth','snr')}
        _,events,inside = transit_masks(frame.time.to_numpy(),peak['period_days'],peak['epoch_bkjd'],peak['duration_days'])
        peak['number_of_transits'] = observed_transits(events, inside, settings['min_points_per_transit'])
        peak['initial_candidate'] = bool(peak['snr'] >= settings['min_snr'] and 0 < peak['depth'] < 1 and peak['number_of_transits'] >= settings['min_transits'])
        peak.update(rank=len(retained)+1, provenance='original_discovery_grid')
        retained.append(peak)
        if len(retained) >= max_peaks: break
    extra = None
    if residual and best['initial_candidate']:
        _,_,mask = transit_masks(frame.time.to_numpy(), best['period_days'],best['epoch_bkjd'],best['duration_days']*1.5)
        remaining = frame.loc[~mask].copy()
        if len(remaining) >= 100:
            try:
                extra,_ = extract_bls(remaining,settings)
                extra.update(provenance='one_search_after_masking_strongest',masked_cadences=int(mask.sum()), mask_duration_factor=1.5,
                             parent_ephemeris={k: best[k] for k in ('period_days','epoch_bkjd','duration_days')})
            except ValueError as exc: extra = {'status':'insufficient_coverage','reason':str(exc)}
    return {'best':best,'peaks':retained,'harmonic_groups':grouped,'residual':extra,'optimization':'none; original BLS settings/grid'}


def smooth_modulation(frame: pd.DataFrame, candidate: dict) -> dict:
    """Flag broad baseline changes that can masquerade as transit-like structure."""
    result = {'status': 'manual_check_required', 'amplitude': None, 'amplitude_to_depth': None,
              'folded_modulation_amplitude': None, 'folded_amplitude_to_depth': None,
              'flagged': False, 'method': 'rolling baseline and transit-excluded phase-bin modulation; heuristic only'}
    if 'flux' not in frame:
        return result
    if 'time' in frame:
        t, y = frame.time.to_numpy(float), frame.flux.to_numpy(float)
    else:
        t, y = np.arange(len(frame), dtype=float), frame.flux.to_numpy(float)
    good = np.isfinite(t) & np.isfinite(y)
    t, y = t[good], y[good]
    if len(y) < 30:
        return result
    window_days = max(0.5, 5 * float(candidate.get('duration_days') or 0.1))
    cadence = np.nanmedian(np.diff(np.sort(t)))
    if not np.isfinite(cadence) or cadence <= 0:
        return result
    window = max(5, int(window_days / cadence))
    smooth = pd.Series(y).rolling(window, center=True, min_periods=max(3, window // 3)).median().to_numpy()
    smooth = smooth[np.isfinite(smooth)]
    if len(smooth) < 10:
        return result
    amplitude = float(np.nanpercentile(smooth, 95) - np.nanpercentile(smooth, 5))
    depth = abs(float(candidate.get('depth') or 0))
    ratio = amplitude / depth if depth > 0 else None
    folded_amp = folded_ratio = None
    if depth > 0 and all(k in candidate for k in ('period_days', 'epoch_bkjd', 'duration_days')):
        period = float(candidate['period_days'])
        phase = (t - float(candidate['epoch_bkjd']) + period / 2) % period - period / 2
        outside = np.abs(phase) > 1.5 * float(candidate['duration_days'])
        if outside.sum() >= 30:
            bins = np.floor((phase[outside] / period + .5) * 40).astype(int).clip(0, 39)
            folded = pd.DataFrame({'bin': bins, 'flux': y[outside]}).groupby('bin').flux.median()
            if len(folded) >= 12:
                folded_amp = float(folded.quantile(.95) - folded.quantile(.05))
                folded_ratio = folded_amp / depth
    result.update(status='screened', amplitude=amplitude, amplitude_to_depth=ratio,
                  folded_modulation_amplitude=folded_amp, folded_amplitude_to_depth=folded_ratio,
                  flagged=bool((ratio is not None and ratio >= 1.2) or (folded_ratio is not None and folded_ratio >= 1.2)),
                  threshold_rolling_amplitudes=1.2, threshold_folded_amplitudes=1.2, window_days=window_days)
    return result


def red_noise(frame: pd.DataFrame, candidate: dict | None = None) -> dict:
    """Compare elapsed-time binned scatter outside transits within continuous segments."""
    result = {'status': 'manual_check_required', 'ratio': None, 'bin_cadences': 10,
              'method': 'elapsed-time bins within cadence-contiguous segments; transits masked'}
    if not {'time', 'flux'}.issubset(frame):
        return result
    t, y = frame.time.to_numpy(float), frame.flux.to_numpy(float)
    good = np.isfinite(t) & np.isfinite(y)
    t, y = t[good], y[good]
    quarter = frame.quarter.to_numpy()[good] if 'quarter' in frame else np.zeros(len(t), dtype=int)
    if len(t) < 30:
        return result
    order = np.argsort(t); t, y, quarter = t[order], y[order], quarter[order]
    cadence = float(np.median(np.diff(t)))
    if not np.isfinite(cadence) or cadence <= 0:
        return result
    gap = (np.diff(t) > 2.5 * cadence) | (quarter[1:] != quarter[:-1])
    segment = np.r_[0, np.cumsum(gap)]
    if candidate and all(k in candidate for k in ('period_days', 'epoch_bkjd', 'duration_days')):
        _, _, in_transit = transit_masks(t, candidate['period_days'], candidate['epoch_bkjd'], candidate['duration_days'])
        keep = ~in_transit
        t, y, segment = t[keep], y[keep], segment[keep]
    if len(t) < 30:
        return result
    frame_bins, points = [], []
    for sid in np.unique(segment):
        ix = segment == sid
        if ix.sum() < 20:
            continue
        st, sy = t[ix], y[ix]
        width = 10 * cadence
        bin_id = np.floor((st - st[0]) / width).astype(int)
        tmp = pd.DataFrame({'bin': bin_id, 'flux': sy}).groupby('bin').flux.agg(['mean', 'count'])
        tmp = tmp[tmp['count'] >= 5]
        if len(tmp) < 3:
            continue
        # Remove segment median only; preserve within-segment correlated changes.
        frame_bins.extend((tmp['mean'] - np.median(sy)).to_list())
        points.extend(sy - np.median(sy))
    if len(frame_bins) < 6 or len(points) < 30:
        return result
    bins, pts = np.asarray(frame_bins), np.asarray(points)
    point = float(np.std(pts, ddof=1))
    observed = float(np.std(bins, ddof=1))
    expected = point / np.sqrt(10)
    ratio = observed / expected if expected > 0 else None
    result.update(status='screened', ratio=ratio, flagged=bool(ratio is not None and ratio > 2.0),
                  point_scatter=point, binned_scatter=observed, white_expectation=expected,
                  segments=int(len(np.unique(segment))))
    return result


def quarter_depth_epoch(time, flux, candidate: dict, min_snr: float = 3.0) -> dict:
    """Fit small epoch shifts around the supplied ephemeris and estimate depth."""
    t, y = np.asarray(time, float), np.asarray(flux, float)
    good = np.isfinite(t) & np.isfinite(y)
    t, y = t[good], y[good]
    p, epoch, duration = (float(candidate[k]) for k in ('period_days', 'epoch_bkjd', 'duration_days'))
    empty = {'depth': None, 'depth_error': None, 'snr': None, 'epoch_bkjd': None,
             'epoch_offset_days': None, 'epoch_error_days': None, 'status': 'insufficient_data'}
    if len(t) < 20 or duration <= 0 or p <= 0:
        return empty
    # A local quarter fit may refine the epoch, but cannot slide a box entirely
    # off the supplied transit and reclassify that event as baseline.
    grid = np.linspace(-duration / 2, duration / 2, 41)
    best = None
    for offset in grid:
        _, _, inside = transit_masks(t, p, epoch + offset, duration)
        out = ~inside
        if inside.sum() < 3 or out.sum() < 3:
            continue
        dep = float(np.mean(y[out]) - np.mean(y[inside]))
        err = float(np.sqrt(np.var(y[inside], ddof=1) / inside.sum() + np.var(y[out], ddof=1) / out.sum()))
        snr = dep / err if err > 0 else None
        if snr is not None and (best is None or snr > best['snr']):
            best = {'depth': dep, 'depth_error': err, 'snr': snr,
                    'epoch_bkjd': epoch + float(offset), 'epoch_offset_days': float(offset),
                    'epoch_error_days': max(float(grid[1] - grid[0]), duration / max(abs(snr), 1.0)),
                    'status': 'screened'}
    if best is None:
        return empty
    if best['depth'] <= 0 or best['snr'] < min_snr:
        best['status'] = 'manual_check_required'
        best['reason'] = 'Quarter transit fit is nonpositive or below the minimum depth SNR.'
    return best


def depth_measure(y: np.ndarray, mask: np.ndarray, baseline: np.ndarray) -> dict:
    """Estimate a mean depth with empirical standard error; require enough samples."""
    inside, outside = y[mask], y[baseline]
    if len(inside) < 3 or len(outside) < 3:
        return {'depth': None, 'error': None, 'snr': None}
    depth = float(np.mean(outside) - np.mean(inside))
    err = float(np.sqrt(np.var(inside, ddof=1) / len(inside) + np.var(outside, ddof=1) / len(outside)))
    return {'depth': depth, 'error': err, 'snr': depth / err if err > 0 else None}


def photometric_checks(frame: pd.DataFrame, candidate: dict) -> dict:
    """Compare alternating events and measure a phase-half eclipse independently."""
    p, e, d = (candidate[k] for k in ('period_days', 'epoch_bkjd', 'duration_days'))
    phase, event, inside = transit_masks(frame.time.to_numpy(), p, e, d)
    secondary = np.abs(np.abs(phase) - p / 2) < d / 2
    baseline = (np.abs(phase) > d) & (np.abs(np.abs(phase) - p / 2) > d)
    y = frame.flux.to_numpy()
    odd = depth_measure(y, inside & (event % 2 == 1), baseline)
    even = depth_measure(y, inside & (event % 2 == 0), baseline)
    sigma = None
    if odd['error'] is not None and even['error'] is not None:
        denom = np.hypot(odd['error'], even['error'])
        sigma = abs(odd['depth'] - even['depth']) / denom if denom > 0 else None
    event_checks = []
    for index in np.unique(event[inside]):
        local = (event == index) & (np.abs(phase) > d) & (np.abs(phase) < 3 * d)
        measured = depth_measure(y, inside & (event == index), local)
        event_checks.append({'event': int(index), **measured})
    weights = np.array([max(0., r['snr']) ** 2 for r in event_checks if r['snr'] is not None])
    fraction, remaining = None, None
    if len(weights) >= 3 and weights.sum() > 0:
        fraction = float(weights.max() / weights.sum())
        remaining = float(np.sqrt(max(0., weights.sum() - weights.max())))
    return {'odd': odd, 'even': even, 'odd_even_sigma': sigma,
            'secondary': depth_measure(y, secondary, baseline), 'events': event_checks,
            'max_single_event_fraction': fraction, 'remaining_event_snr': remaining}


def digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def serial(obj):
    if isinstance(obj, dict): return {str(k): serial(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, np.ndarray)): return [serial(v) for v in obj]
    if isinstance(obj, (np.integer,)): return int(obj)
    if isinstance(obj, (np.bool_,)): return bool(obj)
    if isinstance(obj, (float, np.floating)): return float(obj) if np.isfinite(obj) else None
    return obj


def split_quarters(raw, min_cadences=1000, min_quarters=4, inspected=False):
    """Select quarters from time/error/quality coverage only, before flux is searched."""
    valid = raw[np.isfinite(raw.time) & np.isfinite(raw.flux_err) & (raw.flux_err > 0)]
    counts = valid.groupby('quarter').size()
    quarters = sorted(int(q) for q in counts[counts >= min_cadences].index)
    if len(quarters) < min_quarters: raise ValueError('Insufficient usable quarters for untouched evaluation')
    n = max(1, int(np.ceil(len(quarters) * .25)))
    discovery, holdout = quarters[:-n], quarters[-n:]
    plan = {'discovery_quarters': discovery, 'holdout_quarters': holdout,
            'minimum_cadences_per_quarter': min_cadences, 'minimum_usable_quarters': min_quarters,
            'untouched': not inspected, 'coverage_counts': {str(q): int(counts[q]) for q in quarters}}
    return raw[raw.quarter.isin(discovery)].copy(), raw[raw.quarter.isin(holdout)].copy(), plan


def circle_overlap(z, radius):
    """Analytic area occulted on a unit uniform-brightness stellar disk."""
    z = np.asarray(z, float); k = float(radius)
    area = np.zeros_like(z)
    full = z <= abs(1-k); area[full] = np.pi * min(1, k)**2
    partial = (z > abs(1-k)) & (z < 1+k)
    x = z[partial]
    a = np.arccos(np.clip((x*x+1-k*k)/(2*x), -1, 1))
    b = k*k*np.arccos(np.clip((x*x+k*k-1)/(2*x*k), -1, 1))
    c = .5*np.sqrt(np.maximum(0, (-x+1+k)*(x+1-k)*(x-1+k)*(x+1+k)))
    area[partial] = a+b-c
    return area


def physical_transit(time, period, epoch, duration, depth, exposure=.02043365, supersamples=21, impact=.3):
    """Circular opaque sphere, uniform stellar disk (u=0), exposure-integrated.

    Duration is first-to-fourth contact. a/R is derived from duration and b;
    this is a specified physical geometry, not a host-density measurement.
    """
    if not (0 < depth < 1 and 0 < duration < period/2 and 0 <= impact < 1): raise ValueError('Invalid transit geometry')
    k = np.sqrt(depth)
    a = np.sqrt(impact**2 + ((1+k)**2-impact**2)/np.sin(np.pi*duration/period)**2)
    offsets = ((np.arange(supersamples)+.5)/supersamples-.5)*exposure
    t = np.asarray(time, float)[:, None]+offsets
    angle = 2*np.pi*(t-epoch)/period
    z = np.sqrt(a*a*np.sin(angle)**2+impact**2*np.cos(angle)**2)
    occult = circle_overlap(z, k)/np.pi
    occult[np.cos(angle) <= 0] = 0
    return 1-occult.mean(axis=1)


def diagnostics(frame, peak, settings):
    checks = photometric_checks(frame, peak)
    quarters = [{'quarter': int(q), **quarter_depth_epoch(g.time.to_numpy(), g.flux.to_numpy(), peak)} for q,g in frame.groupby('quarter')]
    rn = red_noise(frame, peak)
    missing = [k for k in ('odd_even_sigma', 'max_single_event_fraction', 'remaining_event_snr') if checks.get(k) is None]
    if checks['secondary'].get('snr') is None: missing.append('secondary')
    passed = bool(peak['initial_candidate'] and not missing and checks['odd_even_sigma'] <= settings['odd_even_sigma']
        and (checks['secondary']['snr'] or 0) < settings['secondary_sigma']
        and checks['max_single_event_fraction'] <= settings['max_single_event_fraction']
        and checks['remaining_event_snr'] >= settings['min_remaining_event_snr'])
    repair = settings.get('repair_screen', False)
    modulation = smooth_modulation(frame, peak) if repair else None
    beta = rn.get('ratio')
    corrected_snr = peak['snr'] / max(1., beta) if beta is not None else None
    measurable = sum(e.get('snr') is not None for e in checks['events'])
    reasons = []
    if repair:
        if corrected_snr is None or corrected_snr < 7: reasons.append('correlated_noise_snr')
        if modulation['status'] != 'screened': reasons.append('missing_modulation')
        elif modulation['flagged']: reasons.append('broad_modulation')
        if measurable < 3: reasons.append('sampled_event_coverage')
        passed = passed and not reasons
    return serial({'photometry': checks, 'quarter_consistency': quarters, 'red_noise': rn,
                   'modulation': modulation, 'correlated_noise_snr': corrected_snr,
                   'measurable_events': measurable, 'repair_rejections': reasons,
                   'missing': missing, 'photometry_vetting_pass': passed,
                   'localization': {'status': 'unresolved', 'reason': 'No validated PRF source assignment in photometry gate'}})


def fixed_holdout(frame, frozen, min_transits=3):
    """Measure frozen event template. No period/epoch/duration fitting is permitted."""
    p,e,d = (frozen[k] for k in ('period_days','epoch_bkjd','duration_days'))
    phase, event, inside = transit_masks(frame.time.to_numpy(), p,e,d)
    outside = (np.abs(phase)>d) & (np.abs(phase)<3*d)
    n = observed_transits(event, inside, 2)
    if inside.sum()<3 or outside.sum()<3: return {'status': 'insufficient_coverage','observed_transits': n,'ephemeris_hash': digest(frozen)}
    y = frame.flux.to_numpy(); err = frame.flux_err.to_numpy()
    def mean(mask):
        weights=1/err[mask]**2
        return float(np.sum(y[mask]*weights)/weights.sum()), float(1/np.sqrt(weights.sum()))
    im,ie=mean(inside); om,oe=mean(outside); depth=om-im; error=np.hypot(ie,oe)
    return serial({'status':'evaluated','depth':depth,'error':error,'snr':depth/error,'observed_transits':n,
        'passes_frozen_threshold':bool(n>=min_transits and depth>0 and depth/error>=7),
        'ephemeris_hash':digest(frozen),'refit':False})


def ephemeris_overlap(candidate, record, times, tolerance=.01):
    """Nearby harmonic/epoch screen with actual sampled event overlap."""
    p = float(record['period_days']); ratio = candidate['period_days']/p
    harmonic = next((r for r in (.25,.5,1.,2.,4.) if abs(ratio/r-1)<=tolerance), None)
    _,_,a = transit_masks(times,candidate['period_days'],candidate['epoch_bkjd'],candidate['duration_days'])
    _,_,b = transit_masks(times,p,float(record['epoch_bkjd']),float(record['duration_days']))
    return {'harmonic_ratio':harmonic, 'candidate_in_samples':int(a.sum()), 'overlap_samples':int((a&b).sum()),
            'event_overlap_fraction':float((a&b).sum()/a.sum()) if a.any() else None,
            'possible_ephemeris_contamination':bool(harmonic is not None and a.any() and (a&b).sum()/a.sum()>=.5)}
