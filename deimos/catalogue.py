"""Catalogue check of a finished batch's holdout-supported lead groups .

    python -m deimos.catalogue BATCH_OUTPUT_DIR [--out queue.csv] [--cache DIR]

Matches each supported signal family against KOI cumulative, TCE DR25/DR24, TOI, Kepler EB (J/AJ/151/68), VSX, Gaia DR3 EB
and Gaia DR3 neighbours, adds the Kepler DV source-position (dikco) verdict, and writes one queue row per family.
Catalogue responses are cached under --cache (raw/, query-log.jsonl), so a rerun is offline. Nothing here opens a holdout.
"""
import argparse
import csv
import datetime as dt
import json
import math
import os
import pathlib
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from fractions import Fraction

SUPPORTED = 'supported_fixed_photometry_host_unresolved'
FORBIDDEN = re.compile(r'sealed|access|split|holdout-\d', re.IGNORECASE)
RAT = sorted({Fraction(n, m) for n in range(1, 12) for m in range(1, 12)})
STRICT, LOOSE = 1e-4, 1e-3
NEXT_CHECKS = ('reconcile depth', 'check alias/odd-even', 'resolve the catalogue host', 'run the pixel host test',
               'request external data', 'check the shared systematic')
NEXT_CHECKS = ('reconcile depth', 'check alias/odd-even', 'resolve the catalogue host', 'run the pixel host test',
               'request external data', 'check the shared systematic')
# Systematic-period flag (warning only, never a drop). Thresholds fixed 2026-10-17 23:41Z, before any control-arm
# histogram or real-arm queue row was looked at for this purpose. Calibrated on the inverted arm only.
SYS_ARM = 'inverted'
SYS_ARM_FILE = SYS_ARM + '.json'
SYS_BIN_DEX = 0.01            # log10(period) bin width (2.3 %)
SYS_NEIGHBOURS = (2, 6)       # baseline = median count of the bins 2..6 away on each side (up to 10 bins)
SYS_FACTOR = 4.0              # flag a bin whose count >= 4 x max(baseline, 1) ...
SYS_MIN_COUNT = 8             # ... and >= 8 vetting-pass inverted peaks; a row is flagged within +-1 bin of it
SYS_EPH_DEV = 2e-3            # shared ephemeris: |P/P' - 1| <= 0.2 % ...
SYS_EPH_DAYS = 0.5            # ... and epochs within 0.5 d plus the drift |P/P' - 1| x |T - T'|
SYS_EPH_MIN_SOURCES = 5       # flag when >= 5 other sources carry such an inverted-arm vetting-pass peak


def utc():
    return dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')


def write_csv(path, rows, cols=None):
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = cols or list(dict.fromkeys(k for r in rows for k in r))
    tmp = path.with_name('.' + path.name + '.tmp')
    with open(tmp, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow({c: ('' if r.get(c) is None else json.dumps(r[c], sort_keys=True)
                           if c == 'event_timing' and isinstance(r[c], dict) else r.get(c)) for c in cols})
    os.replace(tmp, path)


def fnum(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def gaia_id(source):
    m = re.search(r'(\d{6,})', str(source))
    return m.group(1) if m else str(source)



def _load_output_json(folder, name):
    if FORBIDDEN.search(name):
        raise ValueError(f'refusing to open {name}')
    with open(pathlib.Path(folder) / name) as fh:
        return json.load(fh)


def peak_summary(p):
    """Flat, numbers-only view of one peak (its diagnostics carried into the queue)."""
    d = p.get('diagnostics') or {}
    ph = d.get('photometry') or {}
    sec = ph.get('secondary') or {}
    odd, even = ph.get('odd') or {}, ph.get('even') or {}
    qc = [q for q in (d.get('quarter_consistency') or []) if isinstance(q, dict)]
    ev = [e for e in (ph.get('events') or []) if isinstance(e, dict)]
    mod = d.get('modulation') or {}
    return {
        'event_timing': event_timing(p),
        'period_days': p.get('period_days'), 'epoch_bkjd': p.get('epoch_bkjd'), 'duration_days': p.get('duration_days'),
        'depth_ppm': None if fnum(p.get('depth')) is None else p['depth'] * 1e6, 'snr': p.get('snr'),
        'correlated_noise_snr': d.get('correlated_noise_snr'), 'n_transits': p.get('number_of_transits'),
        'initial_candidate': p.get('initial_candidate'), 'vetting_pass': d.get('photometry_vetting_pass'),
        'odd_depth_ppm': None if fnum(odd.get('depth')) is None else odd['depth'] * 1e6,
        'even_depth_ppm': None if fnum(even.get('depth')) is None else even['depth'] * 1e6,
        'odd_even_sigma': ph.get('odd_even_sigma'),
        'secondary_depth_ppm': None if fnum(sec.get('depth')) is None else sec['depth'] * 1e6,
        'secondary_snr': sec.get('snr'), 'max_single_event_fraction': ph.get('max_single_event_fraction'),
        'remaining_event_snr': ph.get('remaining_event_snr'),
        'events_total': len(ev), 'events_snr_ge3': sum(1 for e in ev if (fnum(e.get('snr')) or 0) >= 3),
        'events_negative': sum(1 for e in ev if (fnum(e.get('depth')) or 0) < 0),
        'quarters_tested': len(qc), 'quarters_snr_ge3': sum(1 for q in qc if (fnum(q.get('snr')) or 0) >= 3),
        # Depth range over the quarters measured at SNR >= 5 (noisy quarters would dominate a plain min/max).
        'quarter_depth_min_ppm': min((q['depth'] * 1e6 for q in qc if (fnum(q.get('snr')) or 0) >= 5), default=None),
        'quarter_depth_max_ppm': max((q['depth'] * 1e6 for q in qc if (fnum(q.get('snr')) or 0) >= 5), default=None),
        'red_noise_ratio': (d.get('red_noise') or {}).get('ratio'),
        'modulation_amp_to_depth': mod.get('amplitude_to_depth'),
        'folded_modulation_amp_to_depth': mod.get('folded_amplitude_to_depth'),
        'repair_rejections': ';'.join(str(x) for x in (d.get('repair_rejections') or [])) or None,
    }


# =============================================================================== catalogue match
def _viz(src, ra, dec, r):
    return 'https://vizier.cds.unistra.fr/viz-bin/asu-tsv?' + urllib.parse.urlencode(
        {'-source': src, '-c': f'{ra} {dec}', '-c.rs': r, '-out': '**', '-out.max': 500, '-out.add': '_r',
         '-sort': '_r'})


def _nea(q):
    return 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync?' + urllib.parse.urlencode({'query': q, 'format': 'csv'})


DIKCO = ','.join(f'{{p}}dikco_{c}{e}' for c in ('mra', 'mdec', 'msky') for e in ('', '_err'))


def cone_queries(ra, dec):
    """VizieR and archive cone-query endpoints, plus durations and the DV dikco columns."""
    d = 120 / 3600
    dra = d / math.cos(math.radians(dec))
    box = f'ra between {ra - dra} and {ra + dra} and dec between {dec - d} and {dec + d}'
    return {
        'koi-cumulative': _nea('select kepid,kepoi_name,kepler_name,koi_disposition,koi_fpflag_nt,koi_fpflag_ss,'
                               'koi_fpflag_co,koi_fpflag_ec,koi_period,koi_time0bk,koi_duration,koi_depth,koi_comment,'
                               f'ra,dec,{DIKCO.format(p="koi_")} from cumulative where {box}'),
        'tce-dr25': _nea(f'select kepid,tce_plnt_num,tce_period,tce_time0bk,tce_duration,tce_depth,ra,dec,'
                         f'{DIKCO.format(p="tce_")} from q1_q17_dr25_tce where {box}'),
        'tce-dr24': _nea(f'select kepid,tce_plnt_num,tce_period,tce_time0bk,tce_depth,ra,dec from q1_q16_tce where {box}'),
        'toi': _nea(f'select toi,tid,tfopwg_disp,pl_orbper,pl_tranmid,pl_trandep,ra,dec from toi where {box}'),
        'kepler-eb': _viz('J/AJ/151/68', ra, dec, 120),
        'vsx': _viz('B/vsx/vsx', ra, dec, 120),
        'gaia-vari-eb': _viz('I/358/veb', ra, dec, 60),
        'gaia-main': _viz('I/355/gaiadr3', ra, dec, 15),
    }


def fetch_catalogues(sources, raw_dir, log_path):
    """One cone set per Gaia source; a response already saved with HTTP 200 for the same URL is reused."""
    raw_dir.mkdir(parents=True, exist_ok=True)
    seen = {}
    if log_path.exists():
        for ln in open(log_path):
            r = json.loads(ln)
            if r.get('http') == 200:
                seen[(r['file'], r['url'])] = True
    jobs = [(g, k, url) for g, (ra, dec) in sorted(sources.items()) for k, url in cone_queries(ra, dec).items()]
    failures = []
    with open(log_path, 'a') as log:
        for g, k, url in jobs:
            f = raw_dir / f'{g}-{k}.txt'
            if f.exists() and seen.get((f.name, url)):
                continue
            body, code = b'', None
            for _attempt in range(3):
                try:
                    with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'stdlib-urllib'}),
                                                timeout=120) as r:
                        body, code = r.read(), r.status
                    break
                except urllib.error.HTTPError as e:
                    body, code = e.read(), e.code
                    break
                except Exception as e:      # network hiccup: retry, then record the failure
                    body, code = repr(e).encode(), None
                    time.sleep(3)
            f.write_bytes(body)
            log.write(json.dumps({'source': g, 'query': k, 'utc': utc(), 'url': url, 'http': code,
                                  'bytes': len(body), 'file': f.name}) + '\n')
            log.flush()
            if code != 200:
                failures.append((g, k, code))
    return failures


def _tsv(path):
    if not path.exists():
        return []
    lines = [ln for ln in open(path) if not ln.startswith('#') and ln.strip()]
    if len(lines) < 4:
        return []
    head = [c.strip() for c in lines[0].split('\t')]
    return [dict(zip(head, [c.strip() for c in ln.rstrip('\n').split('\t')])) for ln in lines[3:]]


def _csv(path):
    try:
        with open(path) as fh:
            return list(csv.DictReader(fh))
    except (OSError, csv.Error):
        return []


def relation(p, p_ref):
    r = p / p_ref
    best = min(RAT, key=lambda q: abs(r / q - 1))
    return best, r / float(best) - 1


def phase_residual(t, p_ref, t_ref, q):
    """Days from our epoch to the nearest coincidence time of the two ephemerides (unit P_ref / q.denominator)."""
    if None in (t, p_ref, t_ref):
        return None
    unit = p_ref / q.denominator
    x = (t - t_ref) / unit
    return abs(x - round(x)) * unit


def phase_tol(duration, dev):
    """Proposal tolerance only. Drift must never enlarge the identity window."""
    return 0.75 * (duration or 0)


def event_timing(peak):
    """Saved discovery measurements only; no flux read, fit or holdout access."""
    if peak.get('event_timing'):
        value = peak['event_timing']
        return json.loads(value) if isinstance(value, str) else value
    p, t = fnum(peak.get('period_days')), fnum(peak.get('epoch_bkjd'))
    events = ((peak.get('diagnostics') or {}).get('photometry') or {}).get('events') or []
    measured = [e for e in events if fnum(e.get('snr')) is not None]
    times = [t + e['event'] * p for e in measured] if p and t is not None else []
    return {'positive_indices': [e['event'] for e in measured if e['snr'] >= 3 and (e.get('depth') or 0) > 0],
            'start': min(times) if times else None, 'end': max(times) if times else None,
            'measured_events': len(measured), 'evidence': 'saved discovery diagnostics.photometry.events'}


def hydrate_timing(items, path):
    cache = json.loads(path.read_text()) if path.exists() else {}
    for m in items:
        key = f"{m['rank']}:{m['peak_index']}"
        if key in cache:
            m['peak']['event_timing'] = cache[key]


def signal_explained_by(child, parent, catalogue=False):
    """Directional identity proof. Ratios propose; measured dip windows establish coverage.

    Every significant positive child event must overlap a measured positive parent event
    (or a catalogue primary eclipse). A measured parent secondary can explain phase 0.5.
    No accumulated drift is added to the duration window. Missing timing fails open to review.
    """
    cp, ct, cd = (fnum(child.get(k)) for k in ('period_days', 'epoch_bkjd', 'duration_days'))
    pp, pt, pd = (fnum(parent.get(k)) for k in ('period_days', 'epoch_bkjd', 'duration_days'))
    result = {'demonstrated': False, 'reason': 'missing epoch, duration or event evidence'}
    if None in (cp, ct, cd, pp, pt, pd) or min(cp, pp, cd, pd) <= 0:
        return result
    q, dev = relation(cp, pp)
    result.update(relation=str(q), rel_dev=dev)
    if abs(dev) > LOOSE:
        return dict(result, reason='different periods')
    ce, pe = event_timing(child), event_timing(parent)
    times = [ct + i * cp for i in ce['positive_indices']]
    if len(times) < 2 or ce.get('start') is None or ce.get('end') is None:
        return result
    tol = (cd + pd) / 2
    parent_times = sorted(pt + i * pp for i in pe['positive_indices'])
    if not catalogue and len(parent_times) < 2:
        return result
    # A catalogue is an external ephemeris, not evidence for a secondary by itself.
    secondary = not catalogue and (fnum(parent.get('secondary_snr')) or 0) >= 3
    residuals = []
    import bisect
    for t in times:
        if catalogue:
            distance = phase_residual(t, pp, pt, Fraction(1))
        else:
            i = bisect.bisect_left(parent_times, t)
            distance = min((abs(t - parent_times[j]) for j in (i-1, i) if 0 <= j < len(parent_times)), default=math.inf)
            if secondary:
                distance = min(distance, phase_residual(t, pp, pt + pp / 2, Fraction(1)))
        residuals.append(distance)
    # Explicitly constrain the drift over the observed baseline, even for sparse events.
    baseline = ce['end'] - ce['start']
    drift = abs(dev) * baseline
    result.update(events=len(times), explained_events=sum(d <= tol for d in residuals),
                  max_residual_days=max(residuals), tolerance_days=tol,
                  baseline_days=baseline, accumulated_drift_days=drift,
                  evidence=ce.get('evidence'), parent_evidence=pe.get('evidence') if not catalogue else 'catalogue ephemeris')
    result.update(evidence_sha256=ce.get('sha256'), parent_evidence_sha256=pe.get('sha256'))
    ok = all(d <= tol for d in residuals) and drift <= tol
    return dict(result, demonstrated=ok, reason='measured dips covered across baseline' if ok else 'unexplained events or accumulated drift exceeds duration')


def same_ephemeris(a, b, baseline=1400.):
    """Stable identity across changing queue IDs; no rational-ratio inheritance."""
    ap, at, ad = (fnum(a.get(k)) for k in ('period_days', 'epoch_bkjd', 'duration_days'))
    bp, bt = fnum(b.get('period_days')), fnum(b.get('epoch_bkjd'))
    if None in (ap, at, ad, bp, bt) or min(ap, bp, ad) <= 0:
        return False
    if abs(ap / bp - 1) > LOOSE:
        return False
    tol = ad / 2
    if abs(ap / bp - 1) * baseline > tol:
        return False
    # Test both ends, rather than allowing drift to increase the tolerance.
    return all(phase_residual(at + i * ap, bp, bt, Fraction(1)) <= tol
               for i in (0, math.ceil(baseline / ap)))


def dikco_verdict(row, prefix, tra, tdec):
    """Kepler DV difference-image source (catalogued KIC position + dikco offset) against our target.

    Rule of b3b-vetting-20261017T060456Z/closeout/koi_centroid.py, plus the same-star case (the catalogued star
    is our target): on target if the source lies within max(3 sigma, 1") of it."""
    kra, kdec = fnum(row.get('ra')), fnum(row.get('dec'))
    mra, mdec = fnum(row.get(prefix + 'dikco_mra')), fnum(row.get(prefix + 'dikco_mdec'))
    err = fnum(row.get(prefix + 'dikco_msky_err'))
    if None in (kra, kdec):
        return {'dikco_verdict': 'no catalogue position'}
    dx = (tra - kra) * math.cos(math.radians(kdec)) * 3600
    dy = (tdec - kdec) * 3600
    sep = math.hypot(dx, dy)
    if mra is None or mdec is None:
        return {'dikco_verdict': 'no dikco offset', 'cat_star_sep_arcsec': round(sep, 2)}
    to_tgt, to_kic = math.hypot(mra - dx, mdec - dy), math.hypot(mra, mdec)
    if sep < 1.0:
        v = 'source at target' if err and to_tgt < max(3 * err, 1.0) else 'source offset from target'
    elif err and to_tgt < max(3 * err, 1.0) and to_kic > 3 * err:
        v = 'source at target'
    elif err and to_kic < max(3 * err, 1.0) and to_tgt > 3 * err:
        v = 'source at catalogued star'
    else:
        v = 'ambiguous'
    return {'dikco_verdict': v, 'cat_star_sep_arcsec': round(sep, 2), 'dikco_src_to_target_arcsec': round(to_tgt, 2),
            'dikco_src_to_catstar_arcsec': round(to_kic, 2), 'dikco_err_arcsec': err}


def catalogue_refs(raw_dir, g):
    """Every catalogued period near one source (order: KOI, TCE DR25, TCE DR24, Kepler EB, VSX, Gaia EB, TOI)."""
    out = []
    for r in _csv(raw_dir / f'{g}-koi-cumulative.txt'):
        fp = ''.join(r.get(c, '') for c in ('koi_fpflag_nt', 'koi_fpflag_ss', 'koi_fpflag_co', 'koi_fpflag_ec'))
        out.append({'cat': 'KOI', 'ref': f"{r['kepoi_name']} {r.get('kepler_name', '')} {r['koi_disposition']} fp={fp} "
                    f"kic{r['kepid']}", 'P_ref': fnum(r.get('koi_period')), 'T_ref': fnum(r.get('koi_time0bk')),
                    'D_ref': (fnum(r.get('koi_duration')) or 0) / 24 or None,
                    'disposition': r.get('koi_disposition'), 'row': r, 'prefix': 'koi_', 'kepid': r.get('kepid'),
                    'name': r['kepoi_name']})
    for t in ('tce-dr25', 'tce-dr24'):
        for r in _csv(raw_dir / f'{g}-{t}.txt'):
            out.append({'cat': t.upper(), 'ref': f"kic{r['kepid']}.{r['tce_plnt_num']}", 'P_ref': fnum(r.get('tce_period')),
                        'T_ref': fnum(r.get('tce_time0bk')), 'row': r, 'prefix': 'tce_' if t == 'tce-dr25' else None,
                        'D_ref': (fnum(r.get('tce_duration')) or 0) / 24 or None,
                        'kepid': r.get('kepid'), 'name': f"{r['kepid']}.{r['tce_plnt_num']}"})
    for r in _tsv(raw_dir / f'{g}-kepler-eb.txt'):
        bjd0 = fnum(r.get('BJD0'))      # BJD - 2400000; BKJD = BJD - 2454833
        out.append({'cat': 'KeplerEB', 'ref': f"KIC{r['KIC']} r={r['_r']}\" morph={r.get('Morph', '')}",
                    'P_ref': fnum(r.get('Per')), 'T_ref': None if bjd0 is None else bjd0 - 54833.0})
    for r in _tsv(raw_dir / f'{g}-vsx.txt'):
        if fnum(r.get('Period')):
            out.append({'cat': 'VSX', 'ref': f"{r['Name']} {r['Type']} r={r['_r']}\"", 'P_ref': fnum(r['Period']),
                        'T_ref': None})
    for r in _tsv(raw_dir / f'{g}-gaia-vari-eb.txt'):
        f = fnum(r.get('Freq'))
        if f:
            out.append({'cat': 'GaiaEB', 'ref': f"{r['Source']} r={r['_r']}\"", 'P_ref': 1 / f, 'T_ref': None})
    for r in _csv(raw_dir / f'{g}-toi.txt'):
        if fnum(r.get('pl_orbper')):
            out.append({'cat': 'TOI', 'ref': f"TOI {r['toi']} TIC{r['tid']} {r.get('tfopwg_disp', '')}",
                        'P_ref': fnum(r['pl_orbper']), 'T_ref': None})
    return [x for x in out if x['P_ref']]


def gaia_neighbours(raw_dir, g):
    rows = _tsv(raw_dir / f'{g}-gaia-main.txt')
    others = [r for r in rows if r.get('Source') != g and fnum(r.get('_r')) is not None]
    near = min(others, key=lambda r: fnum(r['_r']), default=None)
    bright = min((r for r in others if fnum(r.get('Gmag')) is not None), key=lambda r: fnum(r['Gmag']), default=None)
    tgt = next((r for r in rows if r.get('Source') == g), None)
    # VizieR _r is in arcsec for these cones (units row of the raw responses).
    return {'gaia_rows_15arcsec': len(rows), 'target_G': fnum(tgt.get('Gmag')) if tgt else None,
            'target_RUWE': fnum(tgt.get('RUWE')) if tgt else None,
            'nearest_neighbour': f"{fnum(near['_r']):.2f}\" G{near.get('Gmag', '')}" if near else '',
            'brightest_neighbour': f"G{bright.get('Gmag')} at {fnum(bright['_r']):.2f}\"" if bright else ''}


def match_member(m, refs, tra, tdec):
    """Best-relation classification, plus phase alignment and the dikco verdict for KOI/TCE refs."""
    pk = m['peak']
    P_, T_ = fnum(pk['period_days']), fnum(pk['epoch_bkjd'])
    best, aligned = None, []
    for ref in refs:
        q, dev = relation(P_, ref['P_ref'])
        cand = (ref, q, dev)
        if best is None or abs(dev) < abs(best[2]):
            best = cand
        if abs(dev) <= LOOSE:
            aligned.append(cand)
    res = {'match': 'no_catalogue_period'}
    if best:
        ref, q, dev = best
        ph = phase_residual(T_, ref['P_ref'], ref.get('T_ref'), q)
        res = {'match': 'strict' if abs(dev) <= STRICT else 'loose' if abs(dev) <= LOOSE else 'none',
               'cat': ref['cat'], 'ref': ref['ref'], 'P_ref': ref['P_ref'], 'relation': str(q), 'rel_dev': f'{dev:.1e}',
               'phase_ok': None if ph is None else ph <= phase_tol(fnum(pk['duration_days']), dev)}
    detail = []
    for ref, q, dev in aligned:
        d = {'cat': ref['cat'], 'name': ref.get('name') or ref['ref'], 'relation': str(q), 'rel_dev': dev,
             'class': 'strict' if abs(dev) <= STRICT else 'loose', 'disposition': ref.get('disposition'),
             'period_days': ref['P_ref'], 'epoch_bkjd': ref.get('T_ref'), 'duration_days': ref.get('D_ref')}
        ph = phase_residual(T_, ref['P_ref'], ref.get('T_ref'), q)
        d['phase_residual_days'] = None if ph is None else round(ph, 4)
        d['phase_offset_cycles'] = None if ph is None else round(ph / (ref['P_ref'] / q.denominator), 3)
        d['phase_ok'] = None if ph is None else ph <= phase_tol(fnum(pk['duration_days']), dev)
        # DV locates the catalogue eclipse, which must first explain this signal's dips.
        d['timing_proof'] = signal_explained_by(pk, dict(period_days=ref['P_ref'], epoch_bkjd=ref.get('T_ref'),
                                                        duration_days=ref.get('D_ref')), catalogue=True)
        if d['phase_ok'] is False and d['class'] == 'loose':   # loose and epochs off: treated as chance
            continue
        if ref.get('prefix'):
            d.update(dikco_verdict(ref['row'], ref['prefix'], tra, tdec))
        detail.append(d)
    # Strict refs (any epoch; an offset near 0.5 cycle can be the other eclipse of an EB) and loose refs whose
    # epoch aligns or is unknown, closest period first.
    res['aligned_refs'] = sorted(detail, key=lambda d: abs(d['rel_dev']))
    return res



# =============================================================================== families and next check
def families(groups):
    """Join only aliases whose observed dips are explained by the parent's measured events."""
    fams = []
    for g in groups:
        rep = next(m for m in g['members'] if m['holdout']['status'] == SUPPORTED)
        P_ = rep['peak']['period_days']
        joined = None
        for f in fams:
            if gaia_id(f['source']) != gaia_id(rep['source']):
                continue
            par = f['parent_member']['peak']
            q, dev = relation(P_, par['period_days'])
            proof = signal_explained_by(rep['peak'], par)
            if proof['demonstrated']:
                joined = f
                f['aliases'].append({'lead_group': g['lead_group'], 'relation': str(q), 'rel_dev': dev, 'member': rep,
                                     'timing_proof': proof})
                break
        if joined is None:
            fams.append({'rank': rep['rank'], 'source': rep['source'], 'dependency_group': g['dependency_group'],
                         'parent_lead_group': g['lead_group'], 'parent_member': rep, 'aliases': [],
                         'parent_group_members': len(g['members'])})
    return fams


def proposed_check(row):
    """One suggested follow-up check from a fixed menu; a proposal only."""
    cm, dv = row.get('catalogue_match'), row.get('dikco_verdict') or ''
    if cm in ('strict', 'loose') and dv != 'source at target':
        return 'resolve the catalogue host'
    if cm == 'strict_phase_offset' or row.get('related_families'):
        return 'check alias/odd-even'      # same period at another epoch: possible other eclipse of an EB
    oe, sec = fnum(row.get('odd_even_sigma')) or 0, abs(fnum(row.get('secondary_snr')) or 0)
    rel = [x.split('=', 1)[1] for x in str(row.get('alias_relations') or '').split(';') if '=' in x]
    harmonic = all('/' not in q or q.startswith('1/') for q in rel)
    if not harmonic or oe >= 2.5 or sec >= 3:
        return 'check alias/odd-even'
    hr = fnum(row.get('holdout_depth_ratio'))
    qmin, qmax = fnum(row.get('quarter_depth_min_ppm')), fnum(row.get('quarter_depth_max_ppm'))
    if (hr is not None and not 0.67 <= hr <= 1.5) or (qmin is not None and qmax and qmin < 0.5 * qmax):
        return 'reconcile depth'
    return 'run the pixel host test'


# =============================================================================== systematic-period flag
def systematic_bands(peaks):
    """Over-dense log-period bins of the inverted arm's vetting-pass peaks. No real-arm value enters the rule."""
    import statistics
    counts = {}
    for p in peaks:
        P_ = fnum(p.get('period_days'))
        if p.get('arm') != SYS_ARM or str(p.get('vetting_pass')) != 'True' or not P_ or P_ <= 0:
            continue
        b = math.floor(math.log10(P_) / SYS_BIN_DEX)
        counts[b] = counts.get(b, 0) + 1
    flagged = {}
    if not counts:
        return flagged
    lo, hi = min(counts), max(counts)
    for b in range(lo, hi + 1):
        c = counts.get(b, 0)
        nb = [counts.get(b + sgn * k, 0) for k in range(SYS_NEIGHBOURS[0], SYS_NEIGHBOURS[1] + 1) for sgn in (-1, 1)
              if lo <= b + sgn * k <= hi]
        base = statistics.median(nb) if nb else 0
        if c >= SYS_MIN_COUNT and c >= SYS_FACTOR * max(base, 1):
            flagged[b] = {'bin': b, 'count': c, 'baseline': base, 'P_lo': round(10 ** (b * SYS_BIN_DEX), 3),
                          'P_hi': round(10 ** ((b + 1) * SYS_BIN_DEX), 3)}
    return flagged


def ephemeris_index(peaks):
    return sorted((fnum(p['period_days']), fnum(p['epoch_bkjd']), str(p['rank'])) for p in peaks
                  if p.get('arm') == SYS_ARM and str(p.get('vetting_pass')) == 'True'
                  and fnum(p.get('period_days')) and fnum(p.get('epoch_bkjd')) is not None)


def shared_ephemeris_sources(P_, T_, own_rank, index):
    """Other sources whose inverted-arm vetting-pass peaks share this period (0.2 %) and epoch (0.5 d + drift)."""
    import bisect
    i = bisect.bisect_left(index, (P_ * (1 - SYS_EPH_DEV), -math.inf, ''))
    j = bisect.bisect_right(index, (P_ * (1 + SYS_EPH_DEV), math.inf, '~'))
    out = set()
    for Pc, Tc, rk in index[i:j]:
        if rk == str(own_rank):
            continue
        dev = abs(P_ / Pc - 1)
        if phase_residual(T_, Pc, Tc, Fraction(1)) <= SYS_EPH_DAYS + dev * abs(T_ - Tc):
            out.add(rk)
    return out


def systematic_flag(row, bands, index):
    """Warning text for a queue row, or '' (never a drop)."""
    P_, T_ = fnum(row.get('period_days')), fnum(row.get('epoch_bkjd'))
    if not P_ or P_ <= 0:
        return ''
    why = []
    b = math.floor(math.log10(P_) / SYS_BIN_DEX)
    near = [bands[k] for k in (b - 1, b, b + 1) if k in bands]
    if near:
        top = max(near, key=lambda x: x['count'])
        why.append(f"period band {top['P_lo']}-{top['P_hi']} d: {top['count']} inverted-arm vetted peaks vs "
                   f"neighbour baseline {top['baseline']:g}")
    if T_ is not None:
        n = len(shared_ephemeris_sources(P_, T_, row.get('rank'), index))
        if n >= SYS_EPH_MIN_SOURCES:
            why.append(f'shared ephemeris with {n} other sources in the inverted arm')
    return '; '.join(why)


# =============================================================================== queue
QUEUE_COLS = [
    'section', 'batch', 'family_id', 'peak_index', 'peak_holdout_evaluation',
    'rank', 'gaia_dr3', 'kic', 'origin_kic', 'dependency_group', 'parent_lead_group',
    'alias_lead_groups', 'alias_relations', 'event_timing', 'period_days', 'epoch_bkjd', 'duration_days', 'depth_ppm', 'snr',
    'correlated_noise_snr', 'n_transits', 'odd_even_sigma', 'odd_depth_ppm', 'even_depth_ppm', 'secondary_depth_ppm',
    'secondary_snr', 'max_single_event_fraction', 'remaining_event_snr', 'events_total', 'events_snr_ge3',
    'events_negative', 'quarters_tested', 'quarters_snr_ge3', 'quarter_depth_min_ppm', 'quarter_depth_max_ppm',
    'red_noise_ratio', 'modulation_amp_to_depth', 'folded_modulation_amp_to_depth', 'holdout_status',
    'holdout_depth_ppm', 'holdout_error_ppm', 'holdout_depth_ratio', 'holdout_snr', 'holdout_transits',
    'holdout_odd_even_sigma', 'holdout_secondary_snr', 'catalogue_match', 'catalogue_member_lg', 'catalogue_cat',
    'catalogue_ref', 'catalogue_disposition', 'catalogue_relation', 'catalogue_rel_dev', 'catalogue_phase_ok',
    'member_matches', 'related_families',
    'catalogue_evidence', 'dikco_phase_ok', 'dikco_relation', 'dikco_rel_dev', 'dikco_ref', 'dikco_verdict', 'dikco_src_to_target_arcsec', 'dikco_err_arcsec', 'cat_star_sep_arcsec', 'aligned_refs',
    'target_G', 'target_RUWE', 'nearest_neighbour', 'brightest_neighbour', 'same_period_other_sources',
    'systematic_flag',
    'proposed_next_check']


def build_rows(batch, groups, src_info, raw_dir, timing_path):
    hydrate_timing([m for g in groups for m in g['members']], timing_path)
    rows = []
    by_rank = {}
    for g in groups:
        by_rank.setdefault(next(m for m in g['members'] if m['holdout']['status'] == SUPPORTED)['rank'], []).append(g)
    fid = 0
    for rank in sorted(by_rank, key=lambda r: min(g['lead_group'] for g in by_rank[r])):
        info = src_info[rank]
        g_id = gaia_id(info['source'])
        tra, tdec = float(info['ra']), float(info['dec'])
        refs = catalogue_refs(raw_dir, g_id)
        nb = gaia_neighbours(raw_dir, g_id)
        for fam in families(sorted(by_rank[rank], key=lambda g: g['lead_group'])):
            fid += 1
            par = fam['parent_member']
            pk, ho = par['peak'], par['holdout']
            row = {'section': 'supported', 'batch': batch, 'family_id': f'{batch}:F{fid:03d}', 'rank': rank,
                   'peak_index': par['peak_index'], 'peak_holdout_evaluation': ho.get('status', 'unknown'),
                   'gaia_dr3': g_id, 'kic': info.get('kic'), 'origin_kic': info.get('origin'),
                   'dependency_group': fam['dependency_group'], 'parent_lead_group': fam['parent_lead_group'],
                   'alias_lead_groups': ';'.join(str(a['lead_group']) for a in fam['aliases']),
                   'alias_relations': ';'.join(f"lg{a['lead_group']}={a['relation']}" for a in fam['aliases'])}
            row.update({k: pk.get(k) for k in QUEUE_COLS if k in pk})
            hd = ho.get('diagnostics') or {}
            row.update({'holdout_status': ho.get('status'), 'holdout_depth_ppm': ho.get('depth_ppm'),
                        'holdout_error_ppm': ho.get('error_ppm'), 'holdout_snr': ho.get('snr'),
                        'holdout_transits': ho.get('observed_transits'),
                        'holdout_depth_ratio': None if not (ho.get('depth_ppm') and pk.get('depth_ppm'))
                        else round(ho['depth_ppm'] / pk['depth_ppm'], 3),
                        'holdout_odd_even_sigma': hd.get('odd_even_sigma'), 'holdout_secondary_snr': hd.get('secondary_snr')})
            # Catalogue. Per member, the best |dev| class over all refs (period only) is kept. The family summary needs more: a strict ref counts if it is 1:1
            # or its epoch aligns; a loose ref only if its epoch aligns; a strict n/m relation with no catalogued
            # epoch (VSX, Gaia EB, TOI) is only 'commensurate'.
            cands = [(fam['parent_lead_group'], match_member(par, refs, tra, tdec))]
            cands += [(a['lead_group'], match_member(a['member'], refs, tra, tdec)) for a in fam['aliases']]
            row['member_matches'] = ';'.join(f"lg{_lg}:{c['match']}" + (f" {c['cat']} {c['relation']}" if c.get('cat') else '')
                                               for _lg, c in cands)

            def strength(d):
                if d['class'] == 'strict' and (d['phase_ok'] or (d['phase_ok'] is None and d['relation'] == '1')):
                    return 0, 'strict'
                if d['class'] == 'loose' and d['phase_ok']:
                    return 1, 'loose'
                if d['class'] == 'strict' and d['phase_ok'] is False:
                    return 2, 'strict_phase_offset'
                return 3, 'commensurate'      # n/m relation to a catalogue period without an epoch
            ranked = sorted(((strength(d), _lg, d) for _lg, c in cands for d in c['aligned_refs']),
                            key=lambda t: (t[0][0], t[1] != fam['parent_lead_group'], abs(t[2]['rel_dev'])))
            if ranked:
                (_k, label), lg, d = ranked[0]
                row.update({'catalogue_match': label, 'catalogue_member_lg': lg, 'catalogue_cat': d['cat'],
                            'catalogue_ref': d['name'], 'catalogue_relation': d['relation'],
                            'catalogue_rel_dev': f"{d['rel_dev']:.1e}", 'catalogue_phase_ok': d['phase_ok'],
                            'catalogue_disposition': d.get('disposition')})
            else:
                row['catalogue_match'] = 'none' if any(c['match'] != 'no_catalogue_period' for _lg, c in cands) else 'no_catalogue_period'
            # Kepler DV source position: from the best phase-aligned KOI/TCE, else the best-ranked one.
            kt = [(st, _lg, d) for st, _lg, d in ranked if d.get('dikco_verdict')
                  and _lg == fam['parent_lead_group']]
            pick = next((d for st, _lg, d in kt if d['timing_proof']['demonstrated']
                         and d.get('dikco_verdict') in ('source at target', 'source at catalogued star')
                         and (fnum(d.get('dikco_err_arcsec')) or 0) > 0), None)
            pick = pick or next((d for st, _lg, d in kt if d['phase_ok']), None) or (kt[0][2] if kt else None)
            row['catalogue_evidence'] = json.dumps([dict(lead_group=lg, **d) for lg, c in cands for d in c['aligned_refs']], sort_keys=True)
            if pick:
                row.update({'dikco_phase_ok': pick['phase_ok'], 'dikco_relation': pick['relation'],
                            'dikco_rel_dev': pick['rel_dev']})
                row.update({k: pick.get(k) for k in ('dikco_verdict', 'dikco_src_to_target_arcsec', 'dikco_err_arcsec',
                                                     'cat_star_sep_arcsec')})
                row['dikco_ref'] = f"{pick['cat']} {pick['name']} phase={'ok' if pick['phase_ok'] else 'off'}"
            row['aligned_refs'] = ' | '.join(
                f"lg{_lg}:{d['class']} {d['cat']} {d['name']} {d['relation']} {d['rel_dev']:.1e}"
                + (f" {d.get('disposition')}" if d.get('disposition') else '')
                + (f" phase={'ok' if d['phase_ok'] else 'off'}({d['phase_offset_cycles']}cyc)" if d['phase_ok'] is not None else '')
                + (f" dikco={d['dikco_verdict']}" if d.get('dikco_verdict') else '')
                for _lg, c in cands for d in c['aligned_refs'])
            row.update(nb)
            rows.append(row)
    # Same period on other sources (possible shared signal; not merged).
    for r in rows:
        others = [f"{o['rank']}:lg{o['parent_lead_group']}" for o in rows if o['rank'] != r['rank']
                  and abs(o['period_days'] / r['period_days'] - 1) <= STRICT]
        r['same_period_other_sources'] = ';'.join(others)
    for r in rows:
        rel = []
        for o in rows:
            if o is r or o['rank'] != r['rank']:
                continue
            q, dev = relation(r['period_days'], o['period_days'])
            if abs(dev) <= LOOSE and q.numerator <= 3 and q.denominator <= 3:
                ph = phase_residual(r['epoch_bkjd'], o['period_days'], o['epoch_bkjd'], q)
                rel.append(f"{o['family_id']}({q}, {ph}d; proposed only)")
        r['related_families'] = ';'.join(rel)
    for r in rows:
        r['proposed_next_check'] = proposed_check(r)
    return rows


# =============================================================================== batch output -> queue
def supported_groups(out):
    """Holdout-supported lead groups from reports/ranked-lead-groups.jsonl (no sealed data is read)."""
    path = pathlib.Path(out) / 'reports/ranked-lead-groups.jsonl'
    groups = []
    for ln in path.read_text().splitlines():
        if not ln.strip():
            continue
        g = json.loads(ln)
        if not any(m['holdout']['status'] == SUPPORTED for m in g['members']):
            continue
        members = []
        for m in g['members']:
            h = m.get('holdout') or {}
            ev = h.get('evaluation') or {}
            members.append({'rank': m['rank'], 'source': m['source'], 'group': m['group'],
                            'peak_index': m['peak_index'], 'peak': peak_summary(m['peak']),
                            'holdout': {'status': h.get('status'), 'depth_ppm': None if fnum(ev.get('depth')) is None
                                        else ev['depth'] * 1e6, 'snr': ev.get('snr'), 'error_ppm':
                                        None if fnum(ev.get('error')) is None else ev['error'] * 1e6,
                                        'observed_transits': ev.get('observed_transits'),
                                        'predicted_discovery_depth_snr': h.get('predicted_discovery_depth_snr'),
                                        'refit': h.get('refit'),
                                        'diagnostics': peak_summary({'diagnostics': h.get('diagnostics') or {}})
                                        if h.get('diagnostics') else None}})
        groups.append({'lead_group': g['lead_group'], 'dependency_group': g['dependency_group'], 'members': members})
    return groups


def source_info(out):
    return {int(r['rank']): r for r in (_load_output_json(f, 'receipt.json') for f in sorted(
        (pathlib.Path(out) / 'sources').iterdir()) if f.is_dir() and f.name.isdigit()) if r.get('status') == 'complete'}


def inverted_peaks(out, cache):
    """Vetting-pass inputs of the systematic flag; cached in peaks.csv so reruns need no arm files."""
    path = cache / 'peaks.csv'
    if not path.exists():
        rows = []
        for f in sorted((pathlib.Path(out) / 'sources').glob('[0-9]*/')):
            if (f / SYS_ARM_FILE).exists():
                rc = _load_output_json(f, 'receipt.json')
                for ix, p in enumerate(_load_output_json(f, SYS_ARM_FILE).get('peaks') or []):
                    s = peak_summary(p)
                    rows.append({'rank': rc['rank'], 'arm': SYS_ARM, 'peak_index': ix, **{k: s[k] for k in (
                        'period_days', 'epoch_bkjd', 'vetting_pass')}})
        write_csv(path, rows, ['rank', 'arm', 'peak_index', 'period_days', 'epoch_bkjd', 'vetting_pass'])
    return _csv(path)


def run(out, queue, cache):
    out, cache = pathlib.Path(out), pathlib.Path(cache)
    batch = out.name.removeprefix('out-')
    groups, info = supported_groups(out), source_info(out)
    sources = {}
    for g in groups:
        rep = next(m for m in g['members'] if m['holdout']['status'] == SUPPORTED)
        sources[gaia_id(info[rep['rank']]['source'])] = (float(info[rep['rank']]['ra']), float(info[rep['rank']]['dec']))
    fails = fetch_catalogues(sources, cache / 'raw', cache / 'query-log.jsonl')
    rows = build_rows(batch, groups, info, cache / 'raw', cache / 'event-timing.json')
    peaks = inverted_peaks(out, cache)
    bands, index = systematic_bands(peaks), ephemeris_index(peaks)
    for r in rows:
        r['systematic_flag'] = systematic_flag(r, bands, index)
        if r['systematic_flag']:
            r['proposed_next_check'] = 'check the shared systematic'
    write_csv(queue, rows, QUEUE_COLS)
    return rows, fails


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m deimos.catalogue', description=__doc__.split('\n')[0])
    ap.add_argument('batch_dir', help='finished batch output dir (reports/, sources/)')
    ap.add_argument('--out', default='queue.csv')
    ap.add_argument('--cache', default='catalogue-cache', help='catalogue responses (raw/, query-log.jsonl); reused offline')
    a = ap.parse_args(argv)
    rows, fails = run(a.batch_dir, a.out, a.cache)
    print(f'{len(rows)} families -> {a.out}; {len(fails)} failed catalogue queries')
    return 1 if fails else 0


if __name__ == '__main__':
    raise SystemExit(main())
