"""Incremental reporting.

Each finished source writes its own rows (sources/NNNN/arm-rows.jsonl, peak-rows.jsonl).
Aggregation appends those files to reports/*.jsonl exactly once per source; reports/index.json
records which ranks are in and the committed byte length of each aggregate, so a crash between
append and index update is repaired by truncating back. Nothing re-reads the whole batch except
`report`, which also builds lead groups and CSVs.
"""
import os, json, pathlib, time
import pandas as pd
from .util import save, save_bytes, load, flock, source_key

AGGREGATES = ['attempts.jsonl', 'arm-results.jsonl', 'peaks.jsonl']


def _dumps(row):
    return json.dumps(row, separators=(',', ':'), allow_nan=False, default=str)


def arm_rows(folder, receipt):
    rows = []
    base = dict(rank=receipt['rank'], source=receipt['source'], group=receipt['group'])
    for arm, info in receipt.get('arms', {}).items():
        row = dict(base, arm=arm, source_status=receipt['status'])
        row.update(evidence_status='searched', peaks=info['peaks'], survivors=info['survivors'], seconds=info['seconds'], arm_sha256=info['sha256'])
        row['injection_detected'] = receipt.get('injection_detected') if arm == 'injected' else None
        row['injection_recovered'] = receipt.get('injection_recovered') if arm == 'injected' else None
        rows.append(row)
    return rows


def peak_rows(folder, receipt):
    rows = []
    for arm, info in receipt.get('arms', {}).items():
        path = folder / (arm + '.json')
        if not path.exists():
            continue
        for ix, p in enumerate(load(path)['peaks']):
            item = dict(p); item['peak_rank'] = item.pop('rank', None)
            item.update(rank=receipt['rank'], source=receipt['source'], group=receipt['group'], arm=arm, peak_index=ix)
            rows.append(item)
    return rows


def attempt_row(receipt):
    return {k: v for k, v in receipt.items() if k not in ('arms', 'split')}


def write_source_rows(folder, receipt):
    folder = pathlib.Path(folder)
    save_bytes(folder / 'peak-rows.jsonl', ''.join(_dumps(r) + '\n' for r in peak_rows(folder, receipt)).encode())
    save_bytes(folder / 'arm-rows.jsonl', ''.join(_dumps(r) + '\n' for r in arm_rows(folder, receipt)).encode())


def _index(out):
    p = out / 'reports' / 'index.json'
    return load(p) if p.exists() else dict(ranks=[], lengths={name: 0 for name in AGGREGATES})


def export(out, rebuild=False):
    """Append rows of newly terminal sources. Returns the number of sources added."""
    out = pathlib.Path(out); reports = out / 'reports'; reports.mkdir(parents=True, exist_ok=True)
    with flock(reports / 'report.lock'):
        if rebuild:
            for name in AGGREGATES:
                (reports / name).unlink(missing_ok=True)
            (reports / 'index.json').unlink(missing_ok=True)
        index = _index(out); done = set(index['ranks'])
        for name in AGGREGATES:
            path = reports / name
            if path.exists() and path.stat().st_size != index['lengths'][name]:
                os.truncate(path, index['lengths'][name])   # uncommitted tail from an interrupted export
        new = []
        for receipt_path in sorted((out / 'sources').glob('*/receipt.json')):
            folder = receipt_path.parent
            if int(folder.name) in done:
                continue
            receipt = load(receipt_path)
            if not (folder / 'arm-rows.jsonl').exists():
                write_source_rows(folder, receipt)            # worker died between receipt and rows
            new.append((folder, receipt))
        if not new:
            return 0
        handles = {name: open(reports / name, 'ab') for name in AGGREGATES}
        try:
            for folder, receipt in new:
                handles['attempts.jsonl'].write((_dumps(attempt_row(receipt)) + '\n').encode())
                handles['arm-results.jsonl'].write((folder / 'arm-rows.jsonl').read_bytes())
                handles['peaks.jsonl'].write((folder / 'peak-rows.jsonl').read_bytes())
            for h in handles.values():
                h.flush(); os.fsync(h.fileno())
        finally:
            for h in handles.values():
                h.close()
        index['ranks'] = sorted(done | {int(f.name) for f, _ in new})
        index['lengths'] = {name: (reports / name).stat().st_size for name in AGGREGATES}
        index['wall'] = time.time()
        save(reports / 'index.json', index)
        return len(new)


def read_jsonl(path):
    path = pathlib.Path(path)
    if not path.exists():
        return
    with path.open() as f:
        for line in f:
            if line.strip():
                yield json.loads(line)


def summary(out, batch, status='RUNNING', stop=None):
    """Counts from the small aggregates only (attempts + arm rows); never touches peaks.jsonl."""
    out = pathlib.Path(out); reports = out / 'reports'
    attempts = list(read_jsonl(reports / 'attempts.jsonl'))
    arms = list(read_jsonl(reports / 'arm-results.jsonl'))
    counts = {}
    for spec in batch.arms:
        a = [x for x in arms if x['arm'] == spec['name'] and x['evidence_status'] == 'searched']
        sur = [x for x in a if x['survivors']]
        counts[spec['name']] = dict(completed_arms=len(a), survivor_peaks=sum(x['survivors'] for x in a),
                                    survivor_sources=len(sur), survivor_dependency_groups=len({x['group'] for x in sur}))
    completed = [r for r in attempts if r['status'] == 'complete']
    injected = [r for r in completed if 'injection_detected' in r]
    stages = sorted(int(p.stem.split('-')[1]) for p in out.glob('stage-*-freeze.json'))
    holdouts = [load(p) for p in (out / 'sources').glob('*/holdout-*.json')]
    s = dict(status=status, stop=stop, batch=batch.name, selected=batch.stages[-1], attempted=len(attempts),
             completed=len(completed), failed=sum(r['status'] == 'failed' for r in attempts),
             interrupted=sum(r['status'] == 'interrupted_incomplete' for r in attempts),
             attempted_dependency_groups=len({r['group'] for r in attempts}), counts=counts,
             injection_recovery_denominator=len(injected), injection_detected=sum(bool(r.get('injection_detected')) for r in injected),
             injection_recovered=sum(bool(r.get('injection_recovered')) for r in injected), frozen_stages=stages,
             holdouts_evaluated=len(holdouts), holdout_status_counts=pd.Series([h['status'] for h in holdouts], dtype=object).value_counts().to_dict(),
             wall=time.time())
    return s


def lead_groups(out, holdouts=True):
    """Lead grouping, streaming peaks.jsonl once and reading discovery times only for lead sources."""
    from .core import ephemeris_overlap
    out = pathlib.Path(out); leads = []
    for p in read_jsonl(out / 'reports' / 'peaks.jsonl'):
        if p['arm'] == 'real' and p['diagnostics']['photometry_vetting_pass']:
            folder = out / 'sources' / source_key(p['rank']); ho = folder / f"holdout-{p['peak_index']}.json"
            peak = {k: v for k, v in p.items() if k not in ('rank', 'source', 'group', 'arm', 'peak_index')}
            peak['rank'] = peak.pop('peak_rank')
            leads.append(dict(rank=p['rank'], source=p['source'], group=p['group'], peak_index=p['peak_index'], peak=peak,
                              discovery_status='screen_survivor', holdout=load(ho) if ho.exists() else dict(status='sealed' if holdouts else 'not_applicable_holdouts_disabled'),
                              host_status='unresolved', novelty_status='unresolved'))
    leads.sort(key=lambda x: (-x['peak']['snr'], x['source'], x['peak_index'])); groups = []; times = {}
    for lead in leads:
        if lead['rank'] not in times:
            times[lead['rank']] = pd.read_csv(out / 'sources' / source_key(lead['rank']) / 'discovery.csv.gz', usecols=['time']).time.to_numpy()
        found = None
        for g in groups:
            if g['dependency_group'] == lead['group'] and ephemeris_overlap(lead['peak'], g['representative']['peak'], times[lead['rank']])['possible_ephemeris_contamination']:
                found = g; break
        if found: found['members'].append(lead)
        else: groups.append(dict(lead_group=len(groups) + 1, dependency_group=lead['group'], representative=lead, members=[lead]))
    return groups


def full(out, batch, status='RUNNING', stop=None):
    """Incremental export + summary + lead groups + CSVs. Used by `deimos report` and at finalization."""
    out = pathlib.Path(out); reports = out / 'reports'
    export(out)
    s = summary(out, batch, status, stop)
    groups = lead_groups(out, batch.holdouts)
    with (reports / 'ranked-lead-groups.jsonl').open('w') as f:
        for g in groups:
            f.write(_dumps(g) + '\n')
    pd.DataFrame([dict(lead_group=g['lead_group'], dependency_group=g['dependency_group'], members=len(g['members']), source=g['representative']['source'], period_days=g['representative']['peak']['period_days'], snr=g['representative']['peak']['snr'], holdout_status=g['representative']['holdout']['status'], host_status='unresolved', novelty_status='unresolved') for g in groups]).to_csv(reports / 'ranked-lead-groups.csv', index=False)
    pd.DataFrame(list(read_jsonl(reports / 'attempts.jsonl'))).to_csv(reports / 'attempts.csv', index=False)
    pd.DataFrame(list(read_jsonl(reports / 'arm-results.jsonl'))).to_csv(reports / 'arm-results.csv', index=False)
    s.update(lead_groups=len(groups),
             holdout_supported_groups=sum(any(m['holdout']['status'] == 'supported_fixed_photometry_host_unresolved' for m in g['members']) for g in groups),
             unresolved_lead_groups=sum(any(m['holdout']['status'] in ['sealed', 'not_applicable_holdouts_disabled', 'insufficient_power', 'interrupted_holdout_unresolved', 'holdout_previously_accessed_by_other_batch'] for m in g['members']) for g in groups))
    prov = out / 'provenance.json'
    if prov.exists():
        s['provenance'] = load(prov)
    save(out / 'summary.json', s)
    (out / 'REPORT.md').write_text(
        f"# {batch.name}: paired Kepler Bonus search\n\nStatus: {status}. Stop: {stop}. {s['attempted']}/{s['selected']} sources terminal; "
        f"{s['completed']} complete, {s['failed']} failed, {s['interrupted']} interrupted (never replayed).\n\n```json\n{json.dumps(s['counts'], indent=2)}\n```\n\n"
        f"Injection detected/recovered: {s['injection_detected']}/{s['injection_recovered']} of {s['injection_recovery_denominator']} completed sources. "
        f"Lead groups: {s['lead_groups']}; holdout-supported groups: {s['holdout_supported_groups']}; unresolved: {s['unresolved_lead_groups']}.\n\n"
        "Source-level denominators are correlated; dependency groups are conservative and not certified independent. Real/control contrasts describe "
        "workload, not discovery yield. One flux-level injection per source cannot validate extraction. Host assignment, completeness, FPP and novelty "
        "remain unresolved; no promotion. Rows: reports/*.jsonl (one row per source per arm; peaks.jsonl one row per peak).\n"
        "These are screening outputs, not confirmed planets or validated candidates. Provenance (input sha256, config/code hashes, seeds, package version): provenance.json.\n")
    return s
