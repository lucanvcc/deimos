"""Stage freezes and the once-only sealed holdout .

Holdout rules:
- A holdout is evaluated only for a lead listed in the frozen holdout stage, at the frozen ephemeris, with no refit.
- An access record is written *before* the sealed file is opened; an access without a result is never retried.
- A ledger (<state_dir>/holdout-ledger/, shared by every batch pointing at the same state_dir) keyed by input FITS sha256 refuses any second batch that
  wants the same source's holdout. Sealed holdouts are never cached or copied between batches.
"""
import time, pathlib
import pandas as pd
from .core import diagnostics, fixed_holdout
from .core import clean_frame
from .core import transit_masks
from .util import load, write, sha, flock, source_key


def ledger_dir(batch):
    return batch.state_dir / 'holdout-ledger'
DEFAULT_RANKING_RULE = 'descending discovery SNR, then source designation, then peak index; real arm photometry_vetting_pass only'


def stage_path(out, stage):
    return pathlib.Path(out) / f'stage-{stage}-freeze.json'


def stage_freeze(batch, stage, rows):
    """Freeze the ranked real-arm leads of the first `stage` selected sources. Idempotent."""
    out = batch.output_dir; path = stage_path(out, stage)
    if path.exists(): return load(path)
    if len(rows) != stage: raise ValueError('Cannot freeze stage before every terminal attempt')
    leads = []
    for r in rows:
        p = out / 'sources' / source_key(r['rank']) / 'real.json'
        if not p.exists(): continue
        for ix, peak in enumerate(load(p)['peaks']):
            if peak['diagnostics']['photometry_vetting_pass']: leads.append(dict(rank=r['rank'], source=r['source'], group=r['group'], peak_index=ix, peak=peak))
    leads.sort(key=lambda x: (-x['peak']['snr'], x['source'], x['peak_index']))
    packet = dict(stage=stage, registration_hash=sha(batch.registration_csv), input_receipts={r['rank']: sha(out / 'sources' / source_key(r['rank']) / 'receipt.json') for r in rows},
                  leads=leads, ranking_rule=batch.raw.get('ranking_rule', DEFAULT_RANKING_RULE), wall=time.time())
    write(path, packet); return load(path)


def batch_label(out):
    """Ledger identity of a batch folder: its name plus a short hash of its resolved path (no absolute path is stored)."""
    import hashlib
    out = pathlib.Path(out).resolve()
    return f'{out.name}-{hashlib.sha256(str(out).encode()).hexdigest()[:12]}'


def _claim_ledger(ledger, input_sha256, batch_id, peak_index):
    """Record this batch's access in the machine-wide ledger. False if another batch already accessed this source."""
    ledger.mkdir(parents=True, exist_ok=True)
    entry = ledger / f'{input_sha256}.json'
    with flock(ledger / 'ledger.lock'):
        record = load(entry) if entry.exists() else dict(input_sha256=input_sha256, accesses=[])
        if any(a['batch'] != batch_id for a in record['accesses']):
            return False
        record['accesses'].append(dict(batch=batch_id, peak_index=peak_index, wall=time.time()))
        write(entry, record)
    return True


def holdout(batch, lead, stage):
    cfg = batch.science_config(); out = batch.output_dir
    folder = out / 'sources' / source_key(lead['rank']); res = folder / f"holdout-{lead['peak_index']}.json"; access = folder / f"access-{lead['peak_index']}.json"; fp = stage_path(out, stage)
    if res.exists(): return load(res)
    if access.exists(): return dict(status='interrupted_holdout_unresolved', refit=False)
    packet = load(fp)
    if lead not in packet['leads']: raise ValueError('Lead absent from frozen stage')
    split = load(folder / 'split.json')
    if 'holdout_sha256' not in split: raise ValueError('No sealed holdout exists for this source')
    if sha(folder / 'sealed-holdout.csv.gz') != split['holdout_sha256']: raise ValueError('Sealed holdout drift')
    if not _claim_ledger(ledger_dir(batch), split['input_sha256'], batch_label(out), lead['peak_index']):
        result = dict(status='holdout_previously_accessed_by_other_batch', refit=False, stage_freeze_hash=sha(fp)); write(res, result); return result
    write(access, dict(stage_freeze_hash=sha(fp), wall=time.time(), once_only=True, followup=False))
    h = clean_frame(pd.read_csv(folder / 'sealed-holdout.csv.gz'), cfg['clean']); p = lead['peak']; ep = {k: p[k] for k in ['period_days', 'epoch_bkjd', 'duration_days']}; ev = fixed_holdout(h, ep)
    # Same repair diagnostics at frozen timing; no search/refit.
    proxy = dict(p); proxy['initial_candidate'] = bool(ev.get('passes_frozen_threshold', False)); proxy['snr'] = ev.get('snr', 0) or 0; diag = diagnostics(h, proxy, cfg['validate']) if ev.get('status') == 'evaluated' else None
    _, _, inside = transit_masks(h.time.to_numpy(), ep['period_days'], ep['epoch_bkjd'], ep['duration_days'])
    error = ev.get('error'); pred = p['depth'] / error if error else None; pass_all = bool(ev.get('passes_frozen_threshold') and diag and diag['photometry_vetting_pass'])
    status = 'supported_fixed_photometry_host_unresolved' if pass_all else 'insufficient_power' if pred is None or pred < 7 or ev.get('observed_transits', 0) < 3 else 'failed_frozen_screen'
    result = dict(status=status, evaluation=ev, diagnostics=diag, predicted_discovery_depth_snr=pred, depth_uncertainty=error, covered_in_samples=int(inside.sum()), refit=False, upstream_independence=False, stage_freeze_hash=sha(fp)); write(res, result); return result
