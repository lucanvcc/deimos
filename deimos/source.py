"""One source through every arm; the receipt records the package version and arm seeds."""
import json, time, resource
from .core import split_quarters
from . import __version__, science, report
from .util import write, load, sha, source_key
from .adapter import adapt
from .inputs import find_input
from .search import load_plugin


class Context:
    """Everything a worker needs, built once per process."""

    def __init__(self, batch):
        self.batch = batch
        self.out = batch.output_dir
        self.cfg = batch.science_config()
        self.config_hash = batch.config_hash()
        self.select, plugin_module = load_plugin(batch.search_plugin)
        self.code_hash = science.code_hash(plugin_module)

    def folder(self, rank):
        return self.out / 'sources' / source_key(rank)


def _prepare_from_fits(ctx, row, folder, r):
    b = ctx.batch
    path, digest = find_input(row, b.input_dirs)
    r['transport'] = 'local_file'; r['input_file'] = path.name; r['bytes'] = path.stat().st_size; r['input_sha256'] = digest
    raw, meta = adapt(path); gaia = str(meta['primary_header'].get('GAIAID', '')).replace('Gaia DR3 ', '').strip()
    if gaia and gaia != row['gaia_designation'].replace('Gaia DR3 ', ''): raise ValueError('source GAIAID identity mismatch')
    d, h, plan = split_quarters(raw); r['split'] = plan
    cols = ['time', 'flux', 'flux_err', 'quarter']; d = d[cols].reset_index(drop=True); h = h[cols].reset_index(drop=True)
    d.to_csv(folder / 'discovery.csv.gz', index=False)
    split = dict(split=plan, input_sha256=r['input_sha256'], quarter_audit=meta['quarter_audit'], upstream_independence=False)
    if b.holdouts:
        # Sealed holdout lives only in this batch's source folder; never shared between batches.
        h.to_csv(folder / 'sealed-holdout.csv.gz', index=False); split['holdout_sha256'] = sha(folder / 'sealed-holdout.csv.gz')
    else:
        split['holdout'] = 'not produced (holdouts disabled)'
    write(folder / 'split.json', split)
    del h
    return d


def run_arms(ctx, row, folder, r, d):
    b = ctx.batch; rank = int(row['rank']); seed = science.arm_seed(b.seed_base, rank)
    truth = science.truth_for(rank, d, b.seed_base, b.injection)
    if any(a['name'] == 'injected' for a in b.arms):
        write(folder / 'private-injection-truth.json', truth)
    times = d.time.to_numpy()
    for spec in b.arms:
        arm = spec['name']; armfile = folder / (arm + '.json')
        if armfile.exists(): raise ValueError('unexpected preexisting arm without receipt')
        x = science.arm_frame(d, arm, rank, b.seed_base, truth, b.injection)
        result = ctx.select(x, ctx.cfg); write(armfile, result)
        r['arms'][arm] = dict(peaks=len(result['peaks']), survivors=sum(p['diagnostics']['photometry_vetting_pass'] for p in result['peaks']), seconds=result['seconds'], sha256=sha(armfile), seed=seed)
        if arm == 'injected':
            r['injection_detected'] = any(science.recovered(p, truth, times) and p['initial_candidate'] for p in result['peaks'])
            r['injection_recovered'] = any(science.recovered(p, truth, times) and p['diagnostics']['photometry_vetting_pass'] for p in result['peaks'])


def process(ctx, row):
    rank = int(row['rank']); folder = ctx.folder(rank); folder.mkdir(parents=True, exist_ok=True); receipt = folder / 'receipt.json'
    if receipt.exists(): return load(receipt)
    ident = dict(rank=rank, source=row['gaia_designation'], group=int(row['conservative_group']))
    attempt = folder / 'attempt.json'
    if attempt.exists():
        # Fail closed: an interrupted attempt is never replayed.
        write(receipt, dict(**ident, status='interrupted_incomplete', reason='No replay of outcome-bearing arm after interrupted attempt'))
        report.write_source_rows(folder, load(receipt)); return load(receipt)
    write(attempt, dict(**ident, wall=time.time()))
    start = time.monotonic(); r = dict(**ident, origin=int(row['tpforg']), kic=int(row['kic']), ra=float(row['ra']), dec=float(row['dec']), transport=False, arms={},
                                        code_hash=ctx.code_hash, config_hash=ctx.config_hash, seed_base=ctx.batch.seed_base, package_version=__version__)
    try:
        d = _prepare_from_fits(ctx, row, folder, r)
        run_arms(ctx, row, folder, r, d)
        r['status'] = 'complete'
    except Exception as e:
        r.update(status='failed', error_type=type(e).__name__, error=str(e))
    r.update(seconds=time.monotonic() - start, peak_resident_memory_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    write(receipt, r)
    report.write_source_rows(folder, r)
    print(json.dumps(dict(rank=rank, status=r['status'], seconds=r['seconds'], arms={k: {x: v[x] for x in ('peaks', 'survivors', 'seconds') if x in v} for k, v in r['arms'].items()})), flush=True)
    return r
