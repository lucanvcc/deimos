"""deimos fetch | validate | run | status | report | verify-pins  (worker is internal)."""
import os, sys, json, argparse, collections
from . import __version__
from .manifest import Batch, ManifestError
from .util import load, alive, check_pins, ENV_STATE, ENV_INPUT


def _batch(args):
    if getattr(args, 'input_dir', None):
        os.environ[ENV_INPUT] = os.pathsep.join(os.path.abspath(p) for p in args.input_dir)
    if getattr(args, 'state_dir', None):
        os.environ[ENV_STATE] = os.path.abspath(args.state_dir)
    return Batch(args.batch)


def cmd_validate(args):
    batch = _batch(args); rows = batch.validate()
    missing = []
    from .inputs import find_input, InputNotFound
    for row in rows:
        try:
            find_input(row, batch.input_dirs)
        except InputNotFound:
            missing.append(row['file_name'])
    print(json.dumps(dict(manifest=batch.path.name, sources=len(rows), stages=batch.stages, holdouts=batch.holdouts, arms=[a['name'] for a in batch.arms],
                          config_hash=batch.config_hash(), inputs_found=len(rows) - len(missing), inputs_missing=missing, state_dir='<output_dir>/state (default)' if batch.state_dir == batch.output_dir / 'state' else 'explicit'), indent=2))
    return 0 if not missing else 2


def cmd_fetch(args):
    """The only command that uses the network. Downloads the selected rows' files into the first input dir."""
    from . import inputs
    batch = _batch(args); rows = batch.validate()
    if not batch.input_dirs:
        raise ManifestError('no destination: set "input_dirs" in the manifest or pass --input-dir')
    dest = batch.input_dirs[0]
    if args.dry_run:
        for r in rows:
            print(f"would fetch {r['file_name']} from {inputs.candidate_urls(r['url'])[0]} into {dest.name}/")
        return 0
    result = inputs.fetch(rows, dest, log=lambda m: print(m, file=sys.stderr))
    print(json.dumps(result, indent=2))
    return 1 if result['failed'] else 0


def cmd_run(args):
    from . import coordinator
    batch = _batch(args)
    if args.detach:
        print(json.dumps(coordinator.launch_detached(batch)))
        return 0
    term = coordinator.run(batch)
    print(json.dumps(term))
    failed = [p.parent.name for p in sorted(batch.output_dir.glob('sources/*/receipt.json')) if load(p)['status'] != 'complete']
    if failed:
        print(f'error: {len(failed)} source(s) did not complete ({", ".join(failed)}); see sources/NNNN/receipt.json', file=sys.stderr)
        return 1
    return 0


def cmd_status(args):
    from . import coordinator
    batch = _batch(args); batch.validate(); out = batch.output_dir
    meta = load(out / 'coordinator.json') if (out / 'coordinator.json').exists() else {}
    running = bool(meta) and alive(meta.get('pid')) and not coordinator.lock_free(out)
    receipts = [load(p) for p in sorted((out / 'sources').glob('*/receipt.json'))] if (out / 'sources').exists() else []
    statuses = collections.Counter(r['status'] for r in receipts)
    info = dict(batch=batch.name, coordinator_running=running, terminal=load(out / 'terminal.json') if (out / 'terminal.json').exists() else None,
                selected=batch.stages[-1], terminal_sources=len(receipts), statuses=dict(statuses),
                frozen_stages=sorted(int(p.stem.split('-')[1]) for p in out.glob('stage-*-freeze.json')),
                interrupted_attempts=[p.parent.name for p in out.glob('sources/*/attempt.json') if not (p.parent / 'receipt.json').exists()])
    print(json.dumps(info, indent=2, default=str))
    return 0


def cmd_report(args):
    from . import report
    batch = _batch(args); out = batch.output_dir
    if args.rebuild:
        report.export(out, rebuild=True)
    term = load(out / 'terminal.json') if (out / 'terminal.json').exists() else {}
    s = report.full(out, batch, term.get('status', 'RUNNING'), term.get('stop'))
    print(json.dumps({k: s[k] for k in ('status', 'attempted', 'completed', 'lead_groups', 'counts')}, indent=2))
    return 0


def cmd_verify_pins(args):
    print(json.dumps(dict(ok=True, pins_digest=check_pins())))
    return 0


def cmd_worker(args):
    from .coordinator import worker
    worker(_batch(args), args.rank)
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(prog='deimos', description='Paired-control transit search on local Kepler Bonus light curves (screening only)')
    p.add_argument('--version', action='version', version=f'deimos {__version__}')
    sub = p.add_subparsers(dest='command', required=True)
    loc = argparse.ArgumentParser(add_help=False)
    loc.add_argument('batch', help='path to the manifest (batch.json)')
    loc.add_argument('--input-dir', action='append', help='folder holding the registered FITS files (repeatable; also DEIMOS_INPUT_DIR)')
    loc.add_argument('--state-dir', help='holdout-ledger folder (also DEIMOS_STATE); required when the manifest has holdouts')
    f = sub.add_parser('fetch', parents=[loc], help='download the registered FITS files into the first input dir (the only command that uses the network)')
    f.add_argument('--dry-run', action='store_true', help='list what would be fetched'); f.set_defaults(func=cmd_fetch)
    v = sub.add_parser('validate', parents=[loc], help='check manifest, registration CSV and that every input file is present'); v.set_defaults(func=cmd_validate)
    r = sub.add_parser('run', parents=[loc], help='start or resume a batch (a second coordinator on the same output folder is refused)')
    r.add_argument('--detach', action='store_true', help='run in the background; log to <output_dir>/job.log')
    r.set_defaults(func=cmd_run)
    s = sub.add_parser('status', parents=[loc], help='progress and coordinator liveness'); s.set_defaults(func=cmd_status)
    rp = sub.add_parser('report', parents=[loc], help='incremental export + summary.json, REPORT.md, CSVs and ranked lead groups')
    rp.add_argument('--rebuild', action='store_true', help='rebuild reports/*.jsonl from per-source rows'); rp.set_defaults(func=cmd_report)
    vp = sub.add_parser('verify-pins', help='check the frozen numerical files against frozen-pins.json'); vp.set_defaults(func=cmd_verify_pins)
    w = sub.add_parser('worker', parents=[loc]); w.add_argument('rank', type=int); w.set_defaults(func=cmd_worker)
    args = p.parse_args(argv)
    try:
        return args.func(args)
    except (ManifestError, FileNotFoundError, ValueError, RuntimeError) as e:
        print(f'error: {e}', file=sys.stderr)
        return 1
