"""Batch coordinator: dispatch workers, freeze stages, evaluate holdouts once, export incrementally.

One coordinator per output folder (<output_dir>/coordinator.lock); a resume must match the batch identity.
"""
import os, sys, json, time, subprocess, fcntl, pathlib, platform
from .util import save, load, sha, child_env, check_pins, source_key
from . import report, stages, __version__


def processes():
    return subprocess.run(['ps', '-axo', 'pid=,command='], capture_output=True, text=True).stdout.splitlines()


def orphan_workers(batch):
    tag = f'deimos worker {batch.path} '
    return [line.strip() for line in processes() if tag in line]


class AlreadyRunning(RuntimeError):
    pass


class Locks:
    """Hold the per-batch coordinator lock for the life of the process."""

    def __init__(self, out):
        self.paths = [pathlib.Path(out) / 'coordinator.lock']
        self.handles = []

    def __enter__(self):
        for p in self.paths:
            p.parent.mkdir(parents=True, exist_ok=True)
            h = p.open('a+')
            try:
                fcntl.flock(h, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                h.close(); self.__exit__()
                raise AlreadyRunning(f'another coordinator holds {p}') from None
            self.handles.append(h)
        return self

    def __exit__(self, *exc):
        for h in self.handles:
            fcntl.flock(h, fcntl.LOCK_UN); h.close()
        self.handles = []


def lock_free(out):
    try:
        with Locks(out):
            return True
    except AlreadyRunning:
        return False


def check_identity(batch):
    out = batch.output_dir; out.mkdir(parents=True, exist_ok=True)
    ident = batch.identity(); path = out / 'batch-identity.json'
    if path.exists():
        old = load(path)
        if json.loads(json.dumps(old)) != json.loads(json.dumps(ident, default=str)):
            diff = sorted(k for k in set(old) | set(ident) if json.dumps(old.get(k), default=str) != json.dumps(ident.get(k), default=str))
            raise ValueError(f'batch manifest changed scientific fields since this output_dir started: {diff}. Use a new output_dir.')
    else:
        save(path, json.loads(json.dumps(ident, default=str)))


def write_provenance(batch, ctx):
    """provenance.json: everything needed to rerun and to audit the numbers (no absolute paths)."""
    import numpy, pandas, scipy, astropy
    out = batch.output_dir
    inputs = {}
    for p in sorted((out / 'sources').glob('*/receipt.json')) if (out / 'sources').exists() else []:
        r = load(p)
        if 'input_sha256' in r:
            inputs[source_key(r['rank'])] = dict(file=r.get('input_file'), sha256=r['input_sha256'])
    save(out / 'provenance.json', dict(
        package='deimos', package_version=__version__, frozen_pins_digest=check_pins(), config_hash=ctx.config_hash, code_hash=ctx.code_hash,
        science_config_sha256=sha(batch.science_config_path), registration_sha256=sha(batch.registration_csv), seed_base=batch.seed_base,
        arm_seed_rule='seed_base + rank', stages=batch.stages, holdouts=batch.holdouts,
        python=sys.version.split()[0], platform=platform.platform(),
        libraries=dict(numpy=numpy.__version__, pandas=pandas.__version__, scipy=scipy.__version__, astropy=astropy.__version__),
        inputs=inputs))


def terminal(batch, status, stop):
    out = batch.output_dir; term = dict(status=status, stop=stop, wall=time.time())
    save(out / 'terminal.json', term)
    from .source import Context
    write_provenance(batch, Context(batch))
    report.full(out, batch, status, stop)
    return term


def run(batch):
    out = batch.output_dir
    batch.validate()
    out.mkdir(parents=True, exist_ok=True)
    check_pins()
    with Locks(out):
        save(out / 'coordinator.json', dict(pid=os.getpid(), manifest=batch.path.name, wall=time.time()))
        if orphan_workers(batch):
            raise RuntimeError('worker processes from a previous coordinator are still alive; wait for them to exit:\n' + '\n'.join(orphan_workers(batch)))
        check_identity(batch)
        if (out / 'terminal.json').exists():
            return load(out / 'terminal.json')
        return _loop(batch)


def worker_command(batch, rank):
    return [sys.executable, '-m', 'deimos', 'worker', str(batch.path), str(rank)]


def _loop(batch):
    out = batch.output_dir; rows = batch.selection()
    receipt = lambda rank: out / 'sources' / source_key(rank) / 'receipt.json'
    ranks = [int(r['rank']) for r in rows]
    pending = [r for r in ranks if not receipt(r).exists()]
    env = child_env(batch.threads_per_worker)
    children = {}; last_export = 0; last_status = 0
    save(out / 'launch.json', dict(pid=os.getpid(), workers=batch.workers, threads_per_worker=batch.threads_per_worker,
                                   pending=len(pending), python=sys.version.split()[0], wall=time.time()))
    from .source import Context
    write_provenance(batch, Context(batch))
    try:
        while pending or children:
            finished = False
            for pid, (child, rank) in list(children.items()):
                code = child.poll()
                if code is None: continue
                del children[pid]
                if code or not receipt(rank).exists():
                    raise RuntimeError(f'worker {rank} exited {code} without a receipt; no replay (it will be quarantined as interrupted on resume)')
                finished = True
            if finished:
                _freeze_ready(batch, rows)
            while pending and len(children) < batch.workers:
                rank = pending.pop(0)
                child = subprocess.Popen(worker_command(batch, rank), env=env)
                children[child.pid] = (child, rank)
            if finished or time.monotonic() - last_export > batch.export_every_seconds:
                report.export(out); save(out / 'summary.json', report.summary(out, batch)); last_export = time.monotonic()
            if time.monotonic() - last_status > 10:
                save(out / 'status.json', dict(pid=os.getpid(), workers=[dict(pid=c.pid, rank=r) for c, r in children.values()], pending=len(pending),
                                               wall=time.time()))
                last_status = time.monotonic()
            time.sleep(0.2)
    finally:
        for child, _ in children.values():
            if child.poll() is None:
                child.terminate(); child.wait(timeout=30)
    _freeze_ready(batch, rows)
    return terminal(batch, 'COMPLETE', None)


def _freeze_ready(batch, rows):
    """Freeze every stage whose first-N sources are all terminal; evaluate holdouts once after the holdout stage."""
    out = batch.output_dir
    for stage in batch.stages:
        if stages.stage_path(out, stage).exists():
            continue
        head = rows[:stage]
        paths = [out / 'sources' / source_key(r['rank']) / 'receipt.json' for r in head]
        if not all(p.exists() for p in paths):
            break
        packet = stages.stage_freeze(batch, stage, [load(p) for p in paths])
        print(json.dumps(dict(event='stage_frozen', stage=stage, leads=len(packet['leads']))), flush=True)
    if batch.holdouts and stages.stage_path(out, batch.holdout_stage).exists():
        packet = load(stages.stage_path(out, batch.holdout_stage))
        for lead in packet['leads']:
            stages.holdout(batch, lead, batch.holdout_stage)


def worker(batch, rank):
    from .source import Context, process
    check_pins()
    batch.validate()
    ctx = Context(batch)
    row = next(r for r in batch.selection() if int(r['rank']) == int(rank))
    process(ctx, row)


def launch_detached(batch):
    """Start `run` in its own session."""
    batch.validate()
    out = batch.output_dir; out.mkdir(parents=True, exist_ok=True)
    if not lock_free(out):
        raise AlreadyRunning('a coordinator is already running for this batch (see `status`)')
    log = open(out / 'job.log', 'ab', buffering=0)
    cmd = [sys.executable, '-m', 'deimos', 'run', str(batch.path)]
    job = subprocess.Popen(cmd, env=child_env(batch.threads_per_worker), stdout=log, stderr=log, start_new_session=True)
    return dict(pid=job.pid, log=str(out / 'job.log'))
