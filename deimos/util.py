"""Small shared helpers: atomic JSON, hashing, locks, pin check, child environment.
"""
import os, json, hashlib, pathlib, fcntl, contextlib, subprocess

PACKAGE = pathlib.Path(__file__).resolve().parent
PINS = PACKAGE / 'frozen-pins.json'
ENV_STATE = 'DEIMOS_STATE'          # holdout ledger directory
ENV_INPUT = 'DEIMOS_INPUT_DIR'      # os.pathsep-separated folders holding the registered FITS files


def save(path, obj):
    """Atomic JSON write ."""
    path = pathlib.Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    with temp.open('w') as f:
        json.dump(obj, f, indent=2, default=str, allow_nan=False); f.flush(); os.fsync(f.fileno())
    temp.replace(path)


def save_bytes(path, data):
    path = pathlib.Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + f'.{os.getpid()}.tmp')
    with temp.open('wb') as f:
        f.write(data); f.flush(); os.fsync(f.fileno())
    temp.replace(path)


def write(path, obj):
    from .core import serial
    save(path, serial(obj))


def load(path):
    return json.loads(pathlib.Path(path).read_text())


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False, default=str)


def digest(obj):
    return hashlib.sha256(canonical(obj).encode()).hexdigest()


def source_key(rank):
    return f'{int(rank):04d}'


@contextlib.contextmanager
def flock(path, blocking=True):
    """Exclusive advisory lock. Non-blocking mode raises BlockingIOError if held."""
    path = pathlib.Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | (0 if blocking else fcntl.LOCK_NB))
        try:
            yield handle
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def alive(pid):
    if not pid:
        return False
    out = subprocess.run(['ps', '-p', str(pid), '-o', 'stat='], capture_output=True, text=True).stdout.strip()
    return bool(out) and not out.startswith('Z')


def check_pins():
    """Fail closed if a pinned numerical file differs from its pinned digest."""
    pins = load(PINS)['files']
    for rel, expected in pins.items():
        if sha(PACKAGE / rel) != expected:
            raise ValueError('Frozen input drift: ' + rel)
    return digest(pins)


THREAD_VARS = ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS']


def child_env(threads=1):
    env = os.environ.copy()
    for key in THREAD_VARS:
        env[key] = str(threads)
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    paths = [str(PACKAGE.parent)] + [p for p in env.get('PYTHONPATH', '').split(os.pathsep) if p]
    env['PYTHONPATH'] = os.pathsep.join(dict.fromkeys(paths))
    return env
