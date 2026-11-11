"""The fast BLS plug-in must reproduce astropy's kernel and frozen_select bit for bit."""
import ast, json, pathlib, shutil
import numpy as np, pytest
from astropy.timeseries import BoxLeastSquares
from deimos import fastbls, search, science
from tests.conftest import synthetic_discovery

CFG = json.loads((pathlib.Path(__file__).resolve().parents[1] / 'deimos' / 'science-config.json').read_text())


@pytest.fixture(autouse=True)
def _compiler(tmp_path, monkeypatch):
    if not shutil.which(fastbls._compiler()):
        pytest.skip('no C compiler')
    monkeypatch.setattr(fastbls, 'BUILD_DIR', tmp_path)
    monkeypatch.setattr(fastbls, '_LIB', None)


def _bits(a):
    return np.ascontiguousarray(a, dtype=np.float64).view(np.int64)


def _light_curve(seed, n=3000, gaps=True, transit=True):
    rng = np.random.default_rng(seed)
    t = np.concatenate([50 + 93 * q + np.arange(n // 3) * 0.0204 for q in range(3)])
    if gaps:
        t = np.sort(np.delete(t, rng.choice(len(t), len(t) // 10, replace=False)))
    dy = rng.uniform(3e-4, 3e-3, len(t))
    y = 1 + rng.normal(0, 1e-3, len(t))
    if transit:
        y -= 4e-3 * (np.abs((t - 52.1 + 2.35) % 4.7 - 2.35) < 0.07)
    return t, y, dy


@pytest.mark.parametrize('seed,objective,gaps,transit', [(1, 'snr', True, True), (2, 'snr', False, False), (3, 'likelihood', True, True), (4, 'snr', True, False)])
def test_kernel_bit_identical_to_astropy(seed, objective, gaps, transit):
    t, y, dy = _light_curve(seed, gaps=gaps, transit=transit)
    model = BoxLeastSquares(t, y, dy)
    periods = np.geomspace(0.5, 30, 6000); durations = np.linspace(0.05, 0.3, 8)
    ref = model.power(periods, durations, objective=objective)
    with fastbls.kernel():
        new = model.power(periods, durations, objective=objective)
    assert fastbls._bls_methods.bls_fast is fastbls.ASTROPY_BLS_FAST   # restored after the block
    ok = np.isfinite(ref.power)
    assert np.array_equal(np.isfinite(new.power), ok)
    for field in ('power', 'depth', 'depth_err', 'duration', 'transit_time', 'depth_snr', 'log_likelihood'):
        assert np.array_equal(_bits(getattr(new, field)[ok]), _bits(getattr(ref, field)[ok])), field


def test_kernel_rejects_invalid_inputs_like_astropy():
    t, y, dy = _light_curve(5)
    args = (t - t.min(), y - np.median(y), 1 / dy ** 2)
    with pytest.raises(ValueError):
        fastbls.bls_fast(*args, np.array([0.2, 1.0]), np.array([0.3]), 10, False)   # duration > min period
    with pytest.raises(ValueError):
        fastbls.ASTROPY_BLS_FAST(*args, np.array([0.2, 1.0]), np.array([0.3]), 10, False)


def test_exact_select_packet_identical_to_frozen_select():
    raw = synthetic_discovery(7, quarters=(1, 2, 3, 4), n=1200)
    t = raw.time.to_numpy()
    raw['flux'] *= 1 - 3e-3 * (np.abs((t - 101.1 + 1.6) % 3.2 - 1.6) < 0.05)
    strip = lambda p: {k: v for k, v in p.items() if k != 'seconds'}
    assert strip(fastbls.exact_select(raw, CFG)) == strip(search.frozen_select(raw, CFG))


def test_exact_select_body_is_frozen_select():
    def body(path, name):
        stmts = next(n for n in ast.parse(pathlib.Path(path).read_text()).body if isinstance(n, ast.FunctionDef) and n.name == name).body
        return stmts[1:] if isinstance(stmts[0], ast.Expr) and isinstance(stmts[0].value, ast.Constant) else stmts
    fast = body(fastbls.__file__, 'exact_select')
    assert isinstance(fast[0], ast.With) and ast.unparse(fast[0].items[0].context_expr) == 'kernel()'
    dump = lambda stmts: ast.dump(ast.Module(body=stmts, type_ignores=[]))
    assert dump(fast[0].body) == dump(body(search.__file__, 'frozen_select'))


def test_code_hash_covers_kernel_source():
    assert 'int fast_bls(' in pathlib.Path(fastbls.__file__).read_text()
    assert science.code_hash(fastbls) != science.code_hash(search)
