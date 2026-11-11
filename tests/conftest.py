"""Shared fixtures. Everything runs in tmp_path: no network, no real data, no API keys."""
import json, os, pathlib, socket, sys

import numpy as np, pandas as pd, pytest
from astropy.io import fits

ROOT = pathlib.Path(__file__).resolve().parents[1]
_paths = [str(ROOT)]
for p in _paths:
    if p not in sys.path:
        sys.path.insert(0, p)
# Subprocess workers must import the same sources (and the tests package for the fake plugin).
os.environ['PYTHONPATH'] = os.pathsep.join(_paths + [p for p in os.environ.get('PYTHONPATH', '').split(os.pathsep) if p])

REG_COLUMNS = ['rank', 'gaia_designation', 'conservative_group', 'tpforg', 'kic', 'ra', 'dec', 'file_name', 'url', 'sha256', 'role']
QUARTERS = (1, 2, 3, 4, 5, 6)


@pytest.fixture(autouse=True)
def _isolation(monkeypatch):
    for var in ('DEIMOS_STATE', 'DEIMOS_INPUT_DIR'):
        monkeypatch.delenv(var, raising=False)

    def no_network(*a, **k):
        raise AssertionError('network access attempted in a test')
    monkeypatch.setattr(socket.socket, 'connect', no_network)


def synthetic_discovery(rank, quarters=(1, 2, 3), n=1500):
    rng = np.random.default_rng(rank)
    frames = []
    for i, q in enumerate(quarters):
        t = 100 + 95 * i + np.arange(n) * 0.0204
        frames.append(pd.DataFrame(dict(time=t, flux=1 + rng.normal(0, 1e-3, n), flux_err=np.full(n, 1e-3), quarter=q)))
    return pd.concat(frames, ignore_index=True)


def make_fits(path, gaia='Gaia DR3 1', seed=1, quarters=QUARTERS, n=1500, period=None, depth=0.004, duration=0.15):
    """A synthetic Kepler-Bonus-shaped file (HLSPID KBONUS-BKG, LIGHTCURVE_Qnn extensions) that the adapter accepts."""
    rng = np.random.default_rng(seed)
    primary = fits.PrimaryHDU()
    primary.header['HLSPID'] = 'KBONUS-BKG'; primary.header['GAIAID'] = gaia
    hdus = [primary]
    for i, q in enumerate(quarters):
        t = 100 + 95 * i + np.arange(n) * 0.0204
        flux = 30000 * (1 + rng.normal(0, 1e-3, n))
        if period:
            phase = ((t - 101.0 + 0.5 * period) % period) - 0.5 * period
            flux = flux * np.where(np.abs(phase) < duration / 2, 1 - depth, 1.0)
        err = np.full(n, 30000 * 1e-3)
        cols = [fits.Column(name='TIME', format='D', unit='BJD - 2454833', array=t),
                fits.Column(name='FLUX', format='D', unit='e-/s', array=flux),
                fits.Column(name='FLUX_ERR', format='D', unit='e-/s', array=err),
                fits.Column(name='SAP_QUALITY', format='J', array=np.zeros(n, dtype=np.int32))]
        h = fits.BinTableHDU.from_columns(cols, name=f'LIGHTCURVE_Q{q}')
        h.header['TIMESYS'] = 'TDB'; h.header['BJDREFI'] = 2454833; h.header['BJDREFF'] = 0.0; h.header['TIMEZERO'] = 0.0; h.header['PSFFRAC'] = 0.8
        hdus.append(h)
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return pathlib.Path(path)


def sha_of(path):
    import hashlib
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


@pytest.fixture
def workspace(tmp_path):
    """tmp_path/data holds three synthetic FITS files; tmp_path/registration.csv registers them (with sha256 pins)."""
    data = tmp_path / 'data'; data.mkdir()
    ranks = [11, 12, 13]; rows = []
    for rank in ranks:
        name = f'synthetic_{rank}.fits'
        p = make_fits(data / name, gaia=f'Gaia DR3 {rank}', seed=rank)
        rows.append(dict(rank=rank, gaia_designation=f'Gaia DR3 {rank}', conservative_group=rank % 2, tpforg=rank, kic=rank, ra=1.0, dec=2.0,
                         file_name=name, url='https://example.invalid/uri=mast:HLSP/x', sha256=sha_of(p), role='selected'))
    pd.DataFrame(rows, columns=REG_COLUMNS).to_csv(tmp_path / 'registration.csv', index=False)
    return tmp_path, ranks


def write_manifest(path, ws, out='out', **extra):
    m = dict(name='test-batch', registration_csv=str(ws / 'registration.csv'), output_dir=str(ws / out), input_dirs=[str(ws / 'data')],
             seed_base=2026101457, holdouts=False, workers=2, search_plugin='tests.fakeplugin:select',
             arms=['real', 'inverted', 'block_scrambled', 'injected'])
    m.update(extra)
    path = pathlib.Path(path); path.write_text(json.dumps(m))
    return path
