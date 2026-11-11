import json, subprocess, sys

import numpy as np, pandas as pd, pytest

from deimos import science, report, coordinator, cli, __version__
from deimos.manifest import Batch, DEFAULT_INJECTION
from deimos.inputs import find_input, InputMismatch, InputNotFound
from deimos.source import Context, process
from deimos.util import load, write
from tests.conftest import synthetic_discovery, write_manifest, make_fits, sha_of, ROOT

ARMS = ['real', 'inverted', 'block_scrambled', 'injected']


# --- arm science (by construction and by property) -------------------

def test_transform_properties():
    d = synthetic_discovery(7)
    inv = science.transform(d, 'inverted', 1)
    assert np.allclose(inv.flux + d.flux, 2)
    a, b = science.transform(d, 'block_scrambled', 5), science.transform(d, 'block_scrambled', 5)
    pd.testing.assert_frame_equal(a, b, check_exact=True)
    assert not a.flux.equals(d.flux) and sorted(a.flux) == sorted(d.flux)
    with pytest.raises(ValueError):
        science.transform(d, 'nonsense', 1)


def test_truth_default_grid_matches_previous_for_64_ranks():
    d = synthetic_discovery(3)
    for rank in range(1, 65):
        rng = np.random.default_rng(2026101457 + rank)
        cell = (rank - 1) % 8
        lo, hi = DEFAULT_INJECTION['jitter']
        period = DEFAULT_INJECTION['periods_days'][cell // 4] * rng.uniform(lo, hi)
        depth = DEFAULT_INJECTION['depths'][(cell // 2) % 2] * rng.uniform(lo, hi)
        duration = DEFAULT_INJECTION['durations_days'][cell % 2] * rng.uniform(lo, hi)
        expected = dict(period_days=period, depth=depth, duration_days=duration, epoch_bkjd=float(d.time.min() + rng.uniform(0, period)), cell=cell)
        actual = science.truth_for(rank, d, 2026101457, DEFAULT_INJECTION)
        assert actual == expected and json.dumps(actual).encode() == json.dumps(expected).encode()


def test_truth_expanded_grid_visits_all_60_cells():
    d = synthetic_discovery(3)
    injection = dict(periods_days=[2, 10, 15, 22, 28], depths=[5e-4, 2e-3, 5e-3, 1e-2], durations_days=[0.1, 0.2, 0.3], jitter=[1.0, 1.0])
    truths = [science.truth_for(rank, d, 2026101701, injection) for rank in range(1, 61)]
    assert [t['cell'] for t in truths] == list(range(60))


def test_code_hash_changes_with_plugin():
    import deimos.search as s, tests.fakeplugin as fp
    assert science.code_hash(s) != science.code_hash(fp)


# --- inputs --------------------------------------------------------------------------------------------

def test_find_input_verifies_pin(workspace):
    tmp, ranks = workspace
    row = pd.read_csv(tmp / 'registration.csv').iloc[0].to_dict()
    path, digest = find_input(row, [tmp / 'data'])
    assert digest == row['sha256'] and path.name == row['file_name']
    with pytest.raises(InputMismatch, match='does not match the registration pin'):
        find_input(dict(row, sha256='0' * 64), [tmp / 'data'])
    with pytest.raises(InputNotFound, match='deimos fetch'):
        find_input(dict(row, file_name='absent.fits'), [tmp / 'data'])


def test_missing_input_fails_the_source_clearly_and_validate_reports_it(workspace, capsys):
    tmp, ranks = workspace
    (tmp / 'data' / 'synthetic_12.fits').unlink()
    batch = Batch(write_manifest(tmp / 'm.json', tmp)); batch.validate()
    ctx = Context(batch)
    r = process(ctx, batch.selection()[1])
    assert r['status'] == 'failed' and r['error_type'] == 'InputNotFound' and 'synthetic_12.fits' in r['error']
    capsys.readouterr()
    assert cli.main(['validate', str(tmp / 'm.json')]) == 2
    assert 'synthetic_12.fits' in capsys.readouterr().out


def test_tampered_input_is_refused(workspace):
    tmp, ranks = workspace
    make_fits(tmp / 'data' / 'synthetic_11.fits', gaia='Gaia DR3 11', seed=999)   # different bytes, same name
    batch = Batch(write_manifest(tmp / 'm.json', tmp)); batch.validate()
    r = process(Context(batch), batch.selection()[0])
    assert r['status'] == 'failed' and r['error_type'] == 'InputMismatch'


def test_identity_mismatch_between_registration_and_file_is_refused(workspace):
    tmp, ranks = workspace
    reg = pd.read_csv(tmp / 'registration.csv'); reg.loc[0, 'gaia_designation'] = 'Gaia DR3 999999'; reg.to_csv(tmp / 'registration.csv', index=False)
    batch = Batch(write_manifest(tmp / 'm.json', tmp)); batch.validate()
    r = process(Context(batch), batch.selection()[0])
    assert r['status'] == 'failed' and 'GAIAID identity mismatch' in r['error']


def test_non_kbonus_file_is_refused(workspace):
    from astropy.io import fits
    tmp, ranks = workspace
    p = tmp / 'data' / 'synthetic_11.fits'
    with fits.open(p, mode='update') as h:
        h[0].header['HLSPID'] = 'OTHER'
    reg = pd.read_csv(tmp / 'registration.csv'); reg.loc[0, 'sha256'] = sha_of(p); reg.to_csv(tmp / 'registration.csv', index=False)
    batch = Batch(write_manifest(tmp / 'm.json', tmp)); batch.validate()
    r = process(Context(batch), batch.selection()[0])
    assert r['status'] == 'failed' and 'wrong archive family' in r['error']


# --- per-source processing and reporting ---------------------------------------------------------------

def test_process_writes_arms_and_receipt_with_provenance(workspace):
    tmp, ranks = workspace
    batch = Batch(write_manifest(tmp / 'm.json', tmp)); batch.validate()
    ctx = Context(batch)
    r = process(ctx, batch.selection()[0])
    assert r['status'] == 'complete' and set(r['arms']) == set(ARMS)
    assert r['input_sha256'] == sha_of(tmp / 'data' / 'synthetic_11.fits') and r['package_version'] == __version__
    assert r['arms']['real']['seed'] == 2026101457 + 11 and r['seed_base'] == 2026101457
    assert len(r['code_hash']) == 64 and len(r['config_hash']) == 64
    assert process(ctx, batch.selection()[0]) == load(tmp / 'out' / 'sources' / '0011' / 'receipt.json')   # finished sources are not redone


def test_interrupted_attempt_is_quarantined_not_replayed(workspace):
    tmp, ranks = workspace
    batch = Batch(write_manifest(tmp / 'm.json', tmp)); batch.validate()
    folder = tmp / 'out' / 'sources' / '0011'; folder.mkdir(parents=True)
    write(folder / 'attempt.json', dict(rank=11))
    r = process(Context(batch), batch.selection()[0])
    assert r['status'] == 'interrupted_incomplete' and not (folder / 'real.json').exists()


def test_incremental_export_appends_once_and_repairs_torn_tail(workspace):
    tmp, ranks = workspace
    batch = Batch(write_manifest(tmp / 'm.json', tmp, arms=['real', 'inverted'])); batch.validate()
    ctx = Context(batch); out = tmp / 'out'
    process(ctx, batch.selection()[0])
    assert report.export(out) == 1 and report.export(out) == 0
    rows = list(report.read_jsonl(out / 'reports' / 'arm-results.jsonl'))
    assert len(rows) == 2 and {r['arm'] for r in rows} == {'real', 'inverted'}
    with (out / 'reports' / 'arm-results.jsonl').open('a') as f:
        f.write('{"torn": tru')                                          # crash after append, before the index commit
    for row in batch.selection()[1:]:
        process(ctx, row)
    assert report.export(out) == 2
    rows = list(report.read_jsonl(out / 'reports' / 'arm-results.jsonl'))
    assert len(rows) == 2 * len(ranks) and sorted({r['rank'] for r in rows}) == ranks
    incremental = (out / 'reports' / 'arm-results.jsonl').read_bytes()
    report.export(out, rebuild=True)
    assert sorted(incremental.splitlines()) == sorted((out / 'reports' / 'arm-results.jsonl').read_bytes().splitlines())
    s = report.full(out, batch)
    assert s['attempted'] == len(ranks) and (out / 'reports' / 'ranked-lead-groups.csv').exists()


# --- end to end ------------------------------------------------------------------------------------------

def test_end_to_end_run_resume_identity_and_report(workspace):
    """Real subprocess workers, two at a time; a second run is a no-op resume; a changed manifest is refused."""
    tmp, ranks = workspace
    manifest = write_manifest(tmp / 'm.json', tmp)
    cmd = [sys.executable, '-m', 'deimos']
    p = subprocess.run(cmd + ['run', str(manifest)], capture_output=True, text=True, timeout=300)
    assert p.returncode == 0, p.stderr
    out = tmp / 'out'
    assert load(out / 'terminal.json')['status'] == 'COMPLETE'
    for name in ('summary.json', 'REPORT.md', 'provenance.json', 'stage-3-freeze.json', 'batch-identity.json'):
        assert (out / name).exists(), name
    for name in ('arm-results.csv', 'attempts.csv', 'ranked-lead-groups.csv', 'ranked-lead-groups.jsonl', 'peaks.jsonl'):
        assert (out / 'reports' / name).exists(), name
    assert len(list(report.read_jsonl(out / 'reports' / 'arm-results.jsonl'))) == 4 * len(ranks)
    prov = load(out / 'provenance.json')
    assert prov['package_version'] == __version__ and prov['seed_base'] == 2026101457 and set(prov['inputs']) == {'0011', '0012', '0013'}
    assert prov['inputs']['0011']['sha256'] == sha_of(tmp / 'data' / 'synthetic_11.fits')
    assert load(out / 'summary.json')['provenance']['code_hash'] == prov['code_hash']
    assert 'not confirmed planets' in (out / 'REPORT.md').read_text()
    groups = pd.read_csv(out / 'reports' / 'ranked-lead-groups.csv')
    assert len(groups) > 0 and set(groups.holdout_status) == {'not_applicable_holdouts_disabled'}
    before = {q: q.read_bytes() for q in out.glob('sources/*/*.json')}
    p = subprocess.run(cmd + ['run', str(manifest)], capture_output=True, text=True, timeout=120)
    assert p.returncode == 0, p.stderr
    assert {q: q.read_bytes() for q in out.glob('sources/*/*.json')} == before
    s = subprocess.run(cmd + ['status', str(manifest)], capture_output=True, text=True, timeout=60)
    assert json.loads(s.stdout)['terminal_sources'] == len(ranks)
    m = json.loads(manifest.read_text()); m['seed_base'] += 1; manifest.write_text(json.dumps(m))
    p = subprocess.run(cmd + ['run', str(manifest)], capture_output=True, text=True, timeout=60)
    assert p.returncode != 0 and 'seed_base' in p.stderr
    p = subprocess.run(cmd + ['verify-pins'], capture_output=True, text=True, timeout=60)
    assert p.returncode == 0 and json.loads(p.stdout)['ok'] is True


def test_second_coordinator_on_same_output_is_refused(tmp_path):
    with coordinator.Locks(tmp_path / 'out'):
        assert not coordinator.lock_free(tmp_path / 'out')
        with pytest.raises(coordinator.AlreadyRunning):
            with coordinator.Locks(tmp_path / 'out'):
                pass
    assert coordinator.lock_free(tmp_path / 'out')


def test_cli_errors_are_clear(workspace, capsys):
    tmp, ranks = workspace
    m = write_manifest(tmp / 'm.json', tmp, arms=['bogus'])
    assert cli.main(['validate', str(m)]) == 1
    assert 'unknown arm' in capsys.readouterr().err


# --- the real frozen search on a synthetic injected transit (a few seconds) ----------------------------------

def test_frozen_search_recovers_injected_period(tmp_path):
    from deimos import search
    from deimos.util import load as jl
    cfg = jl(ROOT / 'deimos' / 'science-config.json'); cfg['search'] = dict(cfg['search'], max_period_days=6)
    d = synthetic_discovery(5, quarters=(1, 2, 3, 4))
    from deimos.core import physical_transit
    d['flux'] *= physical_transit(d.time.to_numpy(), 3.3, 101.0, 0.15, 0.004)
    packet = search.frozen_select(d, cfg)
    top = packet['peaks'][0]
    assert abs(top['period_days'] / 3.3 - 1) < 0.01 and top['snr'] > 7
    again = search.frozen_select(d, cfg)
    strip = lambda p: {k: v for k, v in p.items() if k != 'seconds'}
    assert strip(again) == strip(packet)                              # deterministic
