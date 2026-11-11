"""Holdout protection, exercised only with synthetic fixtures (no real sealed data exists in this repository)."""
import json

import pandas as pd
import pytest

from deimos import stages, coordinator
from deimos.manifest import Batch
from deimos.source import Context, process
from deimos.util import load, write, sha
from tests.conftest import synthetic_discovery, write_manifest


def _batch(ws, out, state):
    tmp, ranks = ws
    b = Batch(write_manifest(tmp / f'{out}.json', tmp, out=out, holdouts=True, state_dir=str(state), stages=[3]))
    b.validate()
    return b


def _sealed_fixture(ws, out, state):
    batch = _batch(ws, out, state)
    folder = batch.output_dir / 'sources' / '0011'; folder.mkdir(parents=True)
    synthetic_discovery(99, quarters=(4,)).to_csv(folder / 'sealed-holdout.csv.gz', index=False)
    write(folder / 'split.json', dict(holdout_sha256=sha(folder / 'sealed-holdout.csv.gz'), input_sha256='f' * 64))
    lead = dict(rank=11, source='Gaia DR3 11', group=1, peak_index=0,
                peak=dict(period_days=3.0, epoch_bkjd=100.5, duration_days=0.1, depth=0.001, snr=9.0, initial_candidate=True, diagnostics=dict(photometry_vetting_pass=True)))
    write(stages.stage_path(batch.output_dir, 3), dict(stage=3, leads=[lead]))
    return batch, folder, load(stages.stage_path(batch.output_dir, 3))['leads'][0]


def test_holdout_once_only_and_ledger_blocks_second_batch(workspace, tmp_path):
    state = tmp_path / 'shared-ledger'
    batch, folder, lead = _sealed_fixture(workspace, 'one', state)
    first = stages.holdout(batch, lead, 3)
    assert first['refit'] is False and (folder / 'access-0.json').exists()
    assert json.loads((folder / 'access-0.json').read_text())['once_only'] is True
    assert stages.holdout(batch, lead, 3) == first                       # result is reused, not recomputed
    (folder / 'holdout-0.json').unlink()
    assert stages.holdout(batch, lead, 3)['status'] == 'interrupted_holdout_unresolved'   # access without a result is never retried
    other, ofolder, olead = _sealed_fixture(workspace, 'two', state)
    assert stages.holdout(other, olead, 3)['status'] == 'holdout_previously_accessed_by_other_batch'
    assert not (ofolder / 'access-0.json').exists()                      # the second batch never opened the sealed file
    assert (state / 'holdout-ledger' / ('f' * 64 + '.json')).exists()


def test_separate_state_dirs_do_not_share_the_ledger(workspace, tmp_path):
    """Documented limitation: protection across batches exists only when they share one state_dir."""
    a, fa, la = _sealed_fixture(workspace, 'a', tmp_path / 'ledger-a')
    b, fb, lb = _sealed_fixture(workspace, 'b', tmp_path / 'ledger-b')
    assert stages.holdout(a, la, 3)['status'] != 'holdout_previously_accessed_by_other_batch'
    assert stages.holdout(b, lb, 3)['status'] != 'holdout_previously_accessed_by_other_batch'


def test_holdout_rejects_lead_not_in_freeze(workspace, tmp_path):
    batch, folder, lead = _sealed_fixture(workspace, 'three', tmp_path / 'ledger')
    with pytest.raises(ValueError, match='absent'):
        stages.holdout(batch, dict(lead, peak_index=5), 3)
    assert not (folder / 'access-5.json').exists()


def test_holdout_detects_sealed_file_drift_before_any_access(workspace, tmp_path):
    batch, folder, lead = _sealed_fixture(workspace, 'drift', tmp_path / 'ledger')
    synthetic_discovery(100, quarters=(4,)).to_csv(folder / 'sealed-holdout.csv.gz', index=False)
    with pytest.raises(ValueError, match='Sealed holdout drift'):
        stages.holdout(batch, lead, 3)
    assert not (folder / 'access-0.json').exists()
    assert not (tmp_path / 'ledger' / 'holdout-ledger').exists()


def test_holdout_requires_a_sealed_file_record(workspace, tmp_path):
    batch, folder, lead = _sealed_fixture(workspace, 'nosplit', tmp_path / 'ledger')
    write(folder / 'split.json', dict(input_sha256='f' * 64))
    with pytest.raises(ValueError, match='No sealed holdout'):
        stages.holdout(batch, lead, 3)


def test_discovery_only_batch_never_creates_a_holdout(workspace):
    tmp, ranks = workspace
    batch = Batch(write_manifest(tmp / 'm.json', tmp)); batch.validate()
    ctx = Context(batch)
    for row in batch.selection():
        assert process(ctx, row)['status'] == 'complete'
    assert not list((tmp / 'out').glob('sources/*/sealed-holdout*'))
    assert 'holdout_sha256' not in load(tmp / 'out' / 'sources' / '0011' / 'split.json')
    assert not (tmp / 'out' / 'state').exists()


def test_search_never_reads_the_sealed_file(workspace, tmp_path, monkeypatch):
    """With holdouts on, a source run writes the sealed file but no arm may read it."""
    tmp, ranks = workspace
    batch = Batch(write_manifest(tmp / 'm.json', tmp, holdouts=True, state_dir=str(tmp_path / 'ledger'))); batch.validate()
    opened = []
    real_read_csv = pd.read_csv
    monkeypatch.setattr(pd, 'read_csv', lambda p, *a, **k: (opened.append(str(p)), real_read_csv(p, *a, **k))[1])
    ctx = Context(batch)
    process(ctx, batch.selection()[0])
    sealed = list((tmp / 'out' / 'sources').glob('*/sealed-holdout.csv.gz'))
    assert len(sealed) == 1 and not [p for p in opened if 'holdout' in p]
    assert load(sealed[0].parent / 'split.json')['holdout_sha256'] == sha(sealed[0])


def test_end_to_end_holdout_run_writes_ledger_and_blocks_replication(workspace, tmp_path):
    """Full coordinator with holdouts on synthetic FITS: stage freeze -> one holdout access per lead -> second batch is refused."""
    tmp, ranks = workspace
    state = tmp_path / 'ledger'
    # Make every real-arm peak a lead (fake plugin: photometry_vetting_pass = snr > 7) by using a low-noise fixture where possible;
    # either way the invariants below must hold for whatever leads exist.
    m1 = write_manifest(tmp / 'm1.json', tmp, out='h1', holdouts=True, state_dir=str(state), workers=1)
    coordinator.run(Batch(m1))
    out1 = tmp / 'h1'
    assert load(out1 / 'terminal.json')['status'] == 'COMPLETE'
    frozen = load(out1 / 'stage-3-freeze.json')
    assert frozen['leads'], 'fixture must produce at least one lead so the holdout path is exercised'
    accessed = sorted(out1.glob('sources/*/access-*.json'))
    assert len(accessed) == len(frozen['leads'])
    for a in accessed:
        assert (a.parent / a.name.replace('access', 'holdout')).exists()
    m2 = write_manifest(tmp / 'm2.json', tmp, out='h2', holdouts=True, state_dir=str(state), workers=1)
    coordinator.run(Batch(m2))
    out2 = tmp / 'h2'
    assert not list(out2.glob('sources/*/access-*.json'))
    statuses = {load(p)['status'] for p in out2.glob('sources/*/holdout-*.json')}
    assert statuses <= {'holdout_previously_accessed_by_other_batch'}
    assert len(list(out2.glob('sources/*/holdout-*.json'))) == len(frozen['leads'])
