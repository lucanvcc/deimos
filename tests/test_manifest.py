import json

import pandas as pd
import pytest

from deimos.manifest import Batch, ManifestError
from tests.conftest import write_manifest, REG_COLUMNS


def load(ws, **extra):
    tmp, ranks = ws
    return Batch(write_manifest(tmp / 'm.json', tmp, **extra))


def test_valid_manifest_defaults(workspace):
    b = load(workspace)
    rows = b.validate()
    assert [r['rank'] for r in rows] == [11, 12, 13]
    assert b.stages == [3] and b.holdout_stage == 3 and b.holdouts is False
    assert len(b.config_hash()) == 64


@pytest.mark.parametrize('extra,match', [
    (dict(arms=['real', 'shuffled']), 'unknown arm'),
    (dict(arms=['real', 'real']), 'duplicate arm'),
    (dict(arms=[{'name': 'real', 'fraction': 0.5}]), 'takes only a name'),
    (dict(arms=[]), 'non-empty'),
    (dict(seed_base='abc'), 'seed_base must be an integer'),
    (dict(seed_base=1.5), 'seed_base must be an integer'),
    (dict(workers=0), 'workers'),
    (dict(niceness=5), 'unknown manifest fields'),
    (dict(search_plugin='nocolon'), 'module:function'),
    (dict(budget={'max_seconds': 1}), 'unknown manifest fields'),
    (dict(product_dirs=['x']), 'unknown manifest fields'),
    (dict(search_strategy='coarse_to_fine'), 'unknown manifest fields'),
    (dict(discovery_from_run='x'), 'unknown manifest fields'),
    (dict(sead_base=1), 'unknown manifest fields'),
    (dict(select=dict(top=3)), 'unknown select keys'),
    (dict(holdouts=True), 'explicit state_dir'),
])
def test_manifest_rejections(workspace, extra, match):
    with pytest.raises(ManifestError, match=match):
        load(workspace, **extra).validate()


def test_stage_selection_mismatch(workspace):
    with pytest.raises(ManifestError, match='final stage'):
        load(workspace, stages=[5]).validate()
    with pytest.raises(ManifestError, match='strictly increasing'):
        load(workspace, stages=[2, 2])
    with pytest.raises(ManifestError, match='holdout_stage'):
        load(workspace, stages=[3], holdout_stage=2).validate()
    b = load(workspace, select=dict(min_rank=11, exclude_ranks=[12]), stages=[2])
    assert [r['rank'] for r in b.selection()] == [11, 13]
    with pytest.raises(ManifestError, match='not in the registration'):
        load(workspace, select=dict(ranks=[11, 99]), stages=[2]).validate()


def test_holdouts_with_explicit_state_dir_is_accepted(workspace):
    tmp, _ = workspace
    b = load(workspace, holdouts=True, state_dir=str(tmp / 'ledger'))
    b.validate()
    assert b.state_dir == (tmp / 'ledger').resolve()


def test_state_dir_default_is_under_output_dir(workspace):
    tmp, _ = workspace
    assert load(workspace).state_dir == (tmp / 'out' / 'state').resolve()


def test_missing_and_invalid_files(tmp_path, workspace):
    with pytest.raises(ManifestError, match='manifest not found'):
        Batch(tmp_path / 'nope.json')
    bad = tmp_path / 'bad.json'; bad.write_text('{not json')
    with pytest.raises(ManifestError, match='not valid JSON'):
        Batch(bad)
    arr = tmp_path / 'arr.json'; arr.write_text('[]')
    with pytest.raises(ManifestError, match='JSON object'):
        Batch(arr)
    miss = tmp_path / 'miss.json'; miss.write_text(json.dumps(dict(output_dir='o')))
    with pytest.raises(ManifestError, match='missing fields'):
        Batch(miss)
    with pytest.raises(ManifestError, match='registration_csv not found'):
        Batch(write_manifest(tmp_path / 'm2.json', tmp_path, registration_csv=str(tmp_path / 'absent.csv'))).validate()


def mutate_registration(ws, fn):
    tmp, _ = ws
    reg = pd.read_csv(tmp / 'registration.csv'); fn(reg); reg.to_csv(tmp / 'registration.csv', index=False)


@pytest.mark.parametrize('fn,match', [
    (lambda r: r.drop(columns=['kic'], inplace=True), 'missing columns'),
    (lambda r: r.__setitem__('rank', [1, 1, 2]), 'unique positive'),
    (lambda r: r.__setitem__('rank', [0, 1, 2]), 'unique positive'),
    (lambda r: r.__setitem__('dec', [95.0, 1.0, 1.0]), 'out of range'),
    (lambda r: r.__setitem__('file_name', ['../escape.fits', 'a.fits', 'b.fits']), 'bare \\*.fits'),
    (lambda r: r.__setitem__('file_name', ['a.txt', 'a.fits', 'b.fits']), 'bare \\*.fits'),
    (lambda r: r.__setitem__('sha256', ['xyz', 'a' * 64, 'b' * 64]), '64 lowercase hex'),
    (lambda r: r.drop(r.index, inplace=True), 'no rows'),
])
def test_registration_rejections(workspace, fn, match):
    mutate_registration(workspace, fn)
    with pytest.raises(ManifestError, match=match):
        load(workspace).validate()


def test_science_config_must_have_sections(workspace):
    tmp, _ = workspace
    (tmp / 'cfg.json').write_text(json.dumps({'clean': {}, 'search': {}}))
    with pytest.raises(ManifestError, match='"validate" section'):
        load(workspace, science_config=str(tmp / 'cfg.json')).validate()


def test_config_hash_tracks_scientific_fields(workspace):
    base = load(workspace).config_hash()
    assert load(workspace, injection={'depths': [0.1]}).config_hash() != base
    assert load(workspace, search_plugin='tests.fakeplugin:select2').config_hash() != base
    assert load(workspace, workers=5).config_hash() == base


def test_registration_columns_constant():
    assert REG_COLUMNS[:9] == ['rank', 'gaia_designation', 'conservative_group', 'tpforg', 'kic', 'ra', 'dec', 'file_name', 'url']
