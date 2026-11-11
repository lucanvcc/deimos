import json, shutil

import pytest

from deimos import util

# Digests of the pinned numerical files.
ORIGINAL = {
    'adapter.py': '6f140bc7d3004184b3327fc36f0e856feac9edbe3b07e74b170339cc074307ec',
    'core.py': 'aec735d504e32862ce9ea352a050af5fe6392aa02fa0b15ed8f3be0dc4562e06',
}


def test_pin_file_matches_digests():
    assert util.load(util.PINS)['files'] == ORIGINAL


def test_pinned_files_match_pins():
    util.check_pins()
    for rel, expected in ORIGINAL.items():
        assert util.sha(util.PACKAGE / rel) == expected, rel


def test_pin_drift_is_refused(tmp_path, monkeypatch):
    pins = tmp_path / 'pins.json'
    bad = dict(ORIGINAL, **{'adapter.py': '0' * 64})
    pins.write_text(json.dumps({'files': bad}))
    monkeypatch.setattr(util, 'PINS', pins)
    with pytest.raises(ValueError, match='Frozen input drift: adapter.py'):
        util.check_pins()


def test_edited_core_is_refused(tmp_path, monkeypatch):
    copy = tmp_path / 'pkg'
    shutil.copytree(util.PACKAGE, copy, ignore=shutil.ignore_patterns('__pycache__'))
    (copy / 'core.py').write_text((copy / 'core.py').read_text() + '\n# edited\n')
    monkeypatch.setattr(util, 'PACKAGE', copy)
    monkeypatch.setattr(util, 'PINS', copy / 'frozen-pins.json')
    with pytest.raises(ValueError, match='Frozen input drift: core.py'):
        util.check_pins()


def test_frozen_modules_are_the_only_numerical_imports():
    import deimos.search as s, deimos.science as sc
    assert s.extract_multi.__module__ == 'deimos.core'
    assert sc.physical_transit.__module__ == 'deimos.core'
