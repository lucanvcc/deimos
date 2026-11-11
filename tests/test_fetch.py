"""Fetch (library and CLI), exercised offline through file:// URLs (no network)."""
import csv, json


from deimos import cli, inputs
from tests.conftest import ROOT, make_fits, sha_of, write_manifest


def upstream(tmp_path, ranks=(1, 2)):
    (tmp_path / 'upstream').mkdir()
    rows = []
    for r in ranks:
        src = make_fits(tmp_path / 'upstream' / f'f{r}.fits', gaia=f'Gaia DR3 {r}', seed=r)
        rows.append(dict(rank=r, gaia_designation=f'Gaia DR3 {r}', conservative_group=r, tpforg=r, kic=r, ra=1.0, dec=2.0,
                         file_name=src.name, url=src.as_uri(), sha256=sha_of(src), bytes=src.stat().st_size))
    return rows


def write_reg(path, rows):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def manifest(tmp_path, rows, **extra):
    write_reg(tmp_path / 'registration.csv', rows)
    return write_manifest(tmp_path / 'm.json', tmp_path, input_dirs=[str(tmp_path / 'dest')], **extra)


def test_cli_fetch_verifies_skips_and_honours_select(tmp_path, capsys):
    rows = upstream(tmp_path)
    m = manifest(tmp_path, rows, select=dict(ranks=[2]), stages=[1])
    assert cli.main(['fetch', str(m), '--dry-run']) == 0
    assert not (tmp_path / 'dest').exists() and 'f2.fits' in capsys.readouterr().out
    assert cli.main(['fetch', str(m)]) == 0
    assert json.loads(capsys.readouterr().out)['fetched'] == ['f2.fits']
    assert [p.name for p in (tmp_path / 'dest').iterdir()] == ['f2.fits']       # rank 1 not selected
    assert cli.main(['fetch', str(m)]) == 0
    assert json.loads(capsys.readouterr().out)['skipped'] == ['f2.fits']
    assert cli.main(['validate', str(m)]) == 0


def test_cli_fetch_exit_1_and_cleanup_on_mismatch(tmp_path, capsys):
    rows = upstream(tmp_path); rows[0]['sha256'] = '0' * 64
    m = manifest(tmp_path, rows)
    assert cli.main(['fetch', str(m)]) == 1
    out = json.loads(capsys.readouterr().out)
    assert out['fetched'] == ['f2.fits'] and out['failed'][0]['file'] == 'f1.fits' and 'does not match' in out['failed'][0]['error']
    names = sorted(p.name for p in (tmp_path / 'dest').iterdir())
    assert names == ['f2.fits']                                                  # no .part, nothing unverified


def test_fetch_size_mismatch_and_bad_scheme_and_unreachable(tmp_path):
    rows = upstream(tmp_path)
    r = dict(rows[0], bytes=5)
    res = inputs.fetch([r], tmp_path / 'd1')
    assert res['failed'] and 'size' in res['failed'][0]['error'] and not list((tmp_path / 'd1').iterdir())
    res = inputs.fetch([dict(rows[0], url='ftp://example.invalid/x.fits')], tmp_path / 'd2')
    assert 'unsupported URL scheme' in res['failed'][0]['error']
    res = inputs.fetch([dict(rows[0], url=(tmp_path / 'nope.fits').as_uri())], tmp_path / 'd3')
    assert res['failed'] and not list((tmp_path / 'd3').iterdir())


def test_mirror_fallback_urls():
    u = 'https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:HLSP/kbonus-bkg/lcs/0113/a.fits'
    assert inputs.candidate_urls(u) == [u, 'https://archive.stsci.edu/hlsps/kbonus-bkg/lcs/0113/a.fits']
    assert inputs.candidate_urls('file:///x.fits') == ['file:///x.fits']


def test_fetch_falls_back_to_second_url(tmp_path, monkeypatch):
    row = upstream(tmp_path)[0]
    good = row['url']
    monkeypatch.setattr(inputs, 'candidate_urls', lambda url: [(tmp_path / 'missing.fits').as_uri(), good])
    assert inputs.fetch([row], tmp_path / 'd')['fetched'] == ['f1.fits']


def test_run_exits_nonzero_when_a_source_fails(workspace, capsys):
    tmp, ranks = workspace
    m = write_manifest(tmp / 'm.json', tmp, workers=1)
    (tmp / 'data' / 'synthetic_12.fits').write_bytes(b'tampered')
    assert cli.main(['run', str(m)]) == 1
    assert 'did not complete' in capsys.readouterr().err


def test_packaged_example_registration_is_pinned():
    rows = list(csv.DictReader(open(ROOT / 'examples' / 'kepler-minimal' / 'registration.csv')))
    assert rows and all(len(r['sha256']) == 64 and r['url'].startswith('https://mast.stsci.edu/') for r in rows)
