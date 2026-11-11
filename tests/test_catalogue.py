"""Unit tests for deimos/catalogue.py (no network)."""
from fractions import Fraction

import pytest

from deimos import catalogue as cat


def test_relation_and_phase_residual():
    q, dev = cat.relation(3 * 2.5000, 2.5001)
    assert q == Fraction(3) and abs(dev) < 1e-4
    # our P = 1/3 of the catalogue P: coincidences every P_ref/3
    assert cat.phase_residual(10.0 + 3 * 1.0, 3.0, 10.0, Fraction(1, 3)) == pytest.approx(0.0)
    assert cat.phase_residual(10.5, 3.0, 10.0, Fraction(1)) == pytest.approx(0.5)
    assert cat.phase_tol(0.1, 0.0) == pytest.approx(0.075)
    assert cat.phase_tol(0.2, 1e-4) == pytest.approx(0.15)  # drift never expands identity tolerance


def test_dikco_verdict_neighbour_and_same_star():
    dec = 45.0
    # catalogued star 5" north of the target; DV offset points 5" south (back onto the target)
    row = {'ra': '300.0', 'dec': str(dec + 5 / 3600), 'koi_dikco_mra': '0.0', 'koi_dikco_mdec': '-5.0', 'koi_dikco_msky_err': '0.1'}
    assert cat.dikco_verdict(row, 'koi_', 300.0, dec)['dikco_verdict'] == 'source at target'
    row2 = dict(row, koi_dikco_mdec='0.05')
    assert cat.dikco_verdict(row2, 'koi_', 300.0, dec)['dikco_verdict'] == 'source at catalogued star'
    same = {'ra': '300.0', 'dec': str(dec), 'koi_dikco_mra': '0.2', 'koi_dikco_mdec': '0.1', 'koi_dikco_msky_err': '0.1'}
    assert cat.dikco_verdict(same, 'koi_', 300.0, dec)['dikco_verdict'] == 'source at target'
    far = dict(same, koi_dikco_mdec='4.0')
    assert cat.dikco_verdict(far, 'koi_', 300.0, dec)['dikco_verdict'] == 'source offset from target'
    assert cat.dikco_verdict(dict(same, koi_dikco_mra=''), 'koi_', 300.0, dec)['dikco_verdict'] == 'no dikco offset'


def _member(lg, rank, P, T, snr):
    return {'lead_group': lg, 'dependency_group': 1, 'members': [{
        'rank': rank, 'source': f'Gaia DR3 {rank}', 'peak_index': lg, 'holdout': {'status': cat.SUPPORTED},
        'peak': {'period_days': P, 'epoch_bkjd': T, 'duration_days': 0.1, 'snr': snr,
                 'event_timing': {'positive_indices': list(range(int((1100-T)/P)+1)),
                                  'start': T, 'end': T + int((1100-T)/P)*P}}}]}


def test_families_link_aligned_harmonics_only():
    groups = [_member(1, 7, 2.5, 100.0, 50), _member(2, 7, 7.5, 102.5, 30),   # 3x, epoch on a parent transit
              _member(3, 7, 7.5, 101.25, 20),                                 # 3x, epoch half a cycle off
              _member(4, 8, 2.5, 100.0, 40)]                                  # other source
    fams = cat.families(groups)
    assert [(f['rank'], f['parent_lead_group'], [a['lead_group'] for a in f['aliases']]) for f in fams] == [
        (7, 1, [2]), (7, 3, []), (8, 4, [])]


def test_proposed_check_menu():
    assert cat.proposed_check({'catalogue_match': 'strict', 'dikco_verdict': 'source at catalogued star'}) == 'resolve the catalogue host'
    assert cat.proposed_check({'catalogue_match': 'strict_phase_offset', 'dikco_verdict': 'source at target'}) == 'check alias/odd-even'
    assert cat.proposed_check({'catalogue_match': 'none', 'odd_even_sigma': 2.6}) == 'check alias/odd-even'
    assert cat.proposed_check({'catalogue_match': 'none', 'alias_relations': 'lg5=5/3'}) == 'check alias/odd-even'
    assert cat.proposed_check({'catalogue_match': 'none', 'alias_relations': 'lg5=3;lg6=1/3', 'holdout_depth_ratio': 0.5}) == 'reconcile depth'
    assert cat.proposed_check({'catalogue_match': 'none', 'holdout_depth_ratio': 1.0}) == 'run the pixel host test'
    assert all(cat.proposed_check(r) in cat.NEXT_CHECKS for r in ({}, {'catalogue_match': 'loose'}))


def test_output_reader_refuses_holdout_and_split_files(tmp_path):
    for name in ('sealed-holdout.csv.gz', 'access-0.json', 'split.json', 'holdout-0.json'):
        with pytest.raises(ValueError):
            cat._load_output_json(tmp_path, name)


def _peaks_with_band():
    """Inverted-arm vetted peaks: a smooth log-uniform background plus a 15.3 d pile-up sharing one epoch."""
    peaks, rank = [], 0
    for i in range(400):                          # background: 1-30 d, epochs spread, one per source
        rank += 1
        P_ = 10 ** (0.0 + 1.477 * (i + 0.5) / 400)
        peaks.append({'arm': 'inverted', 'vetting_pass': 'True', 'period_days': P_, 'epoch_bkjd': 131 + (i * 7.3) % P_,
                      'rank': rank})
    for i in range(30):                           # systematic: 15.25-15.35 d, epochs at BKJD 139 +- 0.2
        rank += 1
        peaks.append({'arm': 'inverted', 'vetting_pass': 'True', 'period_days': 15.25 + 0.1 * i / 30,
                      'epoch_bkjd': 139.0 + 0.2 * ((i % 5) - 2) / 2, 'rank': rank})
    # real-arm and failing peaks must be ignored by the calibration
    peaks += [{'arm': 'real', 'vetting_pass': 'True', 'period_days': 3.0, 'epoch_bkjd': 132.0, 'rank': 9000 + k}
              for k in range(50)]
    peaks += [{'arm': 'inverted', 'vetting_pass': 'False', 'period_days': 5.0, 'epoch_bkjd': 132.0, 'rank': 9500 + k}
              for k in range(50)]
    return peaks


def test_systematic_band_is_calibrated_on_vetted_inverted_peaks_only():
    bands = cat.systematic_bands(_peaks_with_band())
    assert bands and all(14.5 < b['P_lo'] < 15.6 for b in bands.values())     # only the pile-up; not 3 d or 5 d
    index = cat.ephemeris_index(_peaks_with_band())
    hit = cat.systematic_flag({'period_days': '15.3', 'epoch_bkjd': '139.05', 'rank': '1'}, bands, index)
    assert 'period band' in hit and 'shared ephemeris' in hit
    off_epoch = cat.systematic_flag({'period_days': '15.3', 'epoch_bkjd': '146.0', 'rank': '1'}, bands, index)
    assert 'period band' in off_epoch and 'shared ephemeris' not in off_epoch
    assert cat.systematic_flag({'period_days': '8.0', 'epoch_bkjd': '135.0', 'rank': '1'}, bands, index) == ''
    assert 'check the shared systematic' in cat.NEXT_CHECKS


def test_shared_ephemeris_counts_other_sources_only():
    index = cat.ephemeris_index(_peaks_with_band())
    own = index[-1][2]
    n_all = len(cat.shared_ephemeris_sources(15.3, 139.0, 'none', index))
    n_ex = len(cat.shared_ephemeris_sources(15.3, 139.0, own, index))
    assert n_all >= cat.SYS_EPH_MIN_SOURCES and n_ex <= n_all


def test_catalogue_response_cache_skips_network(tmp_path, monkeypatch):
    import json
    monkeypatch.setattr(cat, 'cone_queries', lambda *a: {'test': 'https://example.invalid/fixed'})
    raw = tmp_path / 'raw'; raw.mkdir(); (raw / '1-test.txt').write_text('cached')
    log = tmp_path / 'query-log.jsonl'
    log.write_text(json.dumps(dict(http=200, file='1-test.txt', url='https://example.invalid/fixed')) + '\n')
    monkeypatch.setattr(cat.urllib.request, 'urlopen', lambda *a, **k: pytest.fail('cache must avoid network'))
    assert cat.fetch_catalogues({'1': (1, 2)}, raw, log) == []


def test_noninteger_resonance_is_not_collapsed():
    groups = [_member(1, 7, 3, 100, 50), _member(2, 7, 5, 100, 30)]
    assert len(cat.families(groups)) == 2


def timing_peak(period=2., epoch=100., duration=.1, indices=(0, 1, 2, 3), end=None, **extra):
    return dict(period_days=period, epoch_bkjd=epoch, duration_days=duration,
                event_timing=dict(positive_indices=list(indices), start=epoch,
                                  end=end if end is not None else epoch + max(indices)*period), **extra)


def test_distinct_commensurate_signals_have_unexplained_dips():
    parent = timing_peak(period=3, indices=range(101))
    child = timing_peak(period=5, indices=range(61))
    proof = cat.signal_explained_by(child, parent)
    assert proof['relation'] == '5/3'
    assert not proof['demonstrated'] and proof['explained_events'] < proof['events']


def test_true_integer_and_rational_aliases_require_measured_event_coverage():
    parent = timing_peak(period=2, indices=range(101))
    child = timing_peak(period=6, indices=range(34))
    assert cat.signal_explained_by(child, parent)['demonstrated']
    # At a shorter trial period, the extra trial events are null; positive events coincide.
    parent = timing_peak(period=6, indices=range(34))
    child = timing_peak(period=2, indices=range(0, 100, 3))
    assert cat.signal_explained_by(child, parent)['demonstrated']
    # A 5/3 relationship is allowed only where actual positive events establish it.
    parent = timing_peak(period=3, indices=range(101))
    child = timing_peak(period=5, indices=range(0, 61, 3))
    assert cat.signal_explained_by(child, parent)['demonstrated']


def test_period_drift_is_tested_over_observed_baseline_not_added_to_tolerance():
    parent = timing_peak(period=2, indices=range(701))
    child = timing_peak(period=2.001, indices=(0, 1, 2), end=1500.)
    proof = cat.signal_explained_by(child, parent)
    assert proof['explained_events'] == proof['events']
    assert proof['accumulated_drift_days'] > proof['tolerance_days']
    assert not proof['demonstrated']
    assert not cat.same_ephemeris(child, parent)


@pytest.mark.parametrize('missing', ['epoch_bkjd', 'duration_days', 'event_timing'])
def test_missing_identity_evidence_retains_signal(missing):
    parent = timing_peak()
    child = timing_peak()
    child.pop(missing)
    assert not cat.signal_explained_by(child, parent)['demonstrated']


def test_same_period_separate_eclipses_are_not_duplicates():
    first = timing_peak(period=10, epoch=100)
    second = timing_peak(period=10, epoch=102)
    assert not cat.same_ephemeris(first, second)
    assert not cat.signal_explained_by(second, first)['demonstrated']
    assert not cat.signal_explained_by(first, second)['demonstrated']


def test_measured_secondary_can_explain_alias_but_brightening_cannot():
    parent = timing_peak(period=4, indices=range(10), secondary_snr=20)
    child = timing_peak(period=2, indices=range(18))
    assert cat.signal_explained_by(child, parent)['demonstrated']
    parent['secondary_snr'] = -20
    assert not cat.signal_explained_by(child, parent)['demonstrated']


def test_catalogue_primary_must_explain_events_and_have_reliable_source():
    child = timing_peak(period=5, indices=range(20))
    parent = dict(period_days=5., epoch_bkjd=100., duration_days=.1)
    proof = cat.signal_explained_by(child, parent, catalogue=True)
    assert proof['demonstrated']
    wrong = dict(parent, epoch_bkjd=101.)
    assert not cat.signal_explained_by(child, wrong, catalogue=True)['demonstrated']


def test_supported_groups_and_sources_read_only_the_reports(tmp_path):
    import json
    (tmp_path / 'reports').mkdir()
    (tmp_path / 'sources/0001').mkdir(parents=True)
    (tmp_path / 'sources/0001/receipt.json').write_text(json.dumps({'rank': 1, 'status': 'complete', 'ra': 1, 'dec': 2}))
    peak = {'period_days': 2., 'epoch_bkjd': 1., 'duration_days': .1, 'depth': 1e-3, 'diagnostics': {}}
    groups = [{'lead_group': g, 'dependency_group': 1, 'members': [{
        'rank': 1, 'source': 'Gaia DR3 123456', 'group': 1, 'peak_index': g, 'peak': peak,
        'holdout': {'status': st, 'evaluation': {'depth': 5e-4}}}]} for g, st in ((1, cat.SUPPORTED), (2, 'failed'))]
    (tmp_path / 'reports/ranked-lead-groups.jsonl').write_text('\n'.join(json.dumps(g) for g in groups))
    got = cat.supported_groups(tmp_path)
    assert [g['lead_group'] for g in got] == [1]
    assert got[0]['members'][0]['holdout']['depth_ppm'] == pytest.approx(500)
    assert list(cat.source_info(tmp_path)) == [1]
