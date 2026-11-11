import pytest
from deimos.pixel import host_A, read_leads, verdict_of

T, N = 'Gaia DR3 1', 'Gaia DR3 2'
q = lambda *hosts: [{'decision': {'host': h}} for h in hosts]


def test_leads_csv(tmp_path):
    f = tmp_path / 'l.csv'
    f.write_text('rank,lg,role,gaia,origin,peak_index,src\n25901,62,primary,2127936912539397632,10528990,0,x/y\n')
    lead, = read_leads(f)
    assert lead['rank'] == 25901 and lead['gaia'] == '2127936912539397632' and lead['src'].name == 'y'
    f.write_text('rank,lg\n1,2\n')
    with pytest.raises(ValueError): read_leads(f)


def test_decision_rule():
    assert host_A(q(T, None, T)) == T
    assert host_A(q(T, N)) is None  # a quarter assigning another star blocks A
    assert host_A(q(None)) is None
    assert verdict_of(T, T, T) == 'KEEP'
    assert verdict_of(N, N, T) == 'DROP'
    assert verdict_of(T, N, T) == verdict_of(None, None, T) == verdict_of(N, None, T) == 'UNRESOLVED'
