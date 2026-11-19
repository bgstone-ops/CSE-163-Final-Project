"""
Benjamin Stone
CSE 163 AD
Tests the functions in eda.py on small tables where the correct
answers are known.
"""
import pandas as pd

import eda


def make_test_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Returns small drugs, links and trials DataFrames. Drug A has two
    trials (2005 and 2010) and is approved. Drug B's only trial was
    terminated. Drug C's only trial is still recruiting.
    """
    drugs = pd.DataFrame({'molecule_chembl_id': ['A', 'B', 'C'],
                          'max_phase': [4, 2, 2]})
    links = pd.DataFrame({'molecule_chembl_id': ['A', 'A', 'B', 'C'],
                          'nct_id': ['N1', 'N2', 'N3', 'N4']})
    trials = pd.DataFrame({
        'nct_id': ['N1', 'N2', 'N3', 'N4'],
        'start_year': [2010, 2005, 2012, 2015],
        'overall_status': ['COMPLETED', 'COMPLETED', 'TERMINATED',
                           'RECRUITING'],
        'lead_sponsor': ['Pfizer', 'Merck', 'UW', 'UW'],
        'sponsor_class': ['INDUSTRY', 'INDUSTRY', 'OTHER', 'OTHER']
    })
    return drugs, links, trials


def test_clean_company() -> None:
    """Tests the function clean_company."""
    assert eda.clean_company('Pfizer Inc.') == 'PFIZER'
    assert eda.clean_company('PFIZER INC') == 'PFIZER'
    assert eda.clean_company('Merck & Co., Inc.') == 'MERCK &'
    assert eda.clean_company(None) is None


def test_join_drugs_trials() -> None:
    """Tests the function join_drugs_trials."""
    drugs, links, trials = make_test_data()
    joined = eda.join_drugs_trials(drugs, links, trials)
    joined = joined.set_index('molecule_chembl_id')
    assert len(joined) == 3
    assert joined.loc['A', 'first_year'] == 2005
    assert joined.loc['A', 'last_year'] == 2010
    assert joined.loc['A', 'lead_sponsor'] == 'Merck'
    assert joined.loc['B', 'last_status'] == 'TERMINATED'


def test_add_outcome() -> None:
    """Tests the function add_outcome."""
    drugs, links, trials = make_test_data()
    joined = eda.join_drugs_trials(drugs, links, trials)
    result = eda.add_outcome(joined).set_index('molecule_chembl_id')
    assert sorted(result.index) == ['A', 'B']
    assert result.loc['A', 'outcome'] == 'approved'
    assert result.loc['B', 'outcome'] == 'failed'
    assert result.loc['A', 'approved']
    assert not result.loc['B', 'approved']


def test_smallest_expected() -> None:
    """
    Tests the function smallest_expected. Row totals are 10 and 30 and
    column totals are 20 and 20, so the smallest expected count is
    10 * 20 / 40 = 5.
    """
    df = pd.DataFrame({
        'group': ['x'] * 10 + ['y'] * 30,
        'outcome': ['approved'] * 5 + ['failed'] * 5 +
                   ['approved'] * 15 + ['failed'] * 15
    })
    assert eda.smallest_expected(df, 'group') == 5
    df = pd.DataFrame({
        'group': ['x'] * 4 + ['y'] * 36,
        'outcome': ['approved'] * 4 + ['approved'] * 4 + ['failed'] * 32
    })
    assert eda.smallest_expected(df, 'group') == 4 * 8 / 40


def main():
    test_clean_company()
    test_join_drugs_trials()
    test_add_outcome()
    test_smallest_expected()
    print('All tests passed!')


if __name__ == '__main__':
    main()