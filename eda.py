"""
Benjamin Stone
CSE 163 AD
EDA for final project.

Run download_data.py first. This prints a short summary of each
dataset and of the joined data, checks the assumptions of the
statistical tests, and saves 2 plots for each joined dataset (6 in
total) in the plots folder.
"""
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

COMPANY_WORDS = ['INC', 'CORP', 'CORPORATION', 'CO', 'COMPANY', 'LTD',
                 'LLC', 'PLC', 'AG', 'SA', 'THE']
SPONSOR_GROUPS = {'INDUSTRY': 'Industry', 'OTHER': 'Academic/hospital',
                  'NETWORK': 'Academic/hospital', 'NIH': 'Government',
                  'FED': 'Government', 'OTHER_GOV': 'Government'}


def summarize(name: str, df: pd.DataFrame, columns: list[str]) -> None:
    """
    Prints the given name, the number of rows and columns in df, and
    the missing values and a summary of each of the given columns.
    """
    print()
    print('==========', name, '==========')
    print('Rows:', len(df), ' Columns:', len(df.columns))
    print('Missing values:')
    print(df[columns].isna().sum())
    print(df[columns].describe(include='all').round(2))


def save(title: str, filename: str) -> None:
    """
    Adds the given title to the current plot and saves it in the plots
    folder under the given file name.
    """
    plt.title(title)
    plt.savefig('plots/' + filename, bbox_inches='tight')
    plt.close()


def clean_company(name: str | None) -> str | None:
    """
    Takes in a company name and returns it in capitals without periods,
    commas or words like INC, so 'Pfizer Inc.' and 'PFIZER INC' both
    become 'PFIZER'. Returns None if name is missing.
    """
    if not isinstance(name, str):
        return None
    words = name.upper().replace('.', '').replace(',', '').split()
    kept = []
    for word in words:
        if word not in COMPANY_WORDS:
            kept.append(word)
    return ' '.join(kept)


def join_drugs_trials(drugs: pd.DataFrame, links: pd.DataFrame,
                      trials: pd.DataFrame) -> pd.DataFrame:
    """
    Takes in the drugs, drug-trial links and trials DataFrames and
    returns one row per drug with a dated trial: the drug's columns
    plus the start year (first_year) and sponsor of its first trial and
    the start year (last_year) and status (last_status) of its last
    trial.
    """
    df = links.merge(trials, left_on='nct_id', right_on='nct_id')
    df = df.dropna(subset=['start_year']).sort_values('start_year')
    first = df.drop_duplicates(subset='molecule_chembl_id', keep='first')
    last = df.drop_duplicates(subset='molecule_chembl_id', keep='last')
    first = first[['molecule_chembl_id', 'start_year', 'lead_sponsor',
                   'sponsor_class']].rename(columns={'start_year':
                                                     'first_year'})
    last = last[['molecule_chembl_id', 'start_year', 'overall_status']]
    last = last.rename(columns={'start_year': 'last_year',
                                'overall_status': 'last_status'})
    joined = drugs.merge(first, on='molecule_chembl_id')
    return joined.merge(last, on='molecule_chembl_id')


def add_outcome(df: pd.DataFrame, current_year: int = 2026) -> pd.DataFrame:
    """
    Takes in a DataFrame from join_drugs_trials and returns a new
    DataFrame of the drugs whose first trial started from 2000 through
    2022, with an outcome column: 'approved' if max_phase is 4, and
    'failed' if the last trial was stopped, or started 3 or more years
    before current_year and is not still running. Drugs still in
    development are left out. Also adds an approved column that is True
    for approved drugs.
    """
    df = df[(df['first_year'] >= 2000) & (df['first_year'] <= 2022)].copy()
    stopped = df['last_status'].isin(['TERMINATED', 'WITHDRAWN',
                                      'SUSPENDED'])
    running = df['last_status'].isin(['RECRUITING', 'NOT_YET_RECRUITING',
                                      'ACTIVE_NOT_RECRUITING',
                                      'ENROLLING_BY_INVITATION'])
    old = df['last_year'] <= current_year - 3
    df['outcome'] = 'in development'
    df.loc[stopped | (old & ~running), 'outcome'] = 'failed'
    df.loc[df['max_phase'] == 4, 'outcome'] = 'approved'
    df = df[df['outcome'] != 'in development'].copy()
    df['approved'] = df['outcome'] == 'approved'
    return df


def smallest_expected(df: pd.DataFrame, column: str) -> float:
    """
    Takes in a DataFrame with an outcome column and the name of another
    column, and returns the smallest expected count in the chi-square
    table of that column against outcome. Each expected count is the
    row total times the column total divided by the overall total.
    """
    table = pd.crosstab(df[column], df['outcome'])
    row_totals = table.sum(axis=1)
    column_totals = table.sum(axis=0)
    return row_totals.min() * column_totals.min() / row_totals.sum()


def check_chi_square(df: pd.DataFrame, column: str) -> None:
    """
    Prints the table of the given column against outcome and its
    smallest expected count, which must be at least 5 for a chi-square
    test.
    """
    print()
    print('Chi-square check for', column)
    print(pd.crosstab(df[column], df['outcome']))
    smallest = smallest_expected(df, column)
    print('Smallest expected count:', round(smallest, 2),
          '(needs to be at least 5)')


def check_mann_whitney(df: pd.DataFrame, column: str) -> None:
    """
    Prints the number of drugs and the median of the given column for
    approved and failed drugs, to check the groups for a Mann-Whitney U
    test.
    """
    print()
    print('Mann-Whitney check for', column)
    print(df.groupby('outcome')[column].agg(['count', 'median']).round(2))


def main():
    sns.set()
    pd.set_option('display.width', 200)
    os.makedirs('plots', exist_ok=True)
    drugs = pd.read_csv('data/drugs.csv')
    links = pd.read_csv('data/links.csv')
    trials = pd.read_csv('data/trials.csv')
    potency = pd.read_csv('data/potency.csv')
    sec = pd.read_csv('data/sec.csv')
    orphan = pd.read_csv('data/orphan.csv')

    # Each dataset on its own (summaries only)
    summarize('ChEMBL drugs', drugs,
              ['molecule_type', 'max_phase', 'full_mwt', 'alogp', 'psa'])
    print(drugs['molecule_type'].value_counts())
    summarize('ClinicalTrials.gov trials', trials,
              ['start_year', 'overall_status', 'sponsor_class'])
    print(trials['sponsor_class'].value_counts())
    sec['rd_millions'] = sec['rd_expense'] / 1000000
    summarize('SEC R&D spending', sec, ['year', 'rd_millions'])
    summarize('FDA orphan designations', orphan,
              ['generic_name', 'designation_date', 'status'])

    # Joined: ChEMBL + ClinicalTrials.gov
    joined = add_outcome(join_drugs_trials(drugs, links, trials))
    best = potency.groupby('molecule_chembl_id')['pchembl_value'].max()
    joined = joined.merge(best.reset_index(), how='left')
    joined['sponsor_class'] = joined['sponsor_class'].map(SPONSOR_GROUPS)
    summarize('Joined: ChEMBL + ClinicalTrials.gov', joined,
              ['outcome', 'molecule_type', 'first_year', 'pchembl_value'])
    print(joined.groupby('molecule_type')['approved'].mean().round(2))
    sns.catplot(data=joined, y='molecule_type', x='approved', kind='bar',
                errorbar=None)
    plt.xlabel('Share of drugs approved')
    save('Approval rate by drug type', 'joined_rate_type.png')
    sns.catplot(data=joined, x='outcome', y='pchembl_value', kind='box')
    save('Potency of approved and failed drugs', 'joined_potency.png')

    # Joined: + SEC
    joined['company'] = joined['lead_sponsor'].apply(clean_company)
    sec['company'] = sec['entity_name'].apply(clean_company)
    typical = sec.groupby('company')['rd_millions'].median().reset_index()
    joined = joined.merge(typical, how='left')
    industry = joined[joined['sponsor_class'] == 'Industry']
    summarize('Joined: + SEC (industry sponsors)', industry,
              ['outcome', 'rd_millions'])
    print('Industry drugs matched to SEC:',
          industry['rd_millions'].notna().sum(), 'of', len(industry))
    sns.catplot(data=joined, y='sponsor_class', x='approved', kind='bar',
                errorbar=None)
    plt.xlabel('Share of drugs approved')
    save('Approval rate by sponsor type', 'joined_rate_sponsor.png')
    sns.catplot(data=industry, x='outcome', y='rd_millions', kind='box')
    plt.yscale('log')
    save('Sponsor R&D of approved and failed drugs', 'joined_rd.png')

    # Joined: + FDA orphan list
    orphan_names = orphan['generic_name'].str.upper()
    joined['orphan'] = joined['pref_name'].str.upper().isin(orphan_names)
    summarize('Joined: + FDA orphan list', joined, ['orphan'])
    print(joined.groupby('orphan')['approved'].mean().round(2))
    sns.catplot(data=joined, x='orphan', y='approved', kind='bar',
                errorbar=None)
    plt.ylabel('Share of drugs approved')
    save('Approval rate of orphan drugs', 'joined_rate_orphan.png')
    sns.catplot(data=joined, y='molecule_type', hue='orphan', kind='count')
    save('Orphan drugs by drug type', 'joined_orphan_type.png')

    # Statistical test assumption checks
    print()
    print('========== Test assumption checks ==========')
    for column in ['molecule_type', 'sponsor_class', 'orphan']:
        check_chi_square(joined, column)
    for column in ['pchembl_value', 'alogp', 'rd_millions']:
        check_mann_whitney(joined, column)


if __name__ == '__main__':
    main()