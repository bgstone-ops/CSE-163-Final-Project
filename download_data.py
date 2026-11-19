"""
Benjamin Stone
CSE 163 AD
Downloads a sample of the project data and saves each dataset as a
CSV file in the data folder.

Before running this file, export the FDA orphan drug list from
https://www.accessdata.fda.gov/scripts/opdlisting/oopd/index.cfm
(search with no filters, export to Excel) and save it as
data/orphan_designations.xlsx.
"""
import os
import time

import pandas as pd
import requests

CHEMBL = 'https://www.ebi.ac.uk/chembl/api/data/'
TRIALS = 'https://clinicaltrials.gov/api/v2/studies?format=csv&pageSize=1000'
SEC = 'https://data.sec.gov/api/xbrl/frames/us-gaap/' \
      'ResearchAndDevelopmentExpense/USD/CY'

# The SEC asks you to say who you are. Use your own name and email.
HEADERS = {'User-Agent': 'Your Name your_email@uw.edu'}
NUM_DRUGS = 1000


def get_json(url: str) -> dict | None:
    """
    Takes in a url and returns the JSON data at that url as a
    dictionary. Returns None if the page does not exist. If the server
    is busy, waits 10 seconds and tries again, up to 5 times in all,
    then raises an Exception.
    """
    for attempt in range(5):
        response = requests.get(url, headers=HEADERS)
        if response.status_code == 404:
            return None
        if response.status_code == 200 and response.text.startswith('{'):
            return response.json()
        print('  Server busy (status', response.status_code,
              '), trying again in 10 seconds...')
        time.sleep(10)
    raise Exception('Could not download ' + url)


def get_chembl(table: str, id_column: str, ids: list[str], key: str,
               extra: str = '') -> list[dict]:
    """
    Takes in a ChEMBL table name, the name of its id column, a list of
    ids, the name of the list of records in ChEMBL's response, and
    optional extra filters (such as '&assay_type=A'). Returns a list
    of the records in that table for the given ids.
    """
    records = []
    for i in range(0, len(ids), 50):
        id_text = ','.join(ids[i:i + 50])
        url = CHEMBL + table + '.json?limit=1000&' + id_column + \
            '__in=' + id_text + extra
        data = get_json(url)
        records = records + data[key]
    return records


def download_drugs() -> list[str]:
    """
    Saves a random sample of NUM_DRUGS ChEMBL drugs that reached at
    least phase 1 to data/drugs.csv, counting each parent molecule
    once, and returns a list of their ChEMBL ids. The sample is the
    same every time it is run.
    """
    rows = []
    offset = 0
    total = 1
    while offset < total:
        url = CHEMBL + 'molecule.json?max_phase__gte=1&limit=1000&offset=' \
            + str(offset)
        data = get_json(url)
        total = data['page_meta']['total_count']
        for m in data['molecules']:
            props = m['molecule_properties']
            if props is None:
                props = {}
            parent = m['molecule_hierarchy']
            if parent is None:
                parent = {'parent_chembl_id': m['molecule_chembl_id']}
            rows.append({
                'molecule_chembl_id': m['molecule_chembl_id'],
                'parent_chembl_id': parent.get('parent_chembl_id'),
                'pref_name': m['pref_name'],
                'molecule_type': m['molecule_type'],
                'max_phase': m['max_phase'],
                'first_approval': m['first_approval'],
                'full_mwt': props.get('full_mwt'),
                'alogp': props.get('alogp'),
                'psa': props.get('psa'),
                'num_ro5_violations': props.get('num_ro5_violations')
            })
        offset = offset + 1000
    df = pd.DataFrame(rows)
    df = df[df['molecule_chembl_id'] == df['parent_chembl_id']]
    df = df.sample(n=min(NUM_DRUGS, len(df)), random_state=163)
    df.to_csv('data/drugs.csv', index=False)
    return list(df['molecule_chembl_id'])


def download_links(drug_ids: list[str]) -> list[str]:
    """
    Takes in a list of ChEMBL drug ids, saves the ClinicalTrials.gov
    trial ids that ChEMBL lists for each drug to data/links.csv, and
    returns a list of the unique trial ids.
    """
    indications = get_chembl('drug_indication', 'molecule_chembl_id',
                             drug_ids, 'drug_indications')
    rows = []
    for indication in indications:
        for ref in indication['indication_refs']:
            if ref['ref_type'] == 'ClinicalTrials':
                for nct_id in ref['ref_id'].split(','):
                    rows.append({
                        'molecule_chembl_id': indication['molecule_chembl_id'],
                        'nct_id': nct_id.strip()
                    })
    links = pd.DataFrame(rows).drop_duplicates()
    links.to_csv('data/links.csv', index=False)
    return list(links['nct_id'].unique())


def download_targets(drug_ids: list[str]) -> list[dict]:
    """
    Takes in a list of ChEMBL drug ids and saves each drug's targets,
    with each target's type and top two levels of protein class, to
    data/targets.csv. Returns a list of dictionaries, one per drug and
    target, with the keys molecule_chembl_id and target_chembl_id.
    """
    mechanisms = get_chembl('mechanism', 'molecule_chembl_id', drug_ids,
                            'mechanisms')
    mech_rows = []
    target_ids = []
    for m in mechanisms:
        if m['target_chembl_id'] is not None:
            mech_rows.append({'molecule_chembl_id': m['molecule_chembl_id'],
                              'target_chembl_id': m['target_chembl_id']})
            target_ids.append(m['target_chembl_id'])

    targets = get_chembl('target', 'target_chembl_id',
                         list(set(target_ids)), 'targets')
    target_rows = []
    component_ids = []
    for t in targets:
        for component in t['target_components']:
            target_rows.append({'target_chembl_id': t['target_chembl_id'],
                                'target_type': t['target_type'],
                                'component_id': component['component_id']})
            component_ids.append(str(component['component_id']))

    components = get_chembl('target_component', 'component_id',
                            list(set(component_ids)), 'target_components')
    component_rows = []
    class_ids = []
    for c in components:
        for pc in c['protein_classifications']:
            class_id = pc['protein_classification_id']
            component_rows.append({'component_id': c['component_id'],
                                   'protein_class_id': class_id})
            class_ids.append(str(class_id))

    classes = get_chembl('protein_classification', 'protein_class_id',
                         list(set(class_ids)), 'protein_classifications')
    class_rows = []
    for pc in classes:
        description = pc.get('protein_class_desc')
        if description is None:
            description = pc.get('pref_name')
        if description is not None:
            levels = description.split('  ')
            class_l1 = levels[0].strip().capitalize()
            class_l2 = None
            if len(levels) > 1:
                class_l2 = levels[1].strip().capitalize()
            class_rows.append({'protein_class_id': pc['protein_class_id'],
                               'class_l1': class_l1,
                               'class_l2': class_l2})

    df = pd.DataFrame(mech_rows).merge(pd.DataFrame(target_rows))
    df = df.merge(pd.DataFrame(component_rows))
    df = df.merge(pd.DataFrame(class_rows))
    df.to_csv('data/targets.csv', index=False)
    return mech_rows


def download_potency(pairs: list[dict]) -> None:
    """
    Takes in a list of dictionaries with the keys molecule_chembl_id
    and target_chembl_id, and saves every potency result
    (pchembl_value) of each drug against that target to
    data/potency.csv.
    """
    rows = []
    for i in range(len(pairs)):
        drug = pairs[i]['molecule_chembl_id']
        target = pairs[i]['target_chembl_id']
        url = CHEMBL + 'activity.json?limit=1000&pchembl_value__isnull=false' \
            + '&molecule_chembl_id=' + drug + '&target_chembl_id=' + target
        for a in get_json(url)['activities']:
            rows.append({'molecule_chembl_id': drug,
                         'target_chembl_id': target,
                         'pchembl_value': a['pchembl_value'],
                         'standard_type': a['standard_type']})
    pd.DataFrame(rows, columns=['molecule_chembl_id', 'target_chembl_id',
                                'pchembl_value', 'standard_type']).to_csv(
        'data/potency.csv', index=False)


def download_trials(nct_ids: list[str]) -> None:
    """
    Takes in a list of ClinicalTrials.gov trial ids and saves each
    trial's start date and year, status, phases, lead sponsor and
    sponsor type to data/trials.csv.
    """
    tables = []
    for i in range(0, len(nct_ids), 100):
        url = TRIALS + '&filter.ids=' + ','.join(nct_ids[i:i + 100])
        tables.append(pd.read_csv(url, storage_options=HEADERS))
    df = pd.concat(tables)
    df = df[['NCT Number', 'Start Date', 'Study Status', 'Phases',
             'Sponsor', 'Funder Type']].copy()
    df.columns = ['nct_id', 'start_date', 'overall_status', 'phases',
                  'lead_sponsor', 'sponsor_class']
    df['start_year'] = pd.to_numeric(df['start_date'].str[:4],
                                     errors='coerce')
    df.to_csv('data/trials.csv', index=False)


def download_sec() -> None:
    """
    Saves each public company's yearly R&D spending for 2000 through
    2022 to data/sec.csv. Years the SEC has no data for are skipped.
    """
    rows = []
    for year in range(2000, 2023):
        url = SEC + str(year) + '.json'
        data = get_json(url)
        if data is None:
            print('  No SEC data for', year)
        else:
            for company in data['data']:
                rows.append({'entity_name': company['entityName'],
                             'year': year,
                             'rd_expense': company['val']})
    pd.DataFrame(rows).to_csv('data/sec.csv', index=False)


def convert_orphan() -> None:
    """
    Saves the generic name, designation date and status columns of the
    FDA orphan drug list in data/orphan_designations.xlsx to
    data/orphan.csv.
    """
    df = pd.read_excel('data/orphan_designations.xlsx')
    for column in df.columns:
        if 'generic' in column.lower():
            name_column = column
        elif 'date' in column.lower() and 'designat' in column.lower():
            date_column = column
        elif 'status' in column.lower():
            status_column = column
    result = pd.DataFrame({'generic_name': df[name_column],
                           'designation_date': df[date_column],
                           'status': df[status_column]})
    result.to_csv('data/orphan.csv', index=False)


def main():
    os.makedirs('data', exist_ok=True)
    print('Downloading ChEMBL drugs...')
    drug_ids = download_drugs()
    print('Downloading trial links, targets and potency results...')
    nct_ids = download_links(drug_ids)
    pairs = download_targets(drug_ids)
    download_potency(pairs)
    print('Downloading ClinicalTrials.gov trials...')
    download_trials(nct_ids)
    print('Downloading SEC R&D data...')
    download_sec()
    convert_orphan()
    print('Done!')


if __name__ == '__main__':
    main()