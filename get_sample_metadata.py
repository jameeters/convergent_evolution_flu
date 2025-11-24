#! /usr/bin/python3

import requests 
import json
from datetime import date
import pandas as pd

msa_file = 'out/aligned.fa'
muninn_url = 'http://kenny.scripps.edu:8000'
max_collection_span = pd.Timedelta(30, 'd')
out_file = 'out/sample_metadata.tsv'

# get sras
sras = set()
with open(msa_file, 'r') as f:
    for l in f.readlines():
        if l.startswith('>'):
            sra = l.strip().strip('>').strip('_HA')
            sras.add(sra)

# get sample data

query = '|'.join([f'accession={s}' for s in sras])

resp = requests.get(f'{muninn_url}/samples?q={query}')

keep_keys = [
    'accession',
    'collection_start_date',
    'collection_end_date',
    'host',
    'geo_country_name',
    'geo_admin1_name',
    'geo_admin2_name',
    'geo_admin3_name',
]

metadata = pd.DataFrame.from_dict(json.loads(resp.text))[keep_keys]

metadata['collection_start_date'] = pd.to_datetime(metadata['collection_start_date'])
metadata['collection_end_date'] = pd.to_datetime(metadata['collection_end_date'])

metadatas_present_sras = set(metadata['accession'])
missing_metadata_sras = sras - metadatas_present_sras

if len(missing_metadata_sras) > 0:
    print(f'warning: {len(missing_metadata_sras)} SRAs are missing from muninn: {missing_metadata_sras}')


metadata['collection_span'] = metadata['collection_end_date'] - metadata['collection_start_date']


n_dropped_date = metadata[metadata['collection_span'] > pd.Timedelta(days = 0)].shape[0]
if n_dropped_date > 0:
    print(f'warning: {n_dropped_date} rows dropped for collection span > {max_collection_span}')

# filter out wide collection spans
metadata = metadata[metadata['collection_span'] <= max_collection_span]

# calculate "decimal date"
metadata['mid_collection_date'] = metadata['collection_start_date'] + metadata['collection_span'] / 2
metadata['stardate'] = metadata['mid_collection_date'].dt.year + (metadata['mid_collection_date'].dt.day_of_year / 365.25)

with open(out_file, 'w+') as f:
    metadata.to_csv(f, sep='\t', index=False)



