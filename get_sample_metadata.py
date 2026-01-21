#! /usr/bin/python3

import requests 
import json
from datetime import date
import pandas as pd
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("sras_file")
parser.add_argument("output_file")

args = parser.parse_args()

sras_file = args.sras_file
out_file = args.output_file

muninn_url = 'http://kenny.scripps.edu:8000'
# muninn_url = 'https://h5n1.outbreak.info/api'
max_collection_span = pd.Timedelta(30, 'd')


# get sras
sras: list[str]
with open(sras_file, 'r') as f:
    sras = [l.strip() for l in f.readlines()]

# get sample data
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
stride = 250
metadata: pd.DataFrame = pd.DataFrame()
for i in range(0, len(sras), stride):
    chunk = sras[i:i+stride]
    query = '|'.join([f'accession={s}' for s in chunk])

    resp = requests.get(f'{muninn_url}/samples?q={query}')
    print(resp)

    metadata_chunk = pd.DataFrame.from_dict(json.loads(resp.text))[keep_keys]
    metadata = pd.concat([metadata, metadata_chunk], ignore_index=True)

metadata['collection_start_date'] = pd.to_datetime(metadata['collection_start_date'])
metadata['collection_end_date'] = pd.to_datetime(metadata['collection_end_date'])

metadatas_present_sras = set(metadata['accession'])
missing_metadata_sras = set(sras) - metadatas_present_sras

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



