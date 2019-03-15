'''
	Generates metabolite (ECMDB, YMDB) documents for MongoDB
'''


import io
import os
import json
import requests
import requests.exceptions
import warnings
import zipfile
import xmltodict


class MetaboliteNoSQL():
    def __init__(self, verbose, dababase):
        self.verbose = verbose
        self.database = database

        if self.database == 'ecmdb':
            self.domain = 'http://ecmdb.ca'
            self.compound_index = self.domain + \
                '/download/ecmdb.json.zip'  # list of metabolites
            self.compound_url = self.domain + '/compounds/{}.xml'
            self.collection_dir = './cache/ecmdb/'
        else:
            self.domain = 'http://ymdb.ca'
            self.compound_index = self.domain + \
                '/system/downloads/current/ymdb.json.zip'  # list of metabolites
            self.compound_url = self.domain + '/compounds/{}.xml'
            self.collection_dir = './cache/ymdb/'

    # each compound's collection of imformation is written into a json file for mongodb
    def write_to_json(self):

        if self.verbose:
            print('Download list of all compounds: ...')

        response = requests.get(self.compound_index)
        response.raise_for_status()

        if self.verbose:
            print('... Done!')
        if self.verbose:
            print('Unzipping and parsing compound list ...')

        with zipfile.ZipFile(io.BytesIO(response.content), 'r') as zip_file:
            json_name = self.database + '.json'
            with zip_file.open(json_name, 'r') as json_file:
                entries = json.load(json_file)

        if self.verbose:
            print('  found {} compounds'.format(len(entries)))

        if self.database == 'ecmdb':
            entries.sort(key=lambda e: e['m2m_id'])
        else:
            entries.sort(key=lambda e: e['ymdb_id'])

        if self.verbose:
            print('Downloading {} compounds ...'.format(len(entries)))

        for i_entry, entry in enumerate(entries):

            if self.verbose and (i_entry % 10 == 0):
                print('  Downloading compound {} of {}'.format(
                    i_entry + 1, len(entries)))

            # actual xml file for each metabolite
            if self.database == 'ecmdb':
                response = requests.get(
                    self.compound_url.format(entry['m2m_id']))
                try:
                    response.raise_for_status()
                except requests.exceptions.HTTPError:
                    warnings.warn(
                        'Unable to download data for compound {}'.format(entry['m2m_id']))
                    continue
            else:
                response = requests.get(
                    self.compound_url.format(entry['ymdb_id']))
                try:
                    response.raise_for_status()
                except requests.exceptions.HTTPError:
                    warnings.warn(
                        'Unable to download data for compound {}'.format(entry['ymdb_id']))
                    continue

            doc = xmltodict.parse(response.text)

            # delete key "compound" but keep key's value
            new_doc = doc['compound']

            os.makedirs(os.path.dirname(self.collection_dir), exist_ok=True)
            if self.database == 'ecmdb':
                file_name = os.path.join(
                    self.collection_dir + entry['m2m_id'] + '.json')
            else:
                file_name = os.path.join(
                    self.collection_dir + entry['ymdb_id'] + '.json')
            with open(file_name, "w") as f:
                f.write(json.dumps(new_doc))


if __name__ == '__main__':
    ecmdb = MetaboliteNoSQL(True, 'ecmdb')
    ecmdb.write_to_json()

    ymdb = MetaboliteNoSQL(True, 'ymdb')
    ymdb.write_to_json()
