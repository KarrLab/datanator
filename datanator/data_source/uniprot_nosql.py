'''
	Generates uniprot_swiss (reviewed) NoSQL documents
	load documents into MongoDB collections
'''

import io
import math
import os
import json
import pandas
import requests
import requests.exceptions
from pymongo import MongoClient


class UniprotNoSQL():
    def __init__(self, MongoDB, db, max_entries = float('inf')):
        self.url = 'http://www.uniprot.org/uniprot/?fil=reviewed:yes'
        self.MongoDB = MongoDB
        self.db = db
        self.max_entries = max_entries

    # make connections wth mongoDB
    def con_db(self):
        try:
            client = MongoClient(self.MongoDB, 400)  # 400ms max timeout
            client.server_info()
            db = client[self.db]
            collection = db['uniprot']
            return collection
        except pymongo.errors.ConnectionFailure:
            print('Server not available')

    # build dataframe for uniprot_swiss for loading into mongodb
    def get_uniprot(self):
        url = self.url + \
            '&columns=id,entry name,genes(PREFERRED),protein names,sequence,length,mass,ec,database(GeneID),reviewed'
        url += '&format=tab'
        url += '&compress=no'
        if not math.isnan(self.max_entries):
            url += '&limit={}'.format(self.max_entries)

        response = requests.get(url)
        response.raise_for_status()

        data = pandas.read_csv(io.BytesIO(response.content),
                               delimiter='\t', encoding='utf-8')
        data.columns = [
            'uniprot_id', 'entry_name', 'gene_name', 'protein_name', 'canonical_sequence', 'length', 'mass',
            'ec_number', 'entrez_id', 'status',
        ]
        data['entrez_id'] = data['entrez_id'].str.replace(';', '')
        data['mass'] = data['mass'].str.replace(',', '')
        return data

    # load uniprot into MongoDB
    def load_uniprot(self):
        df = self.get_uniprot()
        df_json = json.loads(df.to_json(orient='records'))
        collection = self.con_db()
        collection.insert(df_json)

        return collection

