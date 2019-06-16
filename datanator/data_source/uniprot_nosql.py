'''
	Generates uniprot_swiss (reviewed) NoSQL documents
	load documents into MongoDB collections

:Author: Zhouyang Lian <zhouyang.lian@familian.life>
:Author: Jonathan <jonrkarr@gmail.com>
:Date: 2019-04-02
:Copyright: 2019, Karr Lab
:License: MIT
'''

import io
import math
import os
import json
import pandas
import requests
import requests.exceptions
from datanator.util import mongo_util

class UniprotNoSQL(mongo_util.MongoUtil):
    def __init__(self, MongoDB = None, db = None, max_entries = float('inf'), verbose = False,
        username = None, password = None, authSource = 'admin', replicaSet = None):
        self.url = 'http://www.uniprot.org/uniprot/?fil=reviewed:yes'
        self.MongoDB = MongoDB
        self.db = db
        self.max_entries = max_entries
        super(UniprotNoSQL, self).__init__(MongoDB = MongoDB, db = db, username = username,
                                password = password, authSource = authSource, replicaSet = replicaSet,
                                verbose = verbose)

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
                               delimiter='\t', encoding='utf-8', nrows = self.max_entries)
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
        client, db_obj, collection = self.con_db('uniprot')
        collection.delete_many({})
        collection.insert(df_json)

        return collection

