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
import json
import pandas
import requests
from datanator.util import mongo_util
import datanator.config.core


class UniprotNoSQL(mongo_util.MongoUtil):
    def __init__(self, MongoDB=None, db=None, max_entries=2000000, verbose=False,
         username=None, password=None, authSource='admin', replicaSet=None, collection_str='uniprot'):
        self.url = 'http://www.uniprot.org/uniprot/?fil=reviewed:yes'
        self.MongoDB = MongoDB
        self.db = db
        self.max_entries = max_entries
        self.collection_str = collection_str
        super(UniprotNoSQL, self).__init__(MongoDB=MongoDB, db=db, username=username,
                                 password=password, authSource=authSource, replicaSet=replicaSet,
                                 verbose=verbose)

      # build dataframe for uniprot_swiss for loading into mongodb
    def get_uniprot(self):
        url = self.url + \
            '&columns=id,entry name,genes(PREFERRED),protein names,sequence,length,mass,ec,database(GeneID),reviewed,organism-id,database(KO)'
        url += '&format=tab'
        url += '&compress=no'
        if not math.isnan(self.max_entries):
           url += '&limit={}'.format(self.max_entries)

        response = requests.get(url)
        response.raise_for_status()

        data = pandas.read_csv(io.BytesIO(response.content),
                               delimiter='\t', encoding='utf-8', nrows = 700000)
        data.columns = [
            'uniprot_id', 'entry_name', 'gene_name', 'protein_name', 'canonical_sequence', 'length', 'mass',
            'ec_number', 'entrez_id', 'status', 'ncbi_taxonomy_id', 'ko_number'
        ]
        data['entrez_id'] = data['entrez_id'].str.replace(';', '')
        data['mass'] = data['mass'].str.replace(',', '')
        data['ko_number'] = data['ko_number'].str.replace(';', '')
        return data

    # load uniprot into MongoDB
    def load_uniprot(self):
        df = self.get_uniprot()
        df_json = json.loads(df.to_json(orient='records'))
        client, db_obj, collection = self.con_db(self.collection_str)
        collection.delete_many({})
        collection.insert(df_json)

        return collection

def main():
    db = 'datanator'
    collection_str = 'uniprot'
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    server = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']
    port = datanator.config.core.get_config(
    )['datanator']['mongodb']['port'] 
    manager=UniprotNoSQL(MongoDB = server, db = db, 
        username = username, password = password, collection_str=collection_str)
    manager.load_uniprot()

if __name__ == '__main__':
    main()