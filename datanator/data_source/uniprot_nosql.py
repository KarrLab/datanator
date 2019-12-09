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
import pymongo.errors
from datanator.util import mongo_util
import datanator.config.core


class UniprotNoSQL(mongo_util.MongoUtil):
    def __init__(self, MongoDB=None, db=None, max_entries=float('inf'), verbose=False,
         username=None, password=None, authSource='admin', replicaSet=None, collection_str='uniprot'):
        self.url = 'http://www.uniprot.org/uniprot/?fil=reviewed:yes'
        self.query_url = 'https://www.uniprot.org/uniprot/?query='
        self.MongoDB = MongoDB
        self.db = db
        self.max_entries = max_entries
        self.collection_str = collection_str
        super(UniprotNoSQL, self).__init__(MongoDB=MongoDB, db=db, username=username,
                                 password=password, authSource=authSource, replicaSet=replicaSet,
                                 verbose=verbose, max_entries=max_entries)
        self.client, self.db, self.collection = self.con_db(collection_str)

    # build dataframe for uniprot_swiss for loading into mongodb
    def load_uniprot(self, query=False, msg='', species=None):
        """Build dataframe
        
        Args:
            query (:obj:`bool`, optional): Whether download all reviewed entries of perform individual queries. Defaults to False.
            msg (:obj:`str`, optional): Query message. Defaults to ''.
            species (:obj:`list`, optional): species information to extract from df and loaded into uniprot. Defaults to None.
        """
        fields = '&columns=id,entry name,genes(PREFERRED),protein names,sequence,length,mass,ec,database(GeneID),reviewed,organism-id,database(KO),genes(ALTERNATIVE),genes(ORF),genes(OLN)'
        if not query:
            url = self.url + fields
        else:
            query_msg = msg
            if isinstance(species, list):
                for specie in species:
                    query_msg += '+'+str(specie)
            url = self.query_url + query_msg + '&sort=score' + fields
        url += '&format=tab'
        url += '&compress=no'
        if not math.isnan(self.max_entries):
           url += '&limit={}'.format(self.max_entries)
        
        response = requests.get(url, stream=False)
        response.raise_for_status()

        try:
            data = pandas.read_csv(io.BytesIO(response.content), delimiter='\t', encoding='utf-8')
        except pandas.errors.EmptyDataError:
            return
        data.columns = [
            'uniprot_id', 'entry_name', 'gene_name', 'protein_name', 'canonical_sequence', 'length', 'mass',
            'ec_number', 'entrez_id', 'status', 'ncbi_taxonomy_id', 'ko_number', 'gene_name_alt',
            'gene_name_orf', 'gene_name_oln'
        ]
        data['entrez_id'] = data['entrez_id'].astype(str).str.replace(';', '')

        data['mass'] = data['mass'].str.replace(',', '')

        data['ko_number'] = data['ko_number'].astype(str).str.replace(';', '')
        data['gene_name_oln'] = data['gene_name_oln'].astype(str).str.split(' ')
        data['gene_name_orf'] = data['gene_name_orf'].astype(str).str.split(' ')
        data['gene_name_alt'] = data['gene_name_alt'].astype(str).str.split(' ')
        if species is None:
            self.load_df(data)
        else:
            self.load_df(data.loc[data['ncbi_taxonomy_id'].isin(species)])

    # load pandas.DataFrame into MongoDB
    def load_df(self, df):
        df_json = json.loads(df.to_json(orient='records'))
        try:        
            self.collection.insert(df_json)
        except pymongo.errors.InvalidOperation as e:
            return(str(e))
            

def main():
    db = 'datanator'
    collection_str = 'uniprot'
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    server = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']
    manager=UniprotNoSQL(MongoDB = server, db = db, 
        username = username, password = password, collection_str=collection_str)
    manager.load_uniprot()

if __name__ == '__main__':
    main()