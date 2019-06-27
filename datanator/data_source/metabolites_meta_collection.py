from datanator.core import query_nosql
from datanator.util import chem_util
from datanator.util import file_util
from datanator.util import index_collection
import datanator.config.core
import pymongo
import re

class MetabolitesMeta(query_nosql.QuerySabio):
    ''' meta_loc: database location to save the meta collection
    '''

    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet=None, db=None,
                 verbose=False, max_entries=float('inf'), username = None, 
                 password = None, authSource = 'admin', meta_loc = None):
        self.cache_dirname = cache_dirname
        self.verbose = verbose
        self.MongoDB = MongoDB
        self.replicaSet = replicaSet
        self.max_entries = max_entries
        self.username = username
        self.password = password
        self.authSource = authSource
        self.meta_loc = meta_loc
        super(MetabolitesMeta, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, replicaSet=replicaSet,
                                              db=db, verbose=verbose, max_entries=max_entries, username = username,
                                              password = password, authSource = authSource)
        self.frequency = 50
        self.chem_manager = chem_util.ChemUtil()
        self.file_manager = file_util.FileUtil()

    def load_content(self):
        collection_name = 'metabolites_meta'

        # ecmdb_fields = ['m2m_id', 'inchi', 'synonyms.synonym']
        # self.fill_metabolite_fields(
        #     fields=ecmdb_fields, collection_src='ecmdb', collection_des = collection_name)

        # ymdb_fields = ['ymdb_id', 'inchi', 'synonyms.synonym']
        # self.fill_metabolite_fields(
        #     fields=ymdb_fields, collection_src='ymdb', collection_des = collection_name)

        _, _, collection = self.con_db(collection_name)
        # k = 0
        # for doc in self.doc_feeder(collection_str=collection_name, query={}, projection={'inchi'}):
        #     if k > self.max_entries:
        #         break
        #     kinlaw_id = self.get_kinlawid_by_inchi([doc['inchi']])
        #     rxn_participants = self.find_reaction_participants(kinlaw_id)
        #     collection.update_one({'inchi': doc['inchi']},
        #                           {'$set': {'kinlaw_id': kinlaw_id,
        #                            'reaction_participants': rxn_participants}},
        #                           upsert=False)
        #     k += 1
        i = 0
        cursor = collection.find(filter = {}, projection = {'similar_compounds_bak':1, 'similar_compounds': 1})
        for doc in cursor:
            if i % self.frequency == 0:
                print(i)

            replacement = []
            for k, v in doc['similar_compounds_bak'].items():
                dic = {}
                dic[k] = v
                replacement.append(dic)

            collection.update_one({'_id': doc['_id']},
                                 {'$set': {'similar_compounds': replacement}},
                                upsert=False)
            i += 1
        

    def fill_metabolite_fields(self, fields=None, collection_src=None, collection_des = None):
        '''Fill in values of fields of interest from 
            metabolite collection: ecmdb or ymdb
                Args:
                        fileds: list of fields of interest
                        collection_src: collection in which query will be done
                        collection_des: collection in which result will be updated

        '''
        projection = {}
        for field in fields:
            projection[field] = 1
        projection['_id'] = 0
        _, _, col_src = self.con_db(collection_src)
        _, _, col_des = self.con_db(collection_des)
        cursor = col_src.find(filter={}, projection=projection)
        i = 0
        for doc in cursor:
            if i == self.max_entries:
                break
            if i % self.frequency == 0:
                print('Getting fields of interest from {} document in {}'.format(i, collection_src))
            doc['inchi_deprot'] = self.chem_manager.simplify_inchi(inchi = doc['inchi'])
            doc['inchi_hashed'] = self.chem_manager.hash_inchi(inchi = doc['inchi'])
            doc['inchi_hashed_deprot'] = self.chem_manager.hash_inchi(inchi = doc['inchi_deprot'])
            col_des.update_one({'inchi': doc['inchi']},
                                  { '$set': doc},
                                  upsert=True)
            i += 1



def main():
    db = 'datanator'
    meta_loc = 'datanator'
    username = datanator.config.core.get_config()['datanator']['mongodb']['user']
    password = datanator.config.core.get_config()['datanator']['mongodb']['password']
    MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
    port = datanator.config.core.get_config()['datanator']['mongodb']['port']
    replSet = datanator.config.core.get_config()['datanator']['mongodb']['replSet']
    manager = MetabolitesMeta(cache_dirname=None, MongoDB=MongoDB, replicaSet = replSet, db=db, 
                                verbose=True, max_entries=float('inf'), 
                                username = username, password = password, meta_loc = meta_loc)

    manager.load_content()


if __name__ == '__main__':
    main()
