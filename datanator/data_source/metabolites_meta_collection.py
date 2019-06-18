from datanator.core import query_nosql
from datanator.util import chem_util
from datanator.util import server_util
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
        self.frequency = 100
        self.chem_manager = chem_util.ChemUtil()

    def load_content(self):
        ecmdb_fields = ['m2m_id', 'inchi', 'synonyms.synonym']
        ecmdb_list = self.get_metabolite_fields(
            fields=ecmdb_fields, collection_str='ecmdb')

        ymdb_fields = ['ymdb_id', 'inchi', 'synonyms.synonym']
        ymdb_list = self.get_metabolite_fields(
            fields=ymdb_fields, collection_str='ymdb')

        collection_name = 'metabolites_meta'
        client = pymongo.MongoClient(
            self.MongoDB, replicaSet=self.replicaSet, 
            username = self.username, password = self.password,
            authSource = self.authSource)
        meta_db = client[self.meta_loc]
        collection = meta_db[collection_name]

        i = 0
        for doc in ecmdb_list:
            if i > self.max_entries:
                break
            doc['inchi_deprot'] = self.chem_manager.simplify_inchi(inchi = doc['inchi'])
            doc['inchi_hashed'] = self.chem_manager.hash_inchi(inchi = doc['inchi_deprot'])
            collection.update_one({'inchi': doc['inchi']},
                                  { '$set': doc},
                                  upsert=True)
            i += 1

        j = 0
        for doc in ymdb_list:
            if j > self.max_entries:
                break
            doc['inchi_deprot'] = self.chem_manager.simplify_inchi(inchi = doc['inchi'])
            doc['inchi_hashed'] = self.chem_manager.hash_inchi(inchi = doc['inchi_deprot'])
            collection.update_one({'inchi': doc['inchi']},
                                  { '$set': doc},
                                  upsert=True)
            j += 1

        k = 0
        for doc in self.doc_feeder(collection_str=collection_name, query={}, projection={'inchi'}):
            if k > self.max_entries:
                break
            kinlaw_id = self.get_kinlawid_by_inchi([doc['inchi']])
            rxn_participants = self.find_reaction_participants(kinlaw_id)
            collection.update_one({'inchi': doc['inchi']},
                                  {'$set': {'kinlaw_id': kinlaw_id,
                                   'reaction_participants': rxn_participants}},
                                  upsert=False)
            k += 1
        

        client.close()


    def get_metabolite_fields(self, fields=None, collection_str=None):
        '''Get values of fields of interest from 
            metabolite collection: ecmdb or ymdb
                Args:
                        fileds: list of fields of interest
                        collection_str: collection in which query will be done
                Return:
                        list of dictionaries
                                [{'filed1': value1, 'field2': value2, ... }, ..., {}]
        '''
        projection = {}
        for field in fields:
            projection[field] = 1
        projection['_id'] = 0
        cursor = self.doc_feeder(collection_str=collection_str, query={},
                                 projection=projection)
        i = 0
        dict_list = []
        for doc in cursor:
            if i == self.max_entries:
                break
            if i % self.frequency == 0:
                print('Getting fields of interest from {}'.format(collection_str))
            doc['inchi'] = self.chem_manager.simplify_inchi(inchi=doc['inchi'])
            dict_list.append(doc)
            i += 1

        return dict_list


def main():
    db = 'datanator'
    meta_loc = 'datanator'
    config_file = '/root/host/karr_lab/datanator/.config/config.ini'
    username, password, server, port = server_util.ServerUtil(
        config_file=config_file).get_user_config()
    manager = MetabolitesMeta(cache_dirname=None, MongoDB=server, replicaSet = None, db=db, 
                                verbose=True, max_entries=float('inf'), 
                                username = username, password = password, meta_loc = meta_loc)

    manager.load_content()


if __name__ == '__main__':
    main()
