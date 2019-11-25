from datanator.util import mongo_util, file_util
import datanator.config.core
from pymongo.collation import Collation, CollationStrength
from pymongo import ASCENDING
import os
import tempfile


class RxnAggregate:

    def __init__(self, username=None, password=None, server=None, authSource='admin',
                 src_database='datanator', max_entries=float('inf'), verbose=True,
                 collection='sabio_reaction', destination_database='datanator', cache_dir=None):
        '''
                Args:
                        src_database (:obj: `str`): name of database in which source collections reside
                        destination_database (:obj: `str`): name of database to put the aggregated collection
        '''
        self.client, self.db, self.col = mongo_util.MongoUtil(MongoDB=server, username=username,
                                                              password=password, authSource=authSource,
                                                              db=destination_database).con_db(collection)
        self.mongo_manager = mongo_util.MongoUtil(MongoDB=server, username=username,
                                                  password=password, authSource=authSource, db=src_database)
        self.file_manager = file_util.FileUtil()
        self.collation = Collation(locale='en', strength=CollationStrength.SECONDARY)
        self.verbose = verbose
        self.max_entries = max_entries

    def fill_collection(self):
        projection = {'_id': 0,'resource': 1, 'reaction_participant.substrate_aggregate': 1,
                    'reaction_participant.product_aggregate': 1, 'kinlaw_id': 1}
        _, _, collection = self.mongo_manager.con_db('sabio_rk_old')
        docs = collection.find({})
        count = collection.count_documents({})
        start = 58496
        for i, doc in enumerate(docs[start:]):
            if self.verbose and i % 100 == 0:
                print('Processing document {} out of {}'.format(i+start, count))
            if doc.get('resource') is None:
                continue
            if i == self.max_entries:
                break
            kinlaw_id = doc['kinlaw_id']
            rxn_id = self.get_rxn_id(doc)
            reactants = self.create_reactants(doc)
            self.col.update_one({'rxn_id': rxn_id},
                                {'$addToSet': {'kinlaw_id': kinlaw_id},
                                '$set': {'substrates': reactants['substrate_aggregate'],
                                        'products': reactants['product_aggregate']}}, upsert=True)
            if i == 0:
                self.col.create_index([("rxn_id", ASCENDING)], background=True)

    def get_rxn_id(self, doc):
        resource = doc['resource']
        sr = self.file_manager.search_dict_list(resource, 'namespace', 'sabiork.reaction')
        _id = sr[0]['id']
        return int(_id)    
    
    def create_reactants(self, doc):
        result = {}
        substrate_aggregate = doc['reaction_participant'][3]['substrate_aggregate']
        product_aggregate = doc['reaction_participant'][4]['product_aggregate']
        result['substrate_aggregate'] = substrate_aggregate
        result['product_aggregate'] = product_aggregate

        return result

def main():
    cache_dirname = tempfile.mkdtemp()
    cache_dir = os.path.join(cache_dirname, 'logs.txt')
    src_db = 'datanator'
    des_db = 'datanator'
    collection_str = 'sabio_reaction_entries'
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    server = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']
    port = datanator.config.core.get_config(
    )['datanator']['mongodb']['port']        
    src = RxnAggregate(username=username, password=password, server=server, 
                        authSource='admin', src_database=src_db,
                        verbose=True, collection=collection_str, destination_database=des_db,
                        cache_dir=cache_dir)
    src.fill_collection()

if __name__ == '__main__':
    main()