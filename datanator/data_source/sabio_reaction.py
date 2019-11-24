from datanator.util import mongo_util, file_util
import datanator.config.core
from pymongo.collation import Collation, CollationStrength
import os


class RxnAggregate:

    def __init__(self, username=None, password=None, server=None, authSource='admin',
                 src_database='datanator', max_entries=float('inf'), verbose=True,
                 collection='sabio_rk_old', destination_database='datanator', cache_dir=None):
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

    def fill_collection(self):
        projection = {'_id': 0,'resource': 1, 'reaction_participant.substrate_aggregate': 1,
                    'reaction_participant.product_aggregate': 1, 'kinlaw_id': 1}
        _, _, collection = self.mongo_manager.con_db(collection)
        docs = collection.find({})

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