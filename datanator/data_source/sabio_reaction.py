from datanator.util import mongo_util
import datanator.config.core
from pymongo.collation import Collation, CollationStrength
import os

class RxnAggregate:

    def __init__(self, username=None, password=None, server=None, authSource='admin',
                 src_database='datanator', max_entries=float('inf'), verbose=True,
                 collection='protein', destination_database='datanator', cache_dir=None):
        '''
                Args:
                        src_database (:obj: `str`): name of database in which source collections reside
                        destination_database (:obj: `str`): name of database to put the aggregated collection
        '''
        self.client, self.db, self.col = mongo_util.MongoUtil(MongoDB=server, username=username,
                                                              password=password, authSource=authSource,
                                                              db=destination_database).con_db(collection)
        self.collation = Collation(locale='en', strength=CollationStrength.SECONDARY)