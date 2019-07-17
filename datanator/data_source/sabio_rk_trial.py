import datanator.config.core
from datanator.util import mongo_util


class SabioRk:

    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet=None, db=None,
                 verbose=False, max_entries=float('inf'), username=None,
                 password=None, authSource='admin'):
        self.cache_dirname = cache_dirname
        self.MongoDB = MongoDB
        self.replicaSet = replicaSet
        self.db = db
        self.verbose = verbose
        self.max_entries = max_entries
        self.client, self.db_obj, self.collection = mongo_util.MongoUtil(
            MongoDB=MongoDB, db=db, username=username, password=password, authSource=authSource)
