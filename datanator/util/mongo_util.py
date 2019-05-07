import pymongo
import wc_utils.quilt
from bson import decode_all

class MongoUtil():

    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet=None, db=None,
                verbose=False, max_entries=float('inf')):
        self.cache_dirname = cache_dirname
        self.MongoDB = MongoDB
        self.db = db
        self.replicaSet = replicaSet
        self.verbose = verbose
        self.max_entries = max_entries

    def con_db(self, collection_str):
        try:
            client = pymongo.MongoClient(
                self.MongoDB, replicaSet=self.replicaSet)  # 400ms max timeout
            db = client[self.db]
            collection = db[collection_str]
            return (client, db, collection)
        except pymongo.errors.ConnectionFailure:
            return ('Server not available')

    def fill_db(self, collection_str, sym_link=False):
        '''Check if collection is already in MongoDB 
            if already in:
                do nothing
            else:
                load data into db from quiltdata (karrlab/datanator_nosql)

            Attributes:
                collection_str: name of collection (e.g. 'ecmdb', 'pax', etc)
                sym_link: whether download should be a sym link
        '''
        _, _, collection = self.con_db(collection_str)
        if collection.find({}).count() != 0:
            return collection
        else:
            manager = wc_utils.quilt.QuiltManager(self.cache_dirname, 'datanator_nosql')
            file = collection_str+'.bson' 
            manager.download(file, sym_link)
            with open((self.cache_dirname+ '/'+file), 'rb') as f:
                collection.insert(decode_all(f.read()))
            return collection