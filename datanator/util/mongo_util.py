import pymongo
import wc_utils.quilt
from bson import decode_all
import hashlib
from genson import SchemaBuilder


class MongoUtil:

    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet=None, db='test',
                 verbose=False, max_entries=float('inf'), username = None, 
                 password = None, authSource = 'admin'):
        self.cache_dirname = cache_dirname
        self.MongoDB = MongoDB
        self.db = db
        self.replicaSet = replicaSet
        self.verbose = verbose
        self.max_entries = max_entries
        self.client = pymongo.MongoClient(
            self.MongoDB, replicaSet=self.replicaSet, 
            username = username, password = password,
            authSource = authSource)  # 400ms max timeout
        self.db_obj = self.client[db]

    def list_all_collections(self):
        '''List all non-system collections within database
        '''

        return self.db_obj.list_collection_names()

    def con_db(self, collection_str):
        try:
            collection = self.db_obj[collection_str]
            return (self.client, self.db_obj, collection)
        except pymongo.errors.ConnectionFailure:
            return ('Server not available')
        except ServerSelectionTimeoutError:
            return ('Server timeout')

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
            manager = wc_utils.quilt.QuiltManager(
                self.cache_dirname, 'datanator_nosql')
            file = collection_str+'.bson'
            manager.download(file, sym_link)
            with open((self.cache_dirname + '/'+file), 'rb') as f:
                collection.insert(decode_all(f.read()))
            return collection

    def print_schema(self, collection_str):
        '''Print out schema of a collection
           removed '_id' from collection due to its object type
           and universality 
        '''
        _, _, collection = self.con_db(collection_str)
        doc = collection.find_one({})
        builder = SchemaBuilder()
        del doc['_id']
        builder.add_object(doc)
        return builder.to_schema()

    def flatten_collection(self, collection_str):
        '''Flatten a collection

            c is ommitted because it does not have a non-object 
            value associated with it
        '''
        _, _, collection = self.con_db(collection_str)

        pipeline = [
            { "$addFields": { "subdoc.a": "$a" } },
            { "$replaceRoot": { "newRoot": "$subdoc" }  }
        ]
        flat_col = collection.aggregate(pipeline)
        return flat_col
