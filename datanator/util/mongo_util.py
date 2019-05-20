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

    def list_all_collections(self):
        '''List all non-system collections within database
        '''
        client = pymongo.MongoClient(
            self.MongoDB, replicaSet=self.replicaSet)  # 400ms max timeout
        db = client[self.db]
        expression = {"name": {"$regex": r"^(?!system\.)"}}
        return db.list_collection_names()

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

    def extract_values(self, obj, key):
        """Pull all values of specified key from nested JSON.
        """
        arr = []

        def extract(obj, arr, key):
            """Recursively search for values of key in JSON tree."""
            if isinstance(obj, dict):
                for k, v in obj.items():
                    if isinstance(v, (dict, list)):
                        extract(v, arr, key)
                    elif k == key:
                        arr.append(v)
            elif isinstance(obj, list):
                for item in obj:
                    extract(item, arr, key)
            return arr

        results = extract(obj, arr, key)

        return results

    def print_schema(self, collection_str):
        '''Print out schema of a collection
        '''
        _, _, collection = self.con_db(collection_str)
        doc = collection.find_one( {} )
        self.print_dict(doc)

    def print_dict(self, dictionary, ident = '', braces=1):
        """ Recursively prints nested dictionaries."""

        for key, value in dictionary.items():
            if isinstance(value, dict):
                print ('%s%s%s%s' %(ident,braces*'[',key,braces*']')) 
                self.print_dict(value, ident+'  ', braces+1)
            else:
                print (ident+'%s = %s' %(key, type(value).__name__))