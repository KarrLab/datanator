import pymongo
import wc_utils.quilt
import json
from bson import decode_all
from genson import SchemaBuilder

class MongoUtil():

    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet=None, db=None,
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
        client = pymongo.MongoClient(
            self.MongoDB, replicaSet=self.replicaSet)  # 400ms max timeout
        db = client[self.db]
        expression = {"name": {"$regex": r"^(?!system\.)"}}
        return db.list_collection_names()

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
           removed '_id' from collection due to its object type
           and universality 
        '''
        _, _, collection = self.con_db(collection_str)
        doc = collection.find_one({})
        builder = SchemaBuilder()
        del doc['_id']
        builder.add_object(doc)
        return builder.to_schema()

    def simplify_inchi(self, inchi=None):
        '''Remove molecules's protonation state
        "InChI=1S/H2O/h1H2" = > "InChI=1S/H2O"
        '''
        # if self.verbose:
        #     print('Parsing inchi by taking out protonation state')
        try:
            inchi_neutral = inchi.split('/h')[0]
            return inchi_neutral
        except AttributeError:
            return None

    def flatten_json(self, nested_json):
        '''
            Flatten json object with nested keys into a single level.
            e.g. 
            {a: b,                      {a: b,  
             c: [                        d: e,
                {d: e},    =>            f: g }
                {f: g}]}
            Args:
                nested_json: A nested json object.
            Returns:
                The flattened json object if successful, None otherwise.
        '''
        out = {}

        def flatten(x, name=''):
            if type(x) is dict:
                for a in x:
                    flatten(x[a], name + a + '_')
            elif type(x) is list:
                i = 0
                for a in x:
                    flatten(a, name + str(i) + '_')
                    i += 1
            else:
                out[name[:-1]] = x

        flatten(nested_json)
        return out

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
