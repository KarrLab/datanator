'''Maps collection in MongoDB into Elasticsearch
   (mongo-connector does not support elasticsearch 6.0 or above)

:Author: Zhouyang Lian <zhouyang.lian@familian.life>
:Author: Jonathan <jonrkarr@gmail.com>
:Date: 2019-04-02
:Copyright: 2019, Karr Lab
:License: MIT
'''

import requests
from elasticsearch import Elasticsearch
from pymongo import MongoClient
import pymongo
from bson import ObjectId
from elasticsearch.serializer import JSONSerializer
from elasticsearch import SerializationError

# default Elasticsearch serializer cannot handle BSON objects
class SetEncoder(JSONSerializer):
    def default(self, o):
        if isinstance(o, ObjectId):
            return str(o)
        return JSONSerializer.default(self, o)

class Mapper():
    def __init__(self, elastic, MongoDB, source, cache_dir, verbose=True, max_entries = float('inf')):
        self.elastic = elastic
        self.MongoDB = MongoDB
        self.source = source
        self.cache_dir = cache_dir
        self.verbose = verbose
        self.max_entries = max_entries

    def elastic_con(self):
        return requests.get(self.elastic)

    def con_db(self):
        try:
            client = MongoClient(self.MongoDB, 400)  # 400ms max timeout
            client.server_info()
            db = client['datanator']
            collection = db[self.source]
            return client, collection
        except pymongo.errors.ConnectionFailure:
            return ('Server not available')

    def map_to(self):
        if self.elastic_con().status_code != 200:
            return ('elasticsearch server not reachable')
        else:
            es = Elasticsearch(self.elastic, serializer=SetEncoder())
            client, collection = self.con_db()
            total_compound = collection.count_documents({})
            if total_compound > self.max_entries:
            	collection = collection.find()[0:self.max_entries+1]
            else:
            	collection = collection.find(no_cursor_timeout=True)
            for i, document in enumerate(collection):
            	if self.verbose and (i % 100 == 0):
            		print ('Indexing {} compound of {}...'.format(i, min(self.max_entries, total_compound)))
            	document['id'] = document.pop('_id')
            	try:
            		es.index(index=self.source, doc_type=self.source, id=i, body=document, ignore=400)
            	except SerializationError as f:
            		print ('compound {} cannot be serialized by indexer'.format(i))
            		continue
            client.close()

