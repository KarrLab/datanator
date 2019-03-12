from pymongo import MongoClient
import os
import json

client = MongoClient('mongodb://localhost:27017/')
metabolite_db = client['metabolites']
ecmdb_collection = metabolite_db['ecmdb']

file_dir = './cache/' + 'ecmdb'
for filename in os.listdir(file_dir):
	with open(os.path.join(file_dir, filename)) as f:
		file_data = json.load(f)
		ecmdb_collection.insert(file_data)

client.close()