'''
	merges ECMDB and YMDB NoSQL documents

:Author: Zhouyang Lian <zhouyang.lian@familian.life>
:Author: Jonathan <jonrkarr@gmail.com>
:Date: 2019-04-02
:Copyright: 2019, Karr Lab
:License: MIT
'''

from pymongo import MongoClient
import os
import json
from bson import ObjectId

# handle ObjectID encoded by pymongo
class JSONEncoder(json.JSONEncoder):
	def default(self, o):
		if isinstance(o, ObjectId):
			return str(o)
		return json.JSONEncoder.default(self, o)

if __name__ == '__main__':

	client = MongoClient('mongodb://localhost:27017/')
	compound_db = client['compounds']
	organisms = ['ecmdb', 'ymdb']
	metabolite_collection = compound_db['metabolites']

	for organism in organisms:
		file_dir = './cache/' + organism
		for filename in os.listdir(file_dir):
			with open(os.path.join(file_dir, filename)) as f:
				file_data = json.load(f)
				metabolite_collection.insert(file_data)

	# aggregate ecmdb and ymdb based on keys
	metabolite_collection.aggregate([
	  {"$group":{"_id":'null', "keys":{"$mergeObjects":"$$ROOT"}}},
	  {"$project":{"keys": { "$map": { "input": { "$objectToArray": "$keys" }, "in": "$$this.k" }}}}
	])

	client.close()

