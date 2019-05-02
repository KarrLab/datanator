'''
Cron jobs for aggregating documents with same keys
'''

from pymongo import MongoClient
import os

class MetaboliteUniprot():
	def __init__(self):
		self.client = MongoClient('mongodb://localhost:27017/')
		self.db = self.client['compounds']
		self.metabolite_collection = self.db['metabolites']
		self.uniprot_collection = self.db['uniprot']

	def aggregate(self):
		self.pipeline = [
			{'$lookup':
				{'from': 'uniprot',
				'localField' : 'enzymes.enzyme.uniprot_id',
				'foreignField' : 'uniprot_id',
				'as' : 'enzyme_detail'}
			}
		]
		self.db.command('aggregate', 'metabolites', pipeline=self.pipeline, explain=True)

if __name__ == '__main__':
	a = MetaboliteUniprot()
	a.aggregate()