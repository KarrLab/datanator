'''
	Generates uniprot_swiss (reviewed) NoSQL documents
	load documents into MongoDB collections
'''

import io
import math
import os
import json
import pandas
import requests
import requests.exceptions
from pymongo import MongoClient

class UniprotNoSQL():
	def __init__(self):
		self.url = 'http://www.uniprot.org/uniprot/?fil=reviewed:yes'
		self.client = MongoClient('mongodb://localhost:27017/')
		self.db = self.client['compounds']
		self.collection = self.db['uniprot']

	# build dataframe for uniprot_swiss for loading into mongodb
	def get_uniprot(self):
		url = self.url + '&columns=id,entry name,genes(PREFERRED),protein names,sequence,length,mass,ec,database(GeneID),reviewed'
		url += '&format=tab'
		url += '&compress=no'

		response = requests.get(url)
		response.raise_for_status()
		
		data = pandas.read_csv(io.BytesIO(response.content), delimiter='\t', encoding='utf-8')
		data.columns = [
			'uniprot_id', 'entry_name', 'gene_name', 'protein_name', 'canonical_sequence', 'length', 'mass',
			'ec_number', 'entrez_id', 'status',
		]
		data['entrez_id'] = data['entrez_id'].str.replace(';', '')
		data['mass'] = data['mass'].str.replace(',', '')
		return data

	# load uniprot into MongoDB
	def load_uniprot(self, df):
		df_json = json.loads(df.to_json(orient='records'))
		self.collection.insert(df_json)

if __name__ == '__main__':
	a = UniprotNoSQL()
	frame = a.get_uniprot()
	a.load_uniprot(frame)