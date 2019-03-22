'''Parse SabioRk json files into MongoDB documents
	(json files acquired by running sqlite_to_json.py)
'''

import json
import pymongo
from pymongo import MongoClient
from pathlib import Path
import re


class SabioRkNoSQL():

    def __init__(self, directory, db, MongoDB):
        '''
                Attributes:
                        directory: temporary os directory
                        db: mongodb database name
                        MongoDB: MongoDB server address and login e.g. 'mongodb://localhost:27017/'

        '''
        self.directory = directory
        self.MongoDB = MongoDB
        self.db = db
        self.collection = 'sabio_rk'

    # make connections wth mongoDB
    def con_db(self):
        try:
            client = MongoClient(self.MongoDB, 400)  # 400ms max timeout
            client.server_info()
            db = client[self.db]
            collection = db[self.collection]
            return collection
        except pymongo.errors.ConnectionFailure:
            print('Server not available')

    # load json files
    def load_json(self):
        pathlist = Path(self.directory).glob('**/*.json')
        file_names = []
        file_dict = {}
        for path in pathlist:
            path_in_str = str(path)
            name = re.findall(r'[^\/]+(?=\.json$)', path_in_str)[0]
            file_names.append(name)
            with open(path_in_str) as f:
                file_dict[name] = json.load(f)
        return (file_names, file_dict)

    '''loads dictionaries as documents (1 kinetic_law per document)
    	Attributes:
    		
    		file_names: list of name of files
    		file_dict: dictionaries of info for each file e.g.
    					{'entry': {...}, 'kinetic_law': {...},  }
    '''
    def make_doc(self, file_names, file_dict):
    	'''
    	Start iteration with kinetic_law

        '''
        for i in range(len(file_names)):
        	vars()[file_names[i] + '_dict'] = file_dict[file_names[i]]

        for i in range(len(kinetic_law_dict)):
        	json_name = 'kinlaw_id_' + str(i)
        	cur_kinlaw_dict = kinetic_law_dict[i]
        	sabio_doc = {}
        	sabio_doc['kinlaw_id'] = i
        	for key in cur_dict:
        		if key != '_id':
        			vars()[key] = cur_kinlaw_dict[key]
        		else:
        			kinlaw_entry_id = cur_dict[key]
        	
        	sabio_doc['compartment'] = []
        	sabio_doc['compound'] = []
        	sabio_doc['compound_structure'] = []
        	sabio_doc['resource'] = []
        	sabio_doc['enzyme'] = []
        	sabio_doc['synonym'] = []
        	sabio_doc['enzyme_subunit'] = []
        	sabio_doc['kinetic_law'] = []
        	sabio_doc['parameter'] = []
        	sabio_doc['reaction_participant']

        	cur_enzyme_dict = next(item for item in enzyme_dict if enzyme_dict['_id'] == enzyme_id)
        	enzyme_mol_weight = cur_enzyme_dict['molecular_weight'] 
        	enzyme_id = entry_dict[enzyme_id - 1]['id']
        	enzyme_name = entry_dict[enzyme_id - 1]['name']
        	enzyme_created = entry_dict[enzyme_id - 1]['created']
        	enzyme_modified = entry_dict[enzyme_id - 1]['modified']

        	substrates = []
        	products = []


