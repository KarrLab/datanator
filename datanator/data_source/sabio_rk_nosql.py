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
            returns a json document in the format:
            {
                'kinlaw_id': 1
                'compartment:' []
                ...
            }
    '''
    def make_doc(self, file_names, file_dict):
    	'''
    	Start iteration with kinetic_law

        '''
        os.makedirs(os.path.dirname(self.directory + 'docs/'), exist_ok=True)
        null = None
        for i in range(len(file_names)):
            '''
            kinetic_law_dict = file_dict['kinetic_law']
            compound_dict = file_dict['compound']
            ...

            '''
        	vars()[file_names[i] + '_list'] = file_dict[file_names[i]]

        for i in range(len(kinetic_law_list)):
        	json_name = 'kinlaw_id_' + str(i)
        	cur_kinlaw_dict = kinetic_law_list[i]
        	sabio_doc = {}
        	sabio_doc['kinlaw_id'] = i
        	for key in cur_kinlaw_dict:
        		if key != '_id':
                    '''
                    ...
                    enzyme_id = cur_kinlaw_dict['enzyme_id']
                    enzyme_compartment_id = cur_kinlaw_dict['enzyme_compartment_id']
                    ...

                    '''
        			vars()[key] = cur_kinlaw_dict[key]
        		else:
        			kinlaw_entry_id = cur_kinlaw_dict[key]
        	
        	sabio_doc['enzyme_compartment'] = [] # done
        	sabio_doc['compound'] = []
        	sabio_doc['compound_structure'] = []
        	sabio_doc['resource'] = [] #done
        	sabio_doc['enzyme'] = [] # done
        	sabio_doc['synonym'] = [] 
        	sabio_doc['enzyme_subunit'] = [] # done
        	sabio_doc['kinetic_law_resource'] = [] # done
        	sabio_doc['parameter'] = []
        	sabio_doc['reaction_participant']

            # every kinlaw entry has one and only one enzyme
            cur_enzyme_dict = next(item for item in enzyme_list if item['_id'] == enzyme_id)
            to_append = {
        	'molecular_weight': cur_enzyme_dict['molecular_weight'] 
        	'id': entry_list[enzyme_id - 1]['id']
        	'name': entry_list[enzyme_id - 1]['name']
        	'created': entry_list[enzyme_id - 1]['created']
        	'modified': entry_list[enzyme_id - 1]['modified']}

            sabio_doc['enzyme'].append(to_append)

            cur_compartment_dict = next(item for item in compartment_list if item['_id'] == enzyme_compartment_id)
            if not cur_compartment_dict:
                to_append = {
                'id': None
                'name': None
                'created': None
                'modified': None
                }
            else:
                to_append = {
                'id': entry_list[enzyme_compartment_id - 1]['id']
                'name': entry_list[enzyme_compartment_id - 1]['name']
                'created': entry_list[enzyme_compartment_id - 1]['created']
                'modified': entry_list[enzyme_compartment_id - 1]['modified']
                }
            sabio_doc['enzyme_compartment'].append(to_append)

            # kinetic_law_resource many-to-one
            cur_kinlaw_resource_dict = next(item for item in kinteic_law_resource_list if item['kinetic_law__id'] == kinlaw_entry_id)
            # entry_resource: same entry id can have multiple resource_id
            cur_entry_resource_list = list(item for item in entry_resource if item['entry__id'] == kinlaw_entry_id)

            to_append = []
            if not cur_kinwlaw_resource_dict:
                to_append.append({
                'namespace': None
                'id': None  
                })
            else:
                resource_id = cur_kinlaw_resource_dict['resource__id']
                cur_resource_dict = next(item for item in resource_list if item['_id'] == resource_id)
                to_append.append({
                'namespace': cur_resource_dict['namespace']
                'id': cur_resource_dict['id']
                })

            if len(cur_entry_resource_list) != 0:
                for j in range(len(cur_entry_resource_list)):
                    resource__id = cur_entry_resource_list[j]['resource__id']
                    cur_resource_dict = next(item for item in resource_list if item['_id'] == resource__id)
                    to_append.append({
                        'namespace': cur_resource_dict['namespace']
                        'id': cur_resource_dict['id']
                        })
            sabio_doc['resource'].append(to_append)

            cur_enzyme_subunit_dict = next(item for item in enzyme_subunit_list if item['enzyme_id'] == enzyme_id)
            if not cur_enzyme_subunit_dict:
                to_append = {
                'id': None
                'name': None
                'coefficient': None
                'sequence': None
                'molecular_weight': None
                'created': None
                'modified': None
                }   
            else:
                entry_id = cur_enzyme_subunit_dict['_id']
                to_append = {
                'id': entry_list[entry_id - 1]['id']
                'name': entry_list[entry_id -1]['name']
                'coefficient': cur_enzyme_subunit_dict['coefficient']
                'sequence': cur_enzyme_subunit_dict['sequence']
                'molecular_weight': cur_enzyme_subunit_dict['molecular_weight']
                'created': entry_list[entry_id -1]['created']
                'modified': entry_list[entry_id -1]['modified']
                }
            sabio_doc['enzyme_subunit'].append(to_append)

            to_append = []
            cur_parameter_list = list(item for item in parameter_list if item['kinetic_law_id'] == kinlaw_entry_id)
            if len(cur_parameter_list) != 0:
                for j in range(len(cur_parameter_list)):
                    entry__id = cur_parameter_list[j]['_id']
                    cur_entry_dict = next(item for item in entry_list if item['_id'] == entry__id)

                    _type = cur_parameter_list[j]['type']

                    compound_id = cur_parameter_list[j]['compound_id']
                    cur_compound_dict = next(item for item in compound_list if item['_id'] == compound_id)
                    compound_structure = []
                    cur_compound_compound_structure_list = list(item for item in compound_compound_structure_list if item['compound__id'] == compound_id)
                    if len(cur_compound_compound_structure_list) != 0:
                        for k in range(len(cur_compound_compound_structure_list)):
                            compound_structure__id = cur_compound_compound_structure_list[k]['compound_structure__id']
                            cur_compound_structure_dict = next(item for item in compound_structure_list if item['_id'] == compound_strcture__id)
                            compound_structure.append({
                                'structure': cur_compound_structure_dict['value'],
                                'format': cur_compound_structure_dict['format'],
                                'inchi': cur_compound_structure_dict['_value_inchi'],
                                'inchi_connectivity': cur_compound_structure_dict['_value_inchi_formula_connectivity']})

                    sabio_doc['compound_structure'].append(compound_structure)
                    synonyms = []
                    cur_entry_synonym_list = list(item for item in entry_list if item['entry__id'] == entry__id)
                    if len(cur_entry_synonym_list) != 0:
                        for k in range(len(cur_entry_synonym_list)):
                            cur_entry_synonym_dict = cur_entry_synonym_list[k]
                            synonym__id = cur_entry_synonym_dict['synonym__id']
                            synonyms.append(synonym_list[synonym__id-1]['name'])

                    #enzyme_id = cur_parameter_list[j]['enzyme_id']
                    compartment_id = cur_parameter_list[j]['compartment_id']
                    para_compartment = entry_list[compartment_id-1]['name']

                    value = cur_parameter_list[j]['value']
                    error = cur_parameter_list[j]['error']
                    units = cur_parameter_list[j]['units']
                    observed_name = cur_parameter_list[j]['observed_name']
                    observed_type = cur_parameter_list[j]['observed_type']
                    observed_value = cur_parameter_list[j]['observed_value']
                    observed_error = cur_parameter_list[j]['observed_error']
                    observed_units = cur_parameter_list[j]['observed_units']

                    to_append.append({
                        'id': cur_entry_dict['id'],
                        'name': cur_entry_dict['name'],
                        'created': cur_entry_dict['created'],
                        'modified': cur_entry_dict['modified'],
                        'type': _type,
                        'compound_id': entry_list[compound_id -1]['id'],
                        'compound_name': entry_list[compound_id-1]['name'],
                        'compound_structure': compound_structure,
                        'compound_created': entry_list[compound_id-1]['created'],
                        'compound_modified': entry_list[compound_id-1]['modified'],
                        'is_name_ambiguous': cur_compound_dict['_is_name_ambiguous'],
                        'compound_synonyms': synonyms,
                        'compartment': para_compartment,
                        'value': value,
                        'error': error,
                        'units': units,
                        'observed_name': observed_name,
                        'observed_type': observed_type,
                        'observed_value': observed_value,
                        'observed_error': observed_error,
                        'observed_units': observed_units
                        })
            sabio_doc['parameter'].append(to_append)


            '''
            single value key/value pairs
            '''
            sabio_doc['enzyme_type'] = cur_kinlaw_dict['enzyme_type']
            sabio_doc['tissue'] = cur_kinlaw_dict['tissue']
            sabio_doc['mechanism'] = cur_kinlaw_dict['mechanism']
            sabio_doc['equation'] = cur_kinlaw_dict['equation']
            sabio_doc['taxon'] = cur_kinlaw_dict['taxon']
            sabio_doc['taxon_wildtype'] = cur_kinlaw_dict['taxon_wildtype']
            sabio_doc['taxon_variant'] = cur_kinlaw_dict['taxon_variant']
            sabio_doc['temperature'] = cur_kinlaw_dict['temperature']
            sabio_doc['ph'] = cur_kinlaw_dict['ph']
            sabio_doc['media'] = cur_kinlaw_dict['media']

            with open(json_name, 'w') as f:
                f.write(json.dumps(sabio_doc))

