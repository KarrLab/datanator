'''Parse SabioRk json files into MongoDB documents
    (json files acquired by running sqlite_to_json.py)
:Author: Zhouyang Lian <zhouyang.lian@familian.life>
:Author: Jonathan <jonrkarr@gmail.com>
:Date: 2019-04-02
:Copyright: 2019, Karr Lab
:License: MIT
'''

import json
import pymongo
import datanator.config.core
from pymongo import MongoClient
from datanator.util import mongo_util
from datanator.util import chem_util
from pathlib import Path
import re
import os
import datetime


class SabioRkNoSQL(mongo_util.MongoUtil):

    def __init__(self, db = None, MongoDB = None, cache_directory=None, 
                verbose=False, max_entries=float('inf'), replicaSet = None,
                username = None, password = None, authSource = 'admin'):
        '''
                Attributes:
                        cache_directory: JSON file (converted from sqlite) directory
                        db: mongodb database name
                        MongoDB: MongoDB server address and login e.g. 'mongodb://mongo:27017/'
        '''
        self.db = db
        self.MongoDB = MongoDB
        self.cache_directory = cache_directory
        self.verbose = verbose
        self.max_entries = max_entries
        self.collection_str = 'sabio_rk_old'
        super(SabioRkNoSQL, self).__init__(cache_dirname=cache_directory, MongoDB=MongoDB, replicaSet=replicaSet, 
                                    db=db, verbose=verbose, max_entries=max_entries, username = username, 
                                    password = password, authSource = authSource)

        self.client, self.db_obj, self.collection = self.con_db(self.collection_str)
        self.chem_manager = chem_util.ChemUtil()

    # load json files
    def load_json(self):
        file_names = []
        file_dict = {}
        directory = '../../datanator/data_source/cache/SabioRk'
        pathlist = Path(directory).glob('**/*.json')
        for path in pathlist:
            path_in_str = str(path)
            name = re.findall(r'[^\/]+(?=\.json$)', path_in_str)[0]
            file_names.append(name)
            with open(path_in_str) as f:
                file_dict[name] = json.load(f)
        return (file_names, file_dict)

    def make_doc(self, file_names, file_dict):

        os.makedirs(os.path.dirname(
          self.cache_directory + 'docs/'), exist_ok=True)
        null = None
        for key, value in file_dict.items():
            # might not be a good idea
            globals()[key + '_list'] = value

        # list of entrie ids that have synonyms
        has_synonym = list(item['entry__id'] for item in entry_synonym_list)
        # list of entrie ids that have compound structure information
        has_structure = list(item['compound__id']
                             for item in compound_compound_structure_list)
        
        for i in range(34170, min(len(kinetic_law_list), self.max_entries)):

            cur_kinlaw_dict = kinetic_law_list[i]
            kinlaw_id = next(
                item['id'] for item in entry_list if item['_id'] == cur_kinlaw_dict['_id'])
            json_name = self.cache_directory+'docs/' + \
                'kinlaw_id_' + str(kinlaw_id) + '.json'

            if self.verbose and (i % 100 == 0):
                print('  Downloading kinlaw_id {} of {}'.format(
                    kinlaw_id, min(len(kinetic_law_list), self.max_entries)))

            sabio_doc = {}
            sabio_doc['kinlaw_id'] = kinlaw_id
            sabio_doc['resource'] = []
            sabio_doc['enzymes'] = [{'enzyme': []}, {
                'compartment': []}, {'subunit': []}]
            sabio_doc['parameter'] = []
            sabio_doc['reaction_participant'] = [
                {'substrate': []}, {'product': []}, {'modifier': []}]

            # every kinlaw entry has no more than one enzyme
            enzyme_id = cur_kinlaw_dict['enzyme_id']
            if enzyme_id != null:
                cur_enzyme_dict = next(
                    item for item in enzyme_list if item['_id'] == enzyme_id)
                cur_enzyme_entry_dict = next(
                    item for item in entry_list if item['_id'] == enzyme_id)
            else:
                sabio_doc['enzymes'][0]['enzyme'].append({
                    'molecular_weight': None,
                    'enzyme_id': None,
                    'enzyme_name': None,
                    'enzyme_type': None,
                    'enzyme_synonym': None,
                    'created': None,
                    'modified': None})

            cur_enzyme_subunit_list = list(
                item for item in enzyme_subunit_list if item['enzyme_id'] == enzyme_id)  # enzyme_id == entry_id
            enzyme_compartment_id = cur_kinlaw_dict['enzyme_compartment_id']
            # enzyme_compartment_id = compartment_id = entry_id
            cur_compartment_list = list(
                item for item in entry_list if item['_id'] == enzyme_compartment_id)

            enzyme_synonym = []
            if enzyme_id in has_synonym:
                synonym_id_list = list(
                    item['synonym__id'] for item in entry_synonym_list if item['entry__id'] == enzyme_id)
                for synonym_id in synonym_id_list:
                    enzyme_synonym.append(synonym_list[synonym_id - 1]['name'])

            sabio_doc['enzymes'][0]['enzyme'].append({
                'molecular_weight': cur_enzyme_dict['molecular_weight'],
                'enzyme_id': cur_enzyme_entry_dict['id'],
                'enzyme_name': cur_enzyme_entry_dict['name'],
                'enzyme_type': cur_kinlaw_dict['enzyme_type'],
                'enzyme_synonym': enzyme_synonym,
                'created': cur_enzyme_entry_dict['created'],
                'modified': cur_enzyme_entry_dict['modified']})

            if len(cur_enzyme_subunit_list) == 0:
                sabio_doc['enzymes'][2]['subunit'].append({
                    'subunit_id': None,
                    'subunit_name': None,
                    'uniprot_id': None,
                    'subunit_coefficient': None,
                    'canonical_sequence': None,
                    'molecular_weight': None,
                    'created': None,
                    'modified': None
                })
            else:
                for j in range(len(cur_enzyme_subunit_list)):
                    cur_enzyme_subunit_dict = cur_enzyme_subunit_list[j]
                    entry_id = cur_enzyme_subunit_dict['_id']
                    resource_id = next(
                        item['resource__id'] for item in entry_resource_list if item['entry__id'] == entry_id)
                    uniprot_id = next(
                        item['id'] for item in resource_list if item['_id'] == resource_id)
                    cur_entry_dict = next(
                        item for item in entry_list if item['_id'] == entry_id)

                    sabio_doc['enzymes'][2]['subunit'].append({
                        'subunit_id': cur_entry_dict['id'],
                        'subunit_name': cur_entry_dict['name'],
                        'uniprot_id': uniprot_id,
                        'subunit_coefficient': cur_enzyme_subunit_dict['coefficient'],
                        'canonical_sequence': cur_enzyme_subunit_dict['sequence'],
                        'molecular_weight': cur_enzyme_subunit_dict['molecular_weight'],
                        'created': cur_entry_dict['created'],
                        'modified': cur_entry_dict['modified']
                    })

            if len(cur_compartment_list) == 0:
                sabio_doc['enzymes'][1]['compartment'].append({
                    'compartment_id': None,
                    'compartment_name': None,
                    'created': None,
                    'modified': None
                })
            else:
                for j in range(len(cur_compartment_list)):
                    cur_compartment_dict = cur_compartment_list[j]
                    sabio_doc['enzymes'][1]['compartment'].append({
                        'compartment_id': cur_compartment_dict['id'],
                        'compartment_name': cur_compartment_dict['name'],
                        'created': cur_compartment_dict['created'],
                        'modified': cur_compartment_dict['modified']
                    })

            # kinetic_law_resource many-to-one
            cur_kinlaw_resource_list = list(
                item for item in kinetic_law_resource_list if item['kinetic_law__id'] == cur_kinlaw_dict['_id'])
            # entry_resource: same entry id can have multiple resource_id
            cur_entry_resource_list = list(
                item for item in entry_resource_list if item['entry__id'] == cur_kinlaw_dict['_id'])

            if len(cur_kinlaw_resource_list) != 0:
                resource_id = cur_kinlaw_resource_list[0]['resource__id']
                cur_resource_dict = next(
                    item for item in resource_list if item['_id'] == resource_id)
                sabio_doc['resource'].append({
                    'namespace': cur_resource_dict['namespace'],
                    'id': cur_resource_dict['id']
                })

            if len(cur_entry_resource_list) != 0:
                for j in range(len(cur_entry_resource_list)):
                    resource__id = cur_entry_resource_list[j]['resource__id']
                    cur_resource_dict = next(
                        item for item in resource_list if item['_id'] == resource__id)
                    sabio_doc['resource'].append({
                        'namespace': cur_resource_dict['namespace'],
                        'id': cur_resource_dict['id']
                    })

            '''Handling reaction_participant document
                substrate, product, modifier
            '''
            cur_reaction_reactant_list = list(
                item for item in reaction_participant_list if item['reactant_kinetic_law_id'] == cur_kinlaw_dict['_id'])
            cur_reaction_product_list = list(
                item for item in reaction_participant_list if item['product_kinetic_law_id'] == cur_kinlaw_dict['_id'])
            cur_reaction_modifier_list = list(
                item for item in reaction_participant_list if item['modifier_kinetic_law_id'] == cur_kinlaw_dict['_id'])

            # substrates
            for j in range(len(cur_reaction_reactant_list)):
                cur_reactant_dict = cur_reaction_reactant_list[j]
                cur_reactant_compound_entry_id = cur_reactant_dict['compound_id']

                resource__id = list(
                    item['resource__id'] for item in entry_resource_list if item['entry__id'] == cur_reactant_compound_entry_id)
                resources = {
                    resource_list[x-1]['namespace']: resource_list[x-1]['id'] for x in resource__id}

                # standardize across collections
                resources['kegg_id'] = resources.pop('kegg.compound', None)
                resources['pubchem_substance_id'] = resources.pop(
                    'pubchem.substance', None)
                resources['pubchem_compound_id'] = resources.pop(
                    'pubchem.compound', None)
                # compound synonym if there is any
                cur_reactant_synonym_id_list = list(
                    item['synonym__id'] for item in entry_synonym_list if item['entry__id'] == cur_reactant_compound_entry_id)
                if len(cur_reactant_synonym_id_list) != 0:
                    reactant_synonyms = [synonym_list[x-1]['name']
                                         for x in cur_reactant_synonym_id_list]
                else:
                    reactant_synonyms = []

                # compound structure if there is any
                cur_compound_structure = []
                if cur_reactant_compound_entry_id in has_structure:
                    structure_id_list = list(
                        item['compound_structure__id'] for item in compound_compound_structure_list if item['compound__id'] == cur_reactant_compound_entry_id)
                    for structure_id in structure_id_list:
                        cur_compound_structure.append({'value': compound_structure_list[structure_id - 1]['value'],
                                                       'format': compound_structure_list[structure_id - 1]['format'],
                                                       'inchi_structure': compound_structure_list[structure_id - 1]['_value_inchi'],
                                                       'inchi_connectivity': compound_structure_list[structure_id - 1]['_value_inchi_formula_connectivity']})

                cur_reactant_compartment_id = cur_reactant_dict.get(
                    'compartment_id')
                if cur_reactant_compartment_id != None:
                    cur_entry_dict = next(
                        item for item in entry_list if item['_id'] == cur_reactant_compartment_id)
                    compartment_name = cur_entry_dict['name']
                    compartment_created = cur_entry_dict['created']
                    compartment_modified = cur_entry_dict['modified']
                    reactant_compartment = {'compartment_name': compartment_name,
                                            'created': compartment_created,
                                            'modified': compartment_modified}
                else:
                    reactant_compartment = {}

                cur_reactant_coefficient = cur_reactant_dict['coefficient']
                cur_reactant_type = cur_reactant_dict['type']
                cur_entry_dict = next(
                    item for item in entry_list if item['_id'] == cur_reactant_compound_entry_id)
                sabio_doc['reaction_participant'][0]['substrate'].append({**{
                    'sabio_compound_id': cur_entry_dict['id'],
                    'substrate_name': cur_entry_dict['name'],
                    'substrate_synonym': reactant_synonyms,
                    'substrate_structure': cur_compound_structure,
                    'substrate_compartment': reactant_compartment,
                    'substrate_coefficient': cur_reactant_dict['coefficient'],
                    'substrate_type': cur_reactant_dict['type'],
                    'created': cur_entry_dict['created'],
                    'modified': cur_entry_dict['modified']
                }, **resources})

            # product
            for j in range(len(cur_reaction_product_list)):
                cur_reactant_dict = cur_reaction_product_list[j]
                cur_reactant_compound_entry_id = cur_reactant_dict['compound_id']

                resource__id = list(
                    item['resource__id'] for item in entry_resource_list if item['entry__id'] == cur_reactant_compound_entry_id)
                resources = {
                    resource_list[x-1]['namespace']: resource_list[x-1]['id'] for x in resource__id}

                # standardize across collections
                resources['kegg_id'] = resources.pop('kegg.compound', None)
                resources['pubchem_compound_id'] = resources.pop(
                    'pubchem.substance', None)
                resources['pubchem_compound_id'] = resources.pop(
                    'pubchem.compound', None)
                # compound synonym if there is any
                cur_reactant_synonym_id_list = list(
                    item['synonym__id'] for item in entry_synonym_list if item['entry__id'] == cur_reactant_compound_entry_id)
                if len(cur_reactant_synonym_id_list) != 0:
                    reactant_synonyms = [synonym_list[x-1]['name']
                                         for x in cur_reactant_synonym_id_list]
                else:
                    reactant_synonyms = []

                # compound structure if there is any
                cur_compound_structure = []
                if cur_reactant_compound_entry_id in has_structure:
                    structure_id_list = list(
                        item['compound_structure__id'] for item in compound_compound_structure_list if item['compound__id'] == cur_reactant_compound_entry_id)
                    for structure_id in structure_id_list:
                        cur_compound_structure.append({'value': compound_structure_list[structure_id - 1]['value'],
                                                       'format': compound_structure_list[structure_id - 1]['format'],
                                                       'inchi_structure': compound_structure_list[structure_id - 1]['_value_inchi'],
                                                       'inchi_connectivity': compound_structure_list[structure_id - 1]['_value_inchi_formula_connectivity']})

                cur_reactant_compartment_id = cur_reactant_dict.get(
                    'compartment_id')
                if cur_reactant_compartment_id != None:
                    cur_entry_dict = next(
                        item for item in entry_list if item['_id'] == cur_reactant_compartment_id)
                    compartment_name = cur_entry_dict['name']
                    compartment_created = cur_entry_dict['created']
                    compartment_modified = cur_entry_dict['modified']
                    reactant_compartment = {'compartment_name': compartment_name,
                                            'created': compartment_created,
                                            'modified': compartment_modified}
                else:
                    reactant_compartment = {}

                cur_reactant_coefficient = cur_reactant_dict['coefficient']
                cur_reactant_type = cur_reactant_dict['type']
                cur_entry_dict = next(
                    item for item in entry_list if item['_id'] == cur_reactant_compound_entry_id)
                sabio_doc['reaction_participant'][1]['product'].append({**{
                    'sabio_compound_id': cur_entry_dict['id'],
                    'product_name': cur_entry_dict['name'],
                    'product_synonym': reactant_synonyms,
                    'product_structure': cur_compound_structure,
                    'product_compartment': reactant_compartment,
                    'product_coefficient': cur_reactant_dict['coefficient'],
                    'product_type': cur_reactant_dict['type'],
                    'created': cur_entry_dict['created'],
                    'modified': cur_entry_dict['modified']
                }, **resources})

            # modifier
            for j in range(len(cur_reaction_modifier_list)):
                cur_modifier_dict = cur_reaction_modifier_list[j]
                cur_reactant_compound_entry_id = cur_modifier_dict['compound_id']

                resource__id = list(
                    item['resource__id'] for item in entry_resource_list if item['entry__id'] == cur_reactant_compound_entry_id)
                resources = {
                    resource_list[x-1]['namespace']: resource_list[x-1]['id'] for x in resource__id}

                # standardize across collections
                resources['kegg_id'] = resources.pop('kegg.compound', None)
                resources['pubchem_compound_id'] = resources.pop(
                    'pubchem.substance', None)
                resources['pubchem_compound_id'] = resources.pop(
                    'pubchem.compound', None)
                # compound synonym if there is any
                cur_reactant_synonym_id_list = list(
                    item['synonym__id'] for item in entry_synonym_list if item['entry__id'] == cur_reactant_compound_entry_id)
                if len(cur_reactant_synonym_id_list) != 0:
                    reactant_synonyms = [synonym_list[x-1]['name']
                                         for x in cur_reactant_synonym_id_list]
                else:
                    reactant_synonyms = []

                cur_reactant_compartment_id = cur_modifier_dict.get(
                    'compartment_id')
                if cur_reactant_compartment_id != None:
                    cur_entry_dict = next(
                        item for item in entry_list if item['_id'] == cur_reactant_compartment_id)
                    compartment_name = cur_entry_dict['name']
                    compartment_created = cur_entry_dict['created']
                    compartment_modified = cur_entry_dict['modified']
                    reactant_compartment = {'compartment_name': compartment_name,
                                            'created': compartment_created,
                                            'modified': compartment_modified}
                else:
                    reactant_compartment = {}

                # compound structure if there is any
                cur_compound_structure = []
                if cur_reactant_compound_entry_id in has_structure:
                    structure_id_list = list(
                        item['compound_structure__id'] for item in compound_compound_structure_list if item['compound__id'] == cur_reactant_compound_entry_id)
                    for structure_id in structure_id_list:
                        cur_compound_structure.append({'value': compound_structure_list[structure_id - 1]['value'],
                                                       'format': compound_structure_list[structure_id - 1]['format'],
                                                       'inchi_structure': compound_structure_list[structure_id - 1]['_value_inchi'],
                                                       'inchi_connectivity': compound_structure_list[structure_id - 1]['_value_inchi_formula_connectivity']})

                cur_reactant_coefficient = cur_reactant_dict['coefficient']
                cur_reactant_type = cur_reactant_dict['type']
                cur_entry_dict = next(
                    item for item in entry_list if item['_id'] == cur_reactant_compound_entry_id)
                sabio_doc['reaction_participant'][2]['modifier'].append({**{
                    'sabio_compound_id': cur_entry_dict['id'],
                    'modifier_name': cur_entry_dict['name'],
                    'modifier_synonym': reactant_synonyms,
                    'modifier_structure': cur_compound_structure,
                    'modifier_compartment': reactant_compartment,
                    'modifier_coefficient': cur_reactant_dict['coefficient'],
                    'modifier_type': cur_reactant_dict['type'],
                    'created': cur_entry_dict['created'],
                    'modified': cur_entry_dict['modified']
                }, **resources})

            '''Handling parameter
            '''
            cur_parameter_list = list(
                item for item in parameter_list if item['kinetic_law_id'] == cur_kinlaw_dict['_id'])
            if len(cur_parameter_list) != 0:
                for j in range(len(cur_parameter_list)):
                    entry__id = cur_parameter_list[j]['_id']
                    cur_entry_dict = next(
                        item for item in entry_list if item['_id'] == entry__id)

                    _type = cur_parameter_list[j]['type']

                    compound_id = cur_parameter_list[j]['compound_id']
                    cur_compound_list = list(
                        item for item in compound_list if item['_id'] == compound_id)
                    if len(cur_compound_list) == 0:
                        cur_compound_dict = {}
                    else:
                        cur_compound_dict = cur_compound_list[0]
                    compound_structure = []
                    cur_compound_compound_structure_list = list(
                        item for item in compound_compound_structure_list if item['compound__id'] == compound_id)
                    if len(cur_compound_compound_structure_list) != 0:
                        for k in range(len(cur_compound_compound_structure_list)):
                            compound_structure__id = cur_compound_compound_structure_list[
                                k]['compound_structure__id']
                            cur_compound_structure_dict = next(
                                item for item in compound_structure_list if item['_id'] == compound_structure__id)
                            compound_structure.append({
                                'structure': cur_compound_structure_dict['value'],
                                'format': cur_compound_structure_dict['format'],
                                'inchi_structure': cur_compound_structure_dict['_value_inchi'],
                                'inchi_connectivity': cur_compound_structure_dict['_value_inchi_formula_connectivity']})

                    synonyms = []
                    cur_entry_synonym_list = list(
                        item for item in entry_synonym_list if item['entry__id'] == entry__id)
                    if len(cur_entry_synonym_list) != 0:
                        for k in range(len(cur_entry_synonym_list)):
                            cur_entry_synonym_dict = cur_entry_synonym_list[k]
                            synonym__id = cur_entry_synonym_dict['synonym__id']
                            synonyms.append(
                                synonym_list[synonym__id-1]['name'])

                    # enzyme_id = cur_parameter_list[j]['enzyme_id']
                    if cur_parameter_list[j]['compartment_id'] != None:
                        compartment_id = cur_parameter_list[j]['compartment_id']
                        para_compartment = next(
                            item for item in entry_list if item['_id'] == compartment_id)['name']
                    else:
                        para_compartment = None

                    value = cur_parameter_list[j]['value']
                    error = cur_parameter_list[j]['error']
                    units = cur_parameter_list[j]['units']
                    observed_name = cur_parameter_list[j]['observed_name']
                    observed_type = cur_parameter_list[j]['observed_type']
                    observed_value = cur_parameter_list[j]['observed_value']
                    observed_error = cur_parameter_list[j]['observed_error']
                    observed_units = cur_parameter_list[j]['observed_units']

                    if compound_id == None:
                        entry_compound_id = None
                        compound_name = None
                        compound_created = None
                        compound_modified = None
                    else:
                        cur_entry_dict = next(
                            item for item in entry_list if item['_id'] == compound_id)
                        entry_compound_id = cur_entry_dict['id']
                        compound_name = cur_entry_dict['name']
                        compound_created = cur_entry_dict['created']
                        compound_modified = cur_entry_dict['modified']

                    sabio_doc['parameter'].append({
                        # 'id': cur_entry_dict['id'],
                        'name': cur_entry_dict['name'],
                        # 'type': _type,
                        'sabio_compound_id': entry_compound_id,
                        'compound_name': compound_name,
                        'compound_structure': compound_structure,
                        'compound_created': compound_created,
                        'compound_modified': compound_modified,
                        'is_name_ambiguous': cur_compound_dict.get('_is_name_ambiguous'),
                        'compound_synonyms': synonyms,
                        'compartment': para_compartment,
                        'value': value,
                        'error': error,
                        'units': units,
                        'observed_name': observed_name,
                        'sbo_type': observed_type,
                        'observed_value': observed_value,
                        'observed_error': observed_error,
                        'observed_units': observed_units
                    })

            sabio_doc['tissue'] = cur_kinlaw_dict['tissue']
            sabio_doc['mechanism'] = cur_kinlaw_dict['mechanism']
            sabio_doc['equation'] = cur_kinlaw_dict['equation']
            sabio_doc['taxon_id'] = cur_kinlaw_dict['taxon']
            sabio_doc['taxon_wildtype'] = cur_kinlaw_dict['taxon_wildtype']
            sabio_doc['taxon_variant'] = cur_kinlaw_dict['taxon_variant']
            sabio_doc['temperature'] = cur_kinlaw_dict['temperature']
            sabio_doc['ph'] = cur_kinlaw_dict['ph']
            sabio_doc['media'] = cur_kinlaw_dict['media']

            with open(json_name, 'w') as f:
                f.write(json.dumps(sabio_doc, indent=4))
            
            self.collection.replace_one(
                {'kinlaw_id': sabio_doc['kinlaw_id']},
                sabio_doc,
                upsert=True
                )

    def add_inchi_hash(self, ids=None):
        '''
            Add inchi key values of _value_inchi in sabio_rk collection
        '''
        
        projection = {'reaction_participant.product': 1, 'reaction_participant.substrate':1, 
                     'kinlaw_id': 1}
        if ids is None:
            query = {}
        else:
            query = {'kinlaw_id': {'$in': ids}}
        cursor = self.collection.find(filter=query, projection=projection)
        count = self.collection.count_documents(query)

        def get_inchi_structure(chem):
            '''Given subsrate or product subdocument from sabio_rk
               find the corresponding inchi
               Args:
                    chem (:obj:`dict`)
                Returns:
                    (:obj:`str`): inchi string
            '''
            try:
                x = chem['substrate_structure'][0].get('inchi_structure', None)
                return x
            except KeyError:
                try:
                    x = chem['product_structure'][0].get('inchi_structure', None)
                    return x
                except IndexError:
                    return None
            except IndexError:
                return None

        def iter_rxnp_subdoc(rxnp, side='substrate'):
            '''Given a substrate or product array from sabio_rk
                append fields of hashed inchi
                Args:
                    rxnp (:obj: `list` of :obj: `dict`)
            '''
            if side == 'substrate':
                key = 'substrate_structure'
            else:
                key = 'product_structure'
            for i, rxn in enumerate(rxnp):
                substrate_inchi = get_inchi_structure(rxn)
                try:
                    hashed_inchi = self.chem_manager.inchi_to_inchikey(substrate_inchi)
                    rxnp[i][key][0]['InChI_Key'] = hashed_inchi
                except AttributeError:
                    rxnp[i][key][0]['InChI_Key'] = None
                except IndexError:
                    rxnp[i][key] = []

            return rxnp

        j = 0
        for doc in cursor:
            if j > self.max_entries:
                break
            if self.verbose == True and j % 100 == 0:
                print('Hashing compounds in kinlaw {} out of {}'.format(j, count))
            substrates = doc['reaction_participant'][0]['substrate']
            products = doc['reaction_participant'][1]['product']
            new_subsrates = iter_rxnp_subdoc(substrates)
            new_products = iter_rxnp_subdoc(products, side='product')

            doc['reaction_participant'][0]['substrate'] = new_subsrates
            doc['reaction_participant'][1]['product'] = new_products

            self.collection.update_one({'kinlaw_id': doc['kinlaw_id']},
                        {'$set': {'reaction_participant.0.substrate': doc['reaction_participant'][0]['substrate'],
                                'reaction_participant.1.product': doc['reaction_participant'][1]['product'],
                                'modified': datetime.datetime.utcnow()} })
            j += 1


def main():
    db = 'datanator'
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    MongoDB = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']
    manager = SabioRkNoSQL(username=username, password=password, verbose=True,
                           db=db, MongoDB=MongoDB, cache_directory='./cache/')
    # file_names, file_dict = manager.load_json()
    # manager.make_doc(file_names, file_dict)
    manager.add_inchi_hash()

if __name__ == '__main__':
    main()