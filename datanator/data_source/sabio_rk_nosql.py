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
from pymongo import MongoClient
from pathlib import Path
import re
import os
import wc_utils.quilt


class SabioRkNoSQL():

    def __init__(self, db, MongoDB, cache_directory=None, quilt_package=None, verbose=False, max_entries=float('inf')):
        '''
                Attributes:
                        cache_directory: JSON file (converted from sqlite) directory
                        db: mongodb database name
                        MongoDB: MongoDB server address and login e.g. 'mongodb://localhost:27017/'
                        quilt_package: quilt package name (e.g. datanator)
                        system_path: subdir of quilt package

        '''
        self.db = db
        self.MongoDB = MongoDB
        self.cache_directory = cache_directory
        self.quilt_package = quilt_package
        self.verbose = verbose
        self.max_entries = max_entries
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
            return ('Server not available')

    # load json files
    def load_json(self):
        file_names = []
        file_dict = {}
        if self.cache_directory == None:
            directory = '../../datanator/data_source/cache/SabioRk'
            pathlist = Path(directory).glob('**/*.json')
            for path in pathlist:
                path_in_str = str(path)
                name = re.findall(r'[^\/]+(?=\.json$)', path_in_str)[0]
                file_names.append(name)
                with open(path_in_str) as f:
                    file_dict[name] = json.load(f)
        else:
            manager = wc_utils.quilt.QuiltManager(
                self.cache_directory, self.quilt_package)
            downloads = manager.download()
            directory = self.cache_directory
            for (dirpath, dirnames, filenames) in os.walk(directory):
                [file_names.append(name[:-5])
                 for name in filenames if name.endswith('json')]
                break
            for name in file_names:
                with open(os.path.join(directory, name + '.json')) as f:
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

        for i in range(min(len(kinetic_law_list), self.max_entries)):

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
                    'id': None,
                    'name': None,
                    'type': None,
                    'synonym': None,
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
                'id': cur_enzyme_entry_dict['id'],
                'name': cur_enzyme_entry_dict['name'],
                'type': cur_kinlaw_dict['enzyme_type'],
                'synonym': enzyme_synonym,
                'created': cur_enzyme_entry_dict['created'],
                'modified': cur_enzyme_entry_dict['modified']})

            if len(cur_enzyme_subunit_list) == 0:
                sabio_doc['enzymes'][2]['subunit'].append({
                    'id': None,
                    'name': None,
                    'uniprot_id': None,
                    'coefficient': None,
                    'sequence': None,
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
                        'id': cur_entry_dict['id'],
                        'name': cur_entry_dict['name'],
                        'uniprot_id': uniprot_id,
                        'coefficient': cur_enzyme_subunit_dict['coefficient'],
                        'sequence': cur_enzyme_subunit_dict['sequence'],
                        'molecular_weight': cur_enzyme_subunit_dict['molecular_weight'],
                        'created': cur_entry_dict['created'],
                        'modified': cur_entry_dict['modified']
                    })

            if len(cur_compartment_list) == 0:
                sabio_doc['enzymes'][1]['compartment'].append({
                    'id': None,
                    'name': None,
                    'created': None,
                    'modified': None
                })
            else:
                for j in range(len(cur_compartment_list)):
                    cur_compartment_dict = cur_compartment_list[j]
                    sabio_doc['enzymes'][1]['compartment'].append({
                        'id': cur_compartment_dict['id'],
                        'name': cur_compartment_dict['name'],
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
                                                       'inchi': compound_structure_list[structure_id - 1]['_value_inchi'],
                                                       'inchi_connectivity': compound_structure_list[structure_id - 1]['_value_inchi_formula_connectivity']})

                cur_reactant_compartment_id = cur_reactant_dict.get(
                    'compartment_id')
                if cur_reactant_compartment_id != None:
                    cur_entry_dict = next(
                        item for item in entry_list if item['_id'] == cur_reactant_compartment_id)
                    compartment_name = cur_entry_dict['name']
                    compartment_created = cur_entry_dict['created']
                    compartment_modified = cur_entry_dict['modified']
                    reactant_compartment = {'name': compartment_name,
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
                    'name': cur_entry_dict['name'],
                    'synonym': reactant_synonyms,
                    'structure': cur_compound_structure,
                    'compartment': reactant_compartment,
                    'coefficient': cur_reactant_dict['coefficient'],
                    'type': cur_reactant_dict['type'],
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
                                                       'inchi': compound_structure_list[structure_id - 1]['_value_inchi'],
                                                       'inchi_connectivity': compound_structure_list[structure_id - 1]['_value_inchi_formula_connectivity']})

                cur_reactant_compartment_id = cur_reactant_dict.get(
                    'compartment_id')
                if cur_reactant_compartment_id != None:
                    cur_entry_dict = next(
                        item for item in entry_list if item['_id'] == cur_reactant_compartment_id)
                    compartment_name = cur_entry_dict['name']
                    compartment_created = cur_entry_dict['created']
                    compartment_modified = cur_entry_dict['modified']
                    reactant_compartment = {'name': compartment_name,
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
                    'name': cur_entry_dict['name'],
                    'synonym': reactant_synonyms,
                    'structure': cur_compound_structure,
                    'compartment': reactant_compartment,
                    'coefficient': cur_reactant_dict['coefficient'],
                    'type': cur_reactant_dict['type'],
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
                    reactant_compartment = {'name': compartment_name,
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
                                                       'inchi': compound_structure_list[structure_id - 1]['_value_inchi'],
                                                       'inchi_connectivity': compound_structure_list[structure_id - 1]['_value_inchi_formula_connectivity']})

                cur_reactant_coefficient = cur_reactant_dict['coefficient']
                cur_reactant_type = cur_reactant_dict['type']
                cur_entry_dict = next(
                    item for item in entry_list if item['_id'] == cur_reactant_compound_entry_id)
                sabio_doc['reaction_participant'][2]['modifier'].append({**{
                    'sabio_compound_id': cur_entry_dict['id'],
                    'name': cur_entry_dict['name'],
                    'synonym': reactant_synonyms,
                    'structure': cur_compound_structure,
                    'compartment': reactant_compartment,
                    'coefficient': cur_reactant_dict['coefficient'],
                    'type': cur_reactant_dict['type'],
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
                                'inchi': cur_compound_structure_dict['_value_inchi'],
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
            sabio_doc['taxon'] = cur_kinlaw_dict['taxon']
            sabio_doc['taxon_wildtype'] = cur_kinlaw_dict['taxon_wildtype']
            sabio_doc['taxon_variant'] = cur_kinlaw_dict['taxon_variant']
            sabio_doc['temperature'] = cur_kinlaw_dict['temperature']
            sabio_doc['ph'] = cur_kinlaw_dict['ph']
            sabio_doc['media'] = cur_kinlaw_dict['media']

            with open(json_name, 'w') as f:
                f.write(json.dumps(sabio_doc, indent=4))

            collection = self.con_db()
            collection.insert(sabio_doc)
