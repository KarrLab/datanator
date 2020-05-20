import pandas as pd
import datetime
from datanator_query_python.query import query_metabolites_meta
from collections import deque
import datanator.config.core
from pymongo.collation import Collation, CollationStrength


class Concentration(query_metabolites_meta.QueryMetabolitesMeta):

    def __init__(self, MongoDB=None, db=None, collection_str=None, username=None,
                 password=None, authSource='admin', readPreference='nearest',
                 verbose=True, max_entries=float('inf')):
        """Init        
        """
        super().__init__(MongoDB=MongoDB, db=db, username=username, password=password,
                         authSource=authSource, readPreference=readPreference)
        self.collation = Collation('en', strength=CollationStrength.SECONDARY)
        self.collection = self.db_obj[collection_str]
        self.taxon_collection = self.client['datanator']['taxon_tree']
        self.meta_collection = self.client['datanator']['metabolites_meta']
        self.e_collection = self.client['datanator']['ecmdb']
        self.y_collection = self.client['datanator']['ymdb']
        self.max_entries = max_entries
        self.verbose = verbose
        self.collection_str = collection_str

    def merge_duplicates(self):
        _, docs = self.get_duplicates(self.collection_str, 'inchikey')
        for doc in docs:
            obj_id = doc['uniqueIds'][0]
            for unique_id in doc['uniqueIds'][1:]:
                concentrations = self.collection.find_one(filter={'_id': unique_id})['concentrations']
                self.collection.update_one({'_id': obj_id},
                                           {'$addToSet': {'concentrations': {'$each': concentrations}}})
                self.collection.delete_one({'_id': unique_id})
        return "Done"

    def edit_nchembio2077(self, start=0):
        """Resolve https://github.com/KarrLab/datanator_rest_api/issues/73
        """
        con_0 = {"concentrations.reference.doi": "10.1038/nchembio.2077"}
        con_1 = {"metabolite": {"$ne": "UDP-D-glucose"}}
        query = {"$and": [con_0, con_1]}
        projection = {"concentrations": 1}
        docs = self.collection.find(filter=query, projection=projection, skip=start)
        count = self.collection.count_documents(query) - start
        for i, doc in enumerate(docs):
            if i == self.max_entries:
                break
            if self.verbose and i % 20 == 0:
                print("Processing doc {} out of {} ...".format(i, count))
            seen_val = deque() # previous concentration value
            next_idx = 0 # next new concentration value's index in conc_obj
            conc_obj = []
            conc_idx = {}  # position lookup table
            for conc in doc['concentrations']:
                val = conc['concentration']
                ref = conc["reference"]
                if ref is None or conc.get('ec_number') is None:
                    conc_obj.append(conc)
                    next_idx += 1
                else:
                    ref = ref.get('doi')
                    if ref == "10.1038/nchembio.2077" and val not in seen_val:  # new concentration
                        conc_obj.append({
                                        "affinities": [{
                                                        "ec_number": conc['ec_number'],
                                                        "k_m": conc.get('k_m'),
                                                        "k_i": conc.get('k_i'),
                                                        "unit": conc.get('unit')
                                                        }],
                                        "concentration": conc['concentration'],
                                        "species_name": conc['species_name'],
                                        "ncbi_taxonomy_id": conc["ncbi_taxonomy_id"],
                                        "last_modified": datetime.datetime.utcnow(),
                                        "reference": {"doi": ref}
                                        })                    
                        conc_idx[val] = next_idx
                        seen_val.append(val)
                        next_idx += 1
                    elif ref =="10.1038/nchembio.2077" and val in seen_val: # existing concentration value
                        idx = conc_idx[val] # where in conc_obj the concentration value is
                        conc_obj[idx]["affinities"].append({
                                                            "ec_number": conc['ec_number'],
                                                            "k_m": conc.get('k_m'),
                                                            "k_i": conc.get('k_i'),
                                                            "unit": conc.get('unit')
                                                            })
                        conc_obj[idx]["last_modified"] = datetime.datetime.utcnow()
                    else:
                        conc_obj.append(conc)
                        next_idx += 1
            self.collection.update_one({'_id': doc['_id']},
                                       {'$set': {"concentrations": conc_obj}})

    def grab_eymdb(self, start=0):
        """Fill collection with concentration values from ec/ymdb.
        Documents in current collection all have inchikey
        """
        query = {'concentrations': {'$ne': None}}
        e_docs = self.e_collection.find(filter=query)
        count_e = self.e_collection.count_documents(query)
        y_docs = self.y_collection.find(filter=query)
        count_y = self.y_collection.count_documents(query)

        for i, e_doc in enumerate(e_docs[start:]):
            if self.verbose and i % 50 == 0:
                print('Processing ecmdb doc {} out of {}'.format(i+start, count_e))
            if e_doc['synonyms']:
                synonyms = e_doc['synonyms']['synonym']
            else:
                synonyms = None
            concentrations = self._flatten_conc_obj(e_doc['concentrations'], 562, 'Escherichia coli')
            self.collection.update_many({'inchikey': e_doc['inchikey']},
                                        {'$set': {'metabolite': e_doc['name'],
                                                    'chebi_id': e_doc['chebi_id'],
                                                    'synonyms': synonyms},
                                        '$addToSet': {'concentrations': {'$each': concentrations}}},
                                        upsert=True)
        for i, y_doc in enumerate(y_docs[start:]):
            if self.verbose and i % 50 == 0:
                print('Processing ymdb doc {} out of {}'.format(i+start, count_y))
            if y_doc['synonyms']:
                synonyms = y_doc['synonyms']['synonym']
            else:
                synonyms = None
            concentrations = self._flatten_conc_obj(y_doc['concentrations'], 4932, 'Saccharomyces cerevisiae')
            self.collection.update_many({'inchikey': y_doc['inchikey']},
                                        {'$set': {'metabolite': y_doc['name'],
                                                    'chebi_id': y_doc['chebi_id'],
                                                    'synonyms': synonyms},
                                        '$addToSet': {'concentrations': {'$each': concentrations}}},
                                        upsert=True)

    def remove_eymdb(self):
        """Remove concentration records pulled from EC/YMDB
        """
        self.collection.update_many({},
                                    {'$pull': {'concentrations': {'growth_meda': {'$exists': True}}}})

    def _flatten_conc_obj(self, obj, ncbi_id, species_name):
        """Flatten concentrations object in ec/ymdb.

        Args:
            obj(:obj:`Obj`): concentrations object.
            ncbi_id(:obj:`int`): NCBI Taxonomy ID of species.
            species_name(:obj:`str`): Name of species.

        Return:
            (:obj:`list` of :obj:`Obj`): [{'growth_media': '', 'concentration': '', ...} ..., {}]
        """
        result = []
        typecheck = list(obj.values())[0]
        _type = isinstance(typecheck, list)
        if not _type: #string
            obj['ncbi_taxonomy_id'] = ncbi_id
            obj['species_name'] = species_name
            return [obj]
        else:
            keys = list(obj.keys())
            vals = list(obj.values())
            count = len(vals[0]) # number of final objects in array
            for i in range(0, count):
                dic = {'ncbi_taxonomy_id': ncbi_id, 'species_name': species_name}
                for key in keys:
                    dic[key] = obj[key][i]
                result.append(dic)
            return result
            


    def fill_collection(self, file_location, sheet_name='Sheet1', start_row=0,
                        use_columns='A:E,G', column_names=['EC_Number', 'metabolite',
                        'concentration', 'k_m', 'abbreviation', 'species_name'],
                        reference={'doi': '10.1038/nchembio.2077'}, unit='M'):
        """
        Fill collection collection_str.

        Args:
            sheet_name(:obj:`str`, optional): Name of sheet in excel.
            start_row (:obj:`int`, optional): Read from csv row. Defaults to 0.
            use_columns(:obj:`str`): Indicates comma separated list of Excel column letters and column ranges (e.g. "A:E" or "A,C,E:F"). Ranges are inclusive of both sides.
            column_names(:obj:`list` of :obj:`str`): Names of columns used.
            reference(:obj:`Obj`): reference information.
            unit(:obj:`str`): Unit used for Km.
        """
        df = pd.read_excel(file_location, sheet_name=sheet_name,
                         header=0, usecols=use_columns,
                         skiprows=start_row)
        df.columns = [x.lower() for x in column_names]
        row_count = len(df.index)
        for i, row in df.iterrows():
            if i == self.max_entries:
                break
            if i % 10 == 0 and self.verbose:
                print("Processing locus {} out {}".format(i, row_count))
            metabolite = row['metabolite']
            con_0 = {'name': metabolite}
            con_1 = {'synonyms': metabolite}
            query = {'$or': [con_0, con_1]}
            projection = {'name': 1, 'chebi_id': 1, 'kegg_id': 1}
            doc = self.meta_collection.find_one(filter=query, projection=projection, collation=self.collation)
            if doc is not None:
                metabolite_name = doc.get('name')
                chebi_id = doc.get('chebi_id')
                kegg_id = doc.get('kegg_id')
                inchikey = doc.get('InChI_Key')
            else:
                metabolite_name = metabolite
                chebi_id = None
                kegg_id = None
                inchikey = None
            row['reference'] = reference
            row['unit'] = unit
            row['last_modified'] = datetime.datetime.utcnow()
            try:
                row['ncbi_taxonomy_id'] = self.taxon_collection.find_one({'tax_name': row['species_name']},
                                                                        projection={'tax_id': 1}, collation=self.collation).get('tax_id')
            except AttributeError:
                continue
            obj = row.to_dict()
            obj.pop('metabolite')
            obj.pop('abbreviation')
            self.collection.update_one({'metabolite': metabolite_name},
                                       {'$addToSet': {'concentrations': obj,
                                                      'synonyms': row['abbreviation']},
                                        '$set': {'chebi_id': chebi_id,
                                                 'kegg_id': kegg_id,
                                                 'inchikey': inchikey}},
                                       collation=self.collation, upsert=True)

    def fill_gerosa_collection(self, file_location, sheet_name='Sheet1', start_row=0,
                                use_columns='A:E,G', column_names=['EC_Number', 'metabolite',
                                'concentration', 'k_m', 'abbreviation', 'species_name'],
                                reference={'doi': '10.1038/nchembio.2077'}, unit='M',
                                nrows=None):
        """
        Fill collection using primary source http://dx.doi.org/10.1016/j.cels.2015.09.008.

        Args:
            sheet_name(:obj:`str`, optional): Name of sheet in excel.
            start_row (:obj:`int`, optional): Read from csv row. Defaults to 0.
            use_columns(:obj:`str`): Indicates comma separated list of Excel column letters and column ranges (e.g. "A:E" or "A,C,E:F"). Ranges are inclusive of both sides.
            column_names(:obj:`list` of :obj:`str`): Names of columns used.
            reference(:obj:`Obj`): reference information.
            unit(:obj:`str`): Unit used for Km.
        """
        df = pd.read_excel(file_location, sheet_name=sheet_name,
                         header=0, usecols=use_columns,
                         skiprows=start_row, nrows=nrows)
        df.columns = [x.lower() for x in column_names]
        row_count = len(df.index)
        for i, row in df.iterrows():
            if i == self.max_entries:
                break
            if i % 10 == 0 and self.verbose:
                print("Processing locus {} out {}".format(i, row_count))
            concentrations = []            
            metabolite = row['metabolite']
            con_0 = {'name': metabolite}
            con_1 = {'synonyms': metabolite}
            query = {'$or': [con_0, con_1]}
            projection = {'name': 1, 'chebi_id': 1, 'kegg_id': 1}
            doc = self.meta_collection.find_one(filter=query, projection=projection, collation=self.collation)
            if doc is not None:
                metabolite_name = doc.get('name')
            else:
                metabolite_name = metabolite

            obj = row.to_dict()
            for i, (_key, val) in enumerate(obj.items()):
                if i != 0 and i < 9:
                    concentrations.append({'reference': reference,
                                           'unit': unit,
                                           'last_modified': datetime.datetime.utcnow(),
                                           'ncbi_taxonomy_id': 562,
                                           'species_name': 'Escherichia coli',
                                           'nutrient': _key,
                                        #    'time_unit': 'h',
                                           'concentration': val,
                                           'std': obj[_key+'_std']})
            self.collection.update_one({'metabolite': metabolite_name},
                                       {'$addToSet': {'concentrations': {'$each': concentrations}}},
                                       collation=self.collation, upsert=True)

    def fill_meta(self):
        """Fill exisiting metabolites with meta information.
        """
        query = {'inchikey': {'$exists': False}}
        count = self.collection.count_documents(query)
        for i, doc in enumerate(self.collection.find(filter=query)):
            if self.verbose and i % 20 == 0:
                print('Processing doc {} out of {}'.format(i, count))
            name = doc['metabolite']
            names = doc.get('synonyms', []) + [name]
            meta_doc = self.get_doc_by_name(names)
            if meta_doc is None:
                print('    '  + str(name))
            else:
                self.collection.update_one({'metabolite': name},
                                            {'$set': {'inchikey': meta_doc['InChI_Key'],
                                                    'kegg_id': meta_doc['kegg_id'],
                                                    'chebi_id': meta_doc['chebi_id']}})


import datanator.config.core
from pathlib import Path

def main():
    db = 'datanator'
    collection_str = 'metabolite_concentrations'
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    MongoDB = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']
    manager = Concentration(MongoDB=MongoDB, db=db, collection_str=collection_str,
                            username=username, password=password)

    # status = manager.merge_duplicates()
    # print(status)

    # file_location = Path('~/karr_lab/datanator/docs/metabolite_concentration/41589_2016_BFnchembio2077_MOESM585_ESM.xlsx').expanduser()
    # manager.fill_collection(file_location, start_row=[0])
    # file_location = Path('~/karr_lab/datanator/docs/metabolite_concentration/41589_2016_BFnchembio2077_MOESM586_ESM.xlsx').expanduser()
    # manager.fill_collection(file_location, start_row=[0])

    # manager.remove_eymdb()
    # manager.grab_eymdb()

    # column_names=['Metabolite', 'Acetate', 'Fructose', 'Galactose', 'Glucose', 'Glycerol', 'Gluconate', 'Pyruvate', 'Succinate',
    #               'Acetate_std', 'Fructose_std', 'Galactose_std', 'Glucose_std', 'Glycerol_std', 'Gluconate_std', 'Pyruvate_std', 'Succinate_std']
    # file_location = Path('~/karr_lab/datanator/docs/metabolite_concentration/mmc2.xlsx').expanduser()
    # manager.fill_gerosa_collection(file_location, start_row=[0], column_names=column_names, unit='µmol * gCDW-1',
    #                                 reference={'doi': '10.1016/j.cels.2015.09.008'},
    #                                 sheet_name='Metabolite concentrations', use_columns='C:K,M:T', nrows=43)

    # column_names=['Metabolite', '-0.766', '-0.366', '0.083', '0.416', '0.750', '1.083', '1.516', '1.866', '2.200',
    #               '-0.766_std', '-0.366_std', '0.083_std', '0.416_std', '0.750_std', '1.083_std', '1.516_std', '1.866_std', '2.200_std']
    # file_location = Path('~/karr_lab/datanator/docs/metabolite_concentration/mmc3.xlsx').expanduser()
    # manager.fill_gerosa_collection(file_location, start_row=[0,1], column_names=column_names, unit='µmol * gCDW-1',
    #                                 reference={'doi': '10.1016/j.cels.2015.09.008'},
    #                                 sheet_name='Metabolite concentrations', use_columns='B:K,M:U', nrows=21)

    # manager.fill_meta()

    manager.edit_nchembio2077()

if __name__ == '__main__':
    main()