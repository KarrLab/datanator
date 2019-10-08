from datanator.util import mongo_util, chem_util, file_util
from datanator_query_python.query import query_nosql
import numpy as np

class QueryMetabolitesMeta(query_nosql.DataQuery):
    '''Queries specific to metabolites_meta collection
    '''

    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet=None, db=None,
                 collection_str='metabolites_meta', verbose=False, max_entries=float('inf'), username=None,
                 password=None, authSource='admin'):
        self.collection_str = collection_str
        self.verbose = verbose
        super(query_nosql.DataQuery, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB,
                                        replicaSet=replicaSet, db=db,
                                        verbose=verbose, max_entries=max_entries, username=username,
                                        password=password, authSource=authSource)
        self.client, self.db_obj, self.collection = self.con_db(
            self.collection_str)
        self.file_manager = file_util.FileUtil()
        self.chem_manager = chem_util.ChemUtil()

    def get_metabolite_synonyms(self, compounds):
        ''' Find synonyms of a compound
            Args:
                compound: name(s) of the compound e.g. "ATP", ["ATP", "Oxygen", ...]
            Returns:
                synonyms: dictionary of synonyms of the compounds
                        {'ATP': [], 'Oxygen': [], ...}
                rxns: dictionary of rxns in which each compound is found
                    {'ATP': [12345,45678,...], 'Oxygen': [...], ...}
        '''
        synonyms = {}
        rxns = {}

        def find_synonyms_of_str(c):
            if len(c) != 0:
                query = {'synonyms': c}
                projection = {'synonyms': 1, '_id': -1, 'kinlaw_id': 1}
                collation = {'locale': 'en', 'strength': 2}
                doc = self.collection.find_one(
                    filter=query, projection=projection, collation=collation)
                synonym = {}
                rxn = {}
                try:
                    synonym[c] = doc['synonyms']
                    rxn[c] = doc['kinlaw_id']
                except TypeError as e:
                    synonym[c] = (c + ' does not exist in ' +
                                  self.collection_str)
                    rxn[c] = (c + ' does not exist in ' + self.collection_str)
                return rxn, synonym
            else:
                return ({'reactions': None}, {'synonyms': None})

        if len(compounds) == 0:
            return ({'reactions': None}, {'synonyms': None})
        elif isinstance(compounds, str):
            rxn, syn = find_synonyms_of_str(compounds)
            synonyms.update(syn)
            rxns.update(rxn)
        else:
            for c in compounds:
                rxn, syn = find_synonyms_of_str(c)
                synonyms.update(syn)
                rxns.update(rxn)
        return rxns, synonyms

    def get_metabolite_inchi(self, compounds):
        '''Given a list of compound name(s)
            Return the corrensponding inchi string
            Args:
                compounds: list of compounds
                ['ATP', '2-Ketobutanoate']
            Return:
                ['....', 'InChI=1S/C4H6O3/c1-2-3(5)4(6)7/...']
        '''
        inchi = []
        projection = {'_id': 0, 'inchi': 1, 'm2m_id': 1, 'ymdb_id': 1}
        collation = {'locale': 'en', 'strength': 2}
        for compound in compounds:
            cursor = self.collection.find_one({'synonyms': compound},
                                              projection=projection, collation=collation)
            inchi.append(
                {"inchi": cursor['inchi'], "m2m_id": cursor.get('m2m_id', None),
                 "ymdb_id": cursor.get('ymdb_id', None)})
        return inchi

    def get_ids_from_hash(self, hashed_inchi):
        ''' Given a hashed inchi string, find its
            corresponding m2m_id and/or ymdb_id
            Args:
                hashed_inchi: string of hashed inchi
            Returns:
                result: dictionary of ids and their keys
                    {'m2m_id': ..., 'ymdb_id': ...}
        '''
        query = {'InChI_Key': hashed_inchi}
        projection = {'_id': 0}
        doc = self.collection.find_one(filter=query, projection=projection)
        result = {}
        result['m2m_id'] = doc.get('m2m_id', None)
        result['ymdb_id'] = doc.get('ymdb_id', None)

        return result

    def get_metabolite_hashed_inchi(self, compounds):
        ''' Given a list of compound name(s)
            Return the corresponding hashed inchi string
            Args:
                compounds: ['ATP', '2-Ketobutanoate']
            Return:
                hashed_inchi: ['3e23df....', '7666ffa....']
        '''
        hashed_inchi = []
        projection = {'_id': 0, 'InChI_Key': 1}
        collation = {'locale': 'en', 'strength': 2}
        for compound in compounds:
            cursor = self.collection.find_one({'synonyms': compound},
                                              projection=projection, collation=collation)
            hashed_inchi.append(cursor['InChI_Key'])
        return hashed_inchi

    def get_metabolite_name_by_hash(self, compounds):
        ''' Given a list of hashed inchi, 
            return a list of name (one of the synonyms)
            for each compound
            Args:
                compounds: list of compounds in inchikey format
            Return:
                result: list of names
                    [name, name, name]
        '''
        result = []
        projection = {'_id': 0, 'synonyms': 1}
        collation = {'locale': 'en', 'strength': 2}
        for compound in compounds:
            cursor = self.collection.find_one({'InChI_Key': compound},
                                              projection=projection)
            if not isinstance(cursor['synonyms'], list):
                cursor['synonyms'] = [cursor['synonyms']]
            result.append(cursor.get('synonyms', ['None']))
            # except TypeError:
            #     result.append(['None'])
        return [x[-1] for x in result]

    def get_metabolite_similar_compounds(self, compounds, num=0, threshold=0):
        ''' Given a list of compound names
            Return the top num number of similar compounds
            with tanimoto score above threshold values
            Args:
                compounds: list of compound names
                num: number of similar compounds to return
                threshold: threshold tanimoto coefficient value
                return_format: return dictionary key format, either
                                hashed inchi or name
            Return:
                result: list of similar compounds and their tanimoto score
                [ {'compound1': score, 'compound2': score, ... 'compound_num': score},
                  {'compound1': score, 'compound2': score, ... 'compound_num': score}, ...]
                    compound(1 - n) will be in name format
                raw: list of similar compounds and their tanimoto score
                [ {'compound1': score, 'compound2': score, ... 'compound_num': score},
                  {'compound1': score, 'compound2': score, ... 'compound_num': score}, ...]
                    compound(1 - n) will be in hashed_inchi format
        '''
        result = []
        raw = []
        hashed_inchi = self.get_metabolite_hashed_inchi(compounds)
        projection = {'_id': 0, 'similar_compounds': 1}

        for item in hashed_inchi:
            cursor = self.collection.find_one(filter={'InChI_Key': item},
                                              projection=projection)
            compounds = cursor['similar_compounds']
            scores = [list(dic.values()) for dic in compounds]
            scores = self.file_manager.unpack_list(scores)
            hashes = [list(dic.keys()) for dic in compounds]
            hashes = self.file_manager.unpack_list(hashes)
            names = self.get_metabolite_name_by_hash(hashes)
            # convert to numpy object for faster calculations
            scores_np = np.asarray(scores)
            indices = np.nonzero(scores_np >= threshold)
            size = indices[0].size
            if size == 0:
                raw.append({'raw': -1})
                result.append({'result': -1})
            elif 0 < size < num:
                first_size = compounds[:size]
                raw = first_size
                replaced = self.file_manager.make_dict(names[:size], scores[:size])
                result.append(replaced)
            else:
                first_num = compounds[:num]
                raw = first_num
                replaced = self.file_manager.make_dict(names[:num], scores[:num])
                result.append(replaced)

        return raw, result
