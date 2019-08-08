from datanator_flask_REST.util import mongo_util, chem_util, file_util
from . import query_nosql
import os
import json

class QueryTaxonTree(query_nosql.DataQuery):
    '''Queries specific to taxon_tree collection
    '''

    def __init__(self, cache_dirname=None, collection_str='taxon_tree', 
                verbose=False, max_entries=float('inf'), username=None, MongoDB=None, 
                password=None, db='datanator', authSource='admin'):
        self.collection_str = collection_str
        super(query_nosql.DataQuery, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB,
                                        db=db, verbose=verbose, max_entries=max_entries, username=username,
                                        password=password, authSource=authSource)
        self.chem_manager = chem_util.ChemUtil()
        self.file_manager = file_util.FileUtil()
        self.client, self.db_obj, self.collection = self.con_db(
            self.collection_str)

    def get_all_species(self):
        ''' Get all organisms in taxon_tree collection
            Return:
                result (:obj: `list` of  :obj: `str`): list of organisms
        '''
        projection = {'tax_name':1}
        mass = self.collection.find({ 'tax_name': {'$exists': True} }, projection=projection)
        
        for thing in mass:
            yield json.dumps({'tax_name': thing['tax_name']})

    def get_name_by_id(self, ids):
        ''' Get organisms' names given their tax_ids
            Args:
                ids: organisms' tax_ids
            Return:
                names: organisms' names
        '''
        names = []
        collation = {'locale': 'en', 'strength': 2}
        projection = {'_id': 0, 'tax_name': 1}
        for _id in ids:
            query = {'tax_id': _id}
            cursor = self.collection.find_one(query, collation=collation,
                                              projection=projection)
            names.append(cursor['tax_name'])
        return names

    def get_anc_by_name(self, names):
        ''' Get organism's ancestor ids by
            using organism's names
            Args:
                names: list of organism's names e.g. Candidatus Diapherotrites
            Return:
                result_id: list of ancestors ids in order of the farthest to the closest
                result_name: list of ancestors' names in order of the farthest to the closest
        '''
        result_id = []
        result_name = []

        collation = {'locale': 'en', 'strength': 2}
        projection = {'_id': 0, 'anc_id': 1, 'anc_name': 1}
        for name in names:
            query = {'tax_name': name}
            cursor = self.collection.find_one(query, collation=collation,
                                              projection=projection)
            result_id.append(cursor['anc_id'])
            result_name.append(cursor['anc_name'])
        return result_id, result_name

    def get_anc_by_id(self, ids):
        ''' Get organism's ancestor ids by
            using organism's ids
            Args:
                ids: list of organism's ids e.g. Candidatus Diapherotrites
            Return:
                result: list of ancestors in order of the farthest to the closest
        '''
        result_name = []
        result_id = []
        projection = {'_id': 0, 'anc_id': 1, 'anc_name': 1}
        for _id in ids:
            query = {'tax_id': _id}
            cursor = self.collection.find_one(query,
                                              projection=projection)
            result_id.append(cursor['anc_id'])
            result_name.append(cursor['anc_name'])
        return result_id, result_name

    def get_common_ancestor(self, org1, org2, org_format='name'):
        ''' Get the closest common ancestor between
            two organisms and their distances to the 
            said ancestor
            Args:
                org1: organism 1
                org2: organism 2
                org_format: the format of organism eg tax_id or tax_name
            Return:
                ancestor: closest common ancestor's name
                distance: each organism's distance to the ancestor
        '''
        if org_format == 'name':
            anc_ids, anc_names = self.get_anc_by_name([org1, org2])
        else:
            anc_ids, anc_names = self.get_anc_by_id([org1, org2])

        if org1 == org2:
            return ('org1', [0, 0])
        org1_anc = anc_ids[0]
        org2_anc = anc_ids[1]

        ancestor = self.file_manager.get_common(org1_anc, org2_anc)
        if ancestor == '':
            return ('No common ancestor', [-1, -1])
        idx_org1 = org1_anc.index(ancestor)
        idx_org2 = org2_anc.index(ancestor)

        distance1 = len(org1_anc) - (idx_org1)
        distance2 = len(org2_anc) - (idx_org2)

        return (ancestor, [distance1, distance2])

    def get_rank(self, ids):
        ''' Given a list of taxon ids, return
            the list of ranks. no rank = '+'
            Args:
                ids: list of taxon ids [1234,2453,431]
            Return:
                ranks: list of ranks ['kingdom', '+', 'phylum']
        '''
        ranks = []
        roi = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
        projection = {'rank': 1}
        for _id in ids:
            query = {'tax_id': _id}
            cursor = self.collection.find_one(filter = query, projection = projection)
            rank = cursor.get('rank', None)
            if rank in roi:
                ranks.append(rank)
            else:
                ranks.append('+')

        return ranks