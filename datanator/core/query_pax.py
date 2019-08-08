from datanator_flask_REST.util import chem_util, file_util
from . import query_nosql


class QueryPax(query_nosql.DataQuery):
    '''Queries specific to pax collection
    '''

    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet=None, db='datanator',
                 collection_str='pax', verbose=False, max_entries=float('inf'), username=None,
                 password=None, authSource='admin'):
        super(query_nosql.DataQuery, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB,
                                                    replicaSet=replicaSet, db=db,
                                                    verbose=verbose, max_entries=max_entries, username=username,
                                                    password=password, authSource=authSource)
        self.chem_manager = chem_util.ChemUtil()
        self.file_manager = file_util.FileUtil()
        self.client, self.db_obj, self.collection = self.con_db(collection_str)

    def get_all_species(self):
        '''
                Get a list of all species in pax collection
                Returns:
                        results (:obj: `list` of :obj: `str`): list of specie names
                                                    with no duplicates
        '''
        results = []
        query = {}
        projection = {'species_name': 1}
        docs = self.collection.find(filter=query, projection=projection)
        for doc in docs:
            results.append(doc['species_name'])
        return list(set(results))
