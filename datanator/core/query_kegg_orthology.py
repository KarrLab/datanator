from datanator.util import mongo_util


class QueryKO:

    def __init__(self, username=None, password=None, server=None, authSource='admin',
                 database='datanator', max_entries=float('inf'), verbose=True):

        mongo_manager = mongo_util.MongoUtil(MongoDB=server, username=username,
                                             password=password, authSource=authSource, db=database)
        self.max_entries = max_entries
        self.verbose = verbose
        self.client, self.db, self.collection = mongo_manager.con_db(
            'kegg_orthology_new')

    def get_ko_by_name(self, name):
        '''
        Get a gene's ko number by its gene name
            Args:
                    name: (:obj: `str`): gene name
            Returns:
                    result: (:obj: `str`): ko number of the gene
        '''
        query = {'gene_name': name}
        projection = {'gene_name': 1, 'kegg_orthology_id': 1}
        collation = {'locale': 'en', 'strength': 2}
        docs = self.collection.find_one(
            filter=query, projection=projection, collation=collation)
        if docs != None:
        	return docs['kegg_orthology_id']
        else:
        	return None

    def get_def_by_ko(self, ko):
        '''
            Get the definition of kegg orthology id
            Args:
                ko (:obj: `list` of :obj: `string`): list of ko id
            Returns:
                result (:obj: `list` of :obj: `dict`): list of information
                [{'ko_number': ..., 'name': [ ...]},
                 {'ko_number': ..., 'name': [ ...]}]
        '''
        result = []
        query = {'kegg_orthology_id': {'$in': ko}}
        projection = {'_id': 0, 'definition.name': 1, 'kegg_orthology_id': 1}
        docs = self.collection.find(filter=query, projection=projection)
        for doc in docs:
            ko_number = doc['kegg_orthology_id']
            names = doc['definition']['name']
            dic = {'ko_number': ko_number, 'ko_name': names}
            result.append(dic)
        return result

