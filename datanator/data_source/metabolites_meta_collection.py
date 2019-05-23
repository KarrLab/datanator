from datanator.core import query_nosql


class MetabolitesMeta(query_nosql.DataQuery):

    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet=None, db=None,
                 collection_str=None, verbose=False, max_entries=float('inf')):
        self.cache_dirname = cache_dirname
        self.MongoDB = MongoDB
        self.replicaSet = replicaSet
        self.db = db
        self.collection_str = collection_str
        self.verbose = verbose
        self.max_entries = max_entries

        super(MetabolitesMeta, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, replicaSet=replicaSet,
                                              db=db, collection_str=collection_str, verbose=verbose, max_entries=max_entries)

    def load_content(self):
        pass

    def find_rxn_id(self, inchi=None):
    	'''Find reactions' kinlaw_id in sabio_rk given inchi structures 
    	'''
    	substrate = 'reaction_participant.substrate.structure.inchi'
    	product = 'reaction_participant.product.structure.inchi'
    	c = 'sabio_rk'
    	query = {'$or': [{substrate : {'$regex' : ".*"+inchi+".*"}},
    					 {product: {'$regex' : ".*"+inchi+".*"}} ]}
    	docs = self.doc_feeder( collection_str=c, query=query )

    def find_metabolite_inchi(self, doc):
        '''Find inchi structure information of metabolites in ecmdb or ymdb
        '''
        try:
        	return self.parse_inchi(doc['inchi'])
        except KeyError:
        	return 'No key named "inchi" in given document'

    def parse_inchi(self, inchi=None):
        '''Remove molecules's protonation state
        "InChI=1S/H2O/h1H2" => "InChI=1S/H2O"
        '''
        if self.verbose:
            print('Parsing inchi by taking out protonation state')
        inchi_neutral = inchi.split('/h')[0]
        return inchi_neutral


def main():

    manager_ecmdb = MetabolitesMeta(cache_dirname=None, MongoDB=None, replicaSet=None, db=None,
                                    collection_str='ecmdb', verbose=False, max_entries=float('inf'))

    manager_ecmdb.load_content()


if __name__ == '__main__':
    main()
