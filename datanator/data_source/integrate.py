import pymongo


class Integrate():
    def __init__(self, cache_dirname, MongoDB, db, verbose=True, max_entries=float('inf')):
        self.e_coli = 562  # Taxon ID For E. Coli
        self.s_cere = 4932  # Taxon ID for Saccharomyces cerevisiae
        self.cache_dirname = cache_dirname
        self.MongoDB = MongoDB
        self.db = db
        self.verbose = verbose
        self.max_entries = max_entries

    '''Making connection to the database
	'''

    def con_db(self):
        try:
            client = pymongo.MongoClient(
                self.MongoDB, 400)  # 400ms max timeout
            client.server_info()
            db = client[self.db]
            return (client, db)
        except pymongo.errors.ConnectionFailure:
            return ('Server not available')

    '''Retrieve collections with the specified Taxon ID
    '''

    def retrieve_sabio(self, collection, taxon):
        return collection.find({'taxon': taxon}).limit(self.max_entries)


    '''Find reaction participant's chebi and return as two lists
    	(substrate, product)
    '''

    def find_sabio_chebi(self, cursor):
            i = 0
            substrate_chebi = []
            product_chebi = []
            for doc in cursor:
                if self.verbose:
                	print (i)
                substrate_list = doc['reaction_participant'][0]['substrate']
                substrate_chebi += [dic['kegg_id'] for dic in substrate_list]
                product_list = doc['reaction_participant'][1]['product']
                product_chebi += [dic['kegg_id'] for dic in product_list]
                i += 1

                if i > self.max_entries:
                    break

            return (substrate_chebi, product_chebi)
