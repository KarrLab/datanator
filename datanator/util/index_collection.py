'''
	Index collections in MongoDB accordingly
'''
from datanator.util import mongo_util
import pymongo

class IndexCollection(mongo_util.MongoUtil):

    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet=None, db=None,
                 verbose=False, max_entries=float('inf')):
        self.cache_dirname = cache_dirname
        self.MongoDB = MongoDB
        self.db = db
        self.replicaSet = replicaSet
        self.verbose = verbose
        self.max_entries = max_entries

    def index_corum(self, collection_str):
        '''Index fields in corum collection
        '''
        collection = self.fill_db(collection_str)
        index1 = pymongo.IndexModel( [("$**", pymongo.TEXT)] , background=False, sparse=True) #index all text fields
        index2 = pymongo.IndexModel( [("PubMed ID", pymongo.ASCENDING)] , background=False, sparse=True)
        index3 = pymongo.IndexModel( [("SWISSPROT organism (NCBI IDs)", pymongo.ASCENDING)] , background=False, sparse=True)
        collection.create_indexes([index1, index2, index3])

    def index_sabio(self, collection_str):
        '''Index relevant fields in sabio_rk collection
        '''

        collection = self.fill_db(collection_str)
        index1 = pymongo.IndexModel( [("$**", pymongo.TEXT)] , background=False, sparse=True) #index all text fields
        index2 = pymongo.IndexModel( [('kinlaw_id', pymongo.ASCENDING)], background=False, sparse=True)
        index3 = pymongo.IndexModel([('enzymes.subunit.coefficient', pymongo.ASCENDING)], background=False, sparse=True)
        index4 = pymongo.IndexModel([('parameter.sabio_compound_id', pymongo.ASCENDING), 
                                     ('parameter.value', pymongo.ASCENDING),
                                     ('parameter.error', pymongo.ASCENDING),
                                     ('parameter.sbo_type', pymongo.ASCENDING),
                                     ('parameter.observed_value', pymongo.ASCENDING),
                                     ('parameter.observed_error', pymongo.ASCENDING)], background=False, sparse=True)
        index5 = pymongo.IndexModel([('reaction_participant.substrate.sabio_compound_id', pymongo.ASCENDING)],
                                    background=False, sparse=False)
        index6 = pymongo.IndexModel([('reaction_participant.product.sabio_compound_id', pymongo.ASCENDING)],
                                    background=False, sparse=False)
        index7 = pymongo.IndexModel([('taxon',pymongo.ASCENDING)], background=False,sparse=False)
        index8 = pymongo.IndexModel([('taxon_wildtype',pymongo.ASCENDING)], background=False,sparse=False)
        index9 = pymongo.IndexModel([('temperature',pymongo.ASCENDING)], background=False,sparse=False)
        index10 = pymongo.IndexModel([('ph',pymongo.ASCENDING)], background=False,sparse=False)

        collection.create_indexes([index1, index2, index3, index4, index5,
                                  index6, index7, index8, index9, index10])
    
    def index_strdb(self,collection_str):
        '''Index relevant fields in string only collections:
                ECMDB, YMDB, and intact_interaction
        '''
        collection = self.fill_db(collection_str)
        index1 = pymongo.IndexModel( [("$**", pymongo.TEXT)] , background=False, sparse=True)
        collection.create_indexes([index1])

    def index_intact_complex(self,collection_str):
        '''Index intact_complex collection
        '''
        collection = self.fill_db(collection_str)
        index1 = pymongo.IndexModel( [("$**", pymongo.TEXT)] , background=False, sparse=True)
        index2 = pymongo.IndexModel( [('ncbi_id', pymongo.ASCENDING)], background=False, sparse=True)
        collection.create_indexes([index1, index2])

    def index_pax(self,collection_str):
        '''Index Pax collection
        '''
        collection = self.fill_db(collection_str)
        index1 = pymongo.IndexModel( [("$**", pymongo.TEXT)] , background=False, sparse=True)
        index2 = pymongo.IndexModel([('ncbi_id', pymongo.ASCENDING)], background=False, sparse=True)
        index3 = pymongo.IndexModel([('weight', pymongo.ASCENDING)], background=False, sparse=True)
        index4 = pymongo.IndexModel([('score', pymongo.ASCENDING)], background=False, sparse=True)
        index5 = pymongo.IndexModel([('coverage', pymongo.ASCENDING)], background=False, sparse=True)
        collection.create_indexes([index1, index2, index3, index4, index5])

    def index_uniprot(self,collection_str):
    	'''Index uniprot collection
    	'''
    	collection = self.fill_db(collection_str)
    	index1 = pymongo.IndexModel( [("$**", pymongo.TEXT)] , background=False, sparse=True)
    	index2 = pymongo.IndexModel([('length', pymongo.ASCENDING)], background=False, sparse=True)
    	collection.create_indexes([index1, index2])

