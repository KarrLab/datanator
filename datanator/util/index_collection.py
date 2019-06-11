'''
    Index collections in MongoDB accordingly
'''
from datanator.util import mongo_util
import pymongo
from datanator.util import server_config

class IndexCollection(mongo_util.MongoUtil):

    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet=None, db=None,
                 verbose=False, max_entries=float('inf'), username = None, 
                 password = None, authSource = 'admin'):
        self.cache_dirname = cache_dirname
        self.MongoDB = MongoDB
        self.db = db
        self.replicaSet = replicaSet
        self.verbose = verbose
        self.max_entries = max_entries
        super(IndexCollection, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, replicaSet=replicaSet, 
                                    db=db, verbose=verbose, max_entries=max_entries, username = username, 
                                    password = password, authSource = authSource)

    def index_corum(self, collection_str):
        '''Index fields in corum collection
        '''
        if self.verbose:
            print('Indexing corum ...')
        collection = self.fill_db(collection_str)
        index1 = pymongo.IndexModel(
            [("$**", pymongo.TEXT)], background=False, sparse=True)  # index all text fields
        index2 = pymongo.IndexModel(
            [("PubMed ID", pymongo.ASCENDING)], background=False, sparse=True)
        index3 = pymongo.IndexModel(
            [("SWISSPROT organism (NCBI IDs)", pymongo.ASCENDING)], background=False, sparse=True)
        collection.create_indexes([index1, index2, index3])

    def index_sabio(self, collection_str='sabio_rk', collation = {}):
        '''Index relevant fields in sabio_rk collection
        '''
        if self.verbose:
            print('Indexing Sabio RK ... ')
        collection = self.fill_db(collection_str)
        index1 = pymongo.IndexModel(
            [("reaction_participant.substrate.structure.inchi", pymongo.TEXT),
            ("reaction_participant.product.structure.inchi", pymongo.TEXT)],
             background=False, sparse=True)  # index inchi fields
        index2 = pymongo.IndexModel(
            [('kinlaw_id', pymongo.ASCENDING)], background=False, sparse=True)
        index3 = pymongo.IndexModel(
            [('enzymes.subunit.coefficient', pymongo.ASCENDING)], background=False, sparse=True)
        index4 = pymongo.IndexModel([('parameter.sabio_compound_id', pymongo.ASCENDING),
                                     ('parameter.value', pymongo.ASCENDING),
                                     ('parameter.error', pymongo.ASCENDING),
                                     ('parameter.sbo_type', pymongo.ASCENDING)], background=False, sparse=True)
        index6 = pymongo.IndexModel([('reaction_participant.product.sabio_compound_id', pymongo.ASCENDING)],
                                    background=False, sparse=True)
        index5 = pymongo.IndexModel([('reaction_participant.substrate.sabio_compound_id', pymongo.ASCENDING)],
                                    background=False, sparse=True)

        index7 = pymongo.IndexModel(
            [('taxon', pymongo.ASCENDING)], background=False, sparse=False)
        index8 = pymongo.IndexModel(
            [('taxon_wildtype', pymongo.ASCENDING)], background=False, sparse=False)
        index9 = pymongo.IndexModel(
            [('temperature', pymongo.ASCENDING)], background=False, sparse=False)
        index10 = pymongo.IndexModel(
            [('ph', pymongo.ASCENDING)], background=False, sparse=False)

        collection.create_indexes([index1, index2, index3, index4, index5,
                                   index6, index7, index8, index9, index10])

    def index_strdb(self, collection_str='ecmdb'):
        '''Index relevant fields in string only collections:
                ecmdb, ymdb, and intact_interaction
        '''
        if self.verbose:
            print('Indexing {} ...'.format(collection_str))
        collection = self.fill_db(collection_str)
        index1 = pymongo.IndexModel(
            [("$**", pymongo.TEXT)], background=False, sparse=True)
        collection.create_indexes([index1])

    def index_intact_complex(self, collection_str='intact_complex'):
        '''Index intact_complex collection
        '''
        if self.verbose:
            print('Indexing intact_complex ... ')
        collection = self.fill_db(collection_str)
        index1 = pymongo.IndexModel(
            [("$**", pymongo.TEXT)], background=False, sparse=True)
        index2 = pymongo.IndexModel(
            [('ncbi_id', pymongo.ASCENDING)], background=False, sparse=True)
        collection.create_indexes([index1, index2])

    def index_pax(self, collection_str='pax'):
        '''Index Pax collection
        '''
        if self.verbose:
            print('Indexing pax ...')
        collection = self.fill_db(collection_str)
        index1 = pymongo.IndexModel(
            [("$**", pymongo.TEXT)], background=False, sparse=True)
        index2 = pymongo.IndexModel(
            [('ncbi_id', pymongo.ASCENDING)], background=False, sparse=True)
        index3 = pymongo.IndexModel(
            [('weight', pymongo.ASCENDING)], background=False, sparse=True)
        index4 = pymongo.IndexModel(
            [('score', pymongo.ASCENDING)], background=False, sparse=True)
        index5 = pymongo.IndexModel(
            [('coverage', pymongo.ASCENDING)], background=False, sparse=True)
        collection.create_indexes([index1, index2, index3, index4, index5])

    def index_uniprot(self, collection_str='uniprot'):
        '''Index uniprot collection
        '''
        if self.verbose:
            print('Indexing uniprot ...')
        collection = self.fill_db(collection_str)
        index1 = pymongo.IndexModel(
            [("$**", pymongo.TEXT)], background=False, sparse=True)
        index2 = pymongo.IndexModel(
            [('length', pymongo.ASCENDING)], background=False, sparse=True)
        collection.create_indexes([index1, index2])


def main():
    config_file = '/root/host/karr_lab/datanator/.config/config_test.ini'
    db = 'datanator'
    username, MongoDB, password, port = server_util.ServerUtil(
            config_file=config_file, username = username, password = password,
                port = port)
    manager = IndexCollection(cache_dirname=None, MongoDB=MongoDB, db=db,
                              verbose=True, max_entries=float('inf'), username = username,
                              password = password)

    manager.index_sabio()
    # manager.index_strdb('ecmdb')
    # manager.index_strdb('ymdb')


if __name__ == '__main__':
    main()
