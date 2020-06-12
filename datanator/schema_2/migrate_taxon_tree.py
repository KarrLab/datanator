from datanator_query_python.config import config as query_config
from datanator_query_python.util import mongo_util
from multiprocessing import Process, Lock, Value
from datanator_query_python.query import query_taxon_tree
from pymongo import UpdateOne
from pymongo.errors import BulkWriteError
from pprint import pprint
import time


class MigrateTaxon(mongo_util.MongoUtil):

    def __init__(self, 
                server=query_config.AtlasConfig.SERVER,
                username=query_config.SchemaMigration.USERNAME,
                password=query_config.SchemaMigration.PASSWORD,
                authSource=query_config.AtlasConfig.AUTHDB,
                replicaSet=query_config.AtlasConfig.REPLSET,
                readPreference=query_config.AtlasConfig.READ_PREFERENCE,
                collection='taxon_tree',
                to_database='datanator-test',
                from_database='datanator',
                max_entries=float('inf')):
        super().__init__(MongoDB=query_config.AtlasConfig.SERVER,
                        username=query_config.SchemaMigration.USERNAME,
                        password=query_config.SchemaMigration.PASSWORD,
                        authSource=query_config.AtlasConfig.AUTHDB,
                        replicaSet=query_config.AtlasConfig.REPLSET,
                        readPreference=query_config.AtlasConfig.READ_PREFERENCE,
                        db=from_database)
        self.from_collection = self.db_obj[collection]
        self.to_collection = self.client[to_database][collection]
        self.max_entries = max_entries

    def index_primary(self, _key, background=True):
        """Index key (single key ascending)

        Args:
            _key(:obj:`str`): Name of key to be indexed.
            background(:obj:`bool`): Building index in the background.
        """
        self.to_collection.create_index(_key, background=background)

    def get_rank(self, ids):
        ''' Given a list of taxon ids, return
            the list of ranks. no rank = '+'

        Args:
            ids(:obj:`list` of :obj:`int`): list of taxon ids [1234,2453,431]

        Return:
            (:obj:`list` of :obj:`str`): list of ranks ['kingdom', '+', 'phylum']
        '''
        ranks = []
        roi = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom']
        projection = {'rank': 1, 'tax_id': 1}
        for _id in ids:
            query = {'tax_id': _id}
            cursor = self.from_collection.find_one(filter=query, projection=projection)
            rank = cursor.get('rank', None)
            tax_id = cursor.get('tax_id', None)
            if rank in roi:
                ranks.append(rank)
            elif tax_id == 131567:
                ranks.append('cellular organisms')
            else:
                ranks.append('+')
        return ranks

    def process_cursors(self, skip=0): # docs, l, counter
        """Process mongodb cursor (for parallel processing)
        Transform data and move to new database

        Args:
            docs(:obj:`pymongo.Cursor`): documents to be processed
            l(:obj:`multiprocess.Lock`): order verbose message.
            counter(:obj:`Counter`): to keep track of processor progress across parallel processes.
        """
        bulk_write = []
        docs = self.to_collection.find(filter={}, projection={'_id': 0},
                                      no_cursor_timeout=True, batch_size=1000)
        for i, doc in enumerate(docs):
            if i == self.max_entries:
                break            
            if i != 0 and i % 100 == 0:
                print("Processing file {}".format(i))
                try:
                    self.to_collection.bulk_write(bulk_write)
                    bulk_write = []
                except BulkWriteError as bwe:
                    pprint(bwe.details)
                    bulk_write = []
            canon_anc_names = []
            canon_anc_ids = []
            anc_ids = doc['anc_id']
            anc_names = doc['anc_name']
            ranks = self.get_rank(anc_ids)
            [canon_anc_names.append(anc) for (anc, rank) in zip(anc_names, ranks) if rank != '+']
            [canon_anc_ids.append(anc) for (anc, rank) in zip(anc_ids, ranks) if rank != '+']
            doc['canon_anc_ids'] = canon_anc_ids
            doc['canon_anc_names'] = canon_anc_names
            doc['schema_version'] = '2'
            bulk_write.append(UpdateOne({'tax_id': doc['tax_id']}, {'$set': doc}, upsert=True))
        if len(bulk_write) != 0:
            try:
                self.to_collection.bulk_write(bulk_write)
            except BulkWriteError as bwe:
                pprint(bwe.details)
            finally:
                print("Done.")            


class Counter(object):
    def __init__(self, initval=0):
        self.val = Value('i', initval)
        self.lock = Lock()

    def increment(self):
        with self.lock:
            self.val.value += 1

    def value(self):
        with self.lock:
            return self.val.value


def main():
    src = MigrateTaxon()
    src.index_primary('tax_id')
    # cursors = src.split_cursors(1, 'datanator', 'taxon_tree')    
    # lock = Lock()
    # counter = Counter(0)    
    # procs = [Process(target=src.process_cursors, args=(docs, lock, counter)) for docs in cursors]
    # for p in procs: p.start()
    # for p in procs: p.join()
    src.process_cursors()

if __name__ == '__main__':
    main()