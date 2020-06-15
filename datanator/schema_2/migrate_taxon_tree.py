from datanator_query_python.config import motor_client_manager
import asyncio
import simplejson as json
from pymongo import UpdateOne
from pymongo.errors import BulkWriteError
from pprint import pprint


class MigrateTaxon:

    def __init__(self, collection="taxon_tree", to_database="datanator-test",
                 from_database="datanator", max_entries=float("inf")):
        self.collection = collection
        self.to_database = to_database
        self.from_collection = motor_client_manager.client.get_database(from_database)[collection]
        self.to_collection = motor_client_manager.client.get_database(to_database)[collection]
        self.max_entries = max_entries

    async def index_primary(self, _key, background=True):
        """Index key (single key ascending)

        Args:
            _key(:obj:`str`): Name of key to be indexed.
            background(:obj:`bool`): Building index in the background.
        """
        yield self.to_collection.create_index(_key, background=background)

    async def get_rank(self, ids, names):
        ''' Given a list of taxon ids, return
            the list of ranks. no rank = '+'

        Args:
            ids(:obj:`list` of :obj:`int`): list of taxon ids [1234,2453,431]
            names(:obj:`list` of :obj:`str`): list of taxon names.

        Return:
            (:obj:`tuple`): canon_anc_id, canon_anc_name
        '''
        canon_anc_id = []
        canon_anc_name = []
        roi = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom']
        projection = {'rank': 1}
        for _id, name in zip(ids, names):
            if _id == 131567:
                canon_anc_id.append(_id)
                canon_anc_name.append(name)
                continue
            query = {'tax_id': _id}
            doc = await self.from_collection.find_one(filter=query, projection=projection)
            rank = doc.get('rank', None)
            if rank in roi:
                canon_anc_id.append(_id)
                canon_anc_name.append(name)
        return canon_anc_id, canon_anc_name

    async def process_cursors(self, skip=0): # docs, l, counter
        """Process mongodb cursor (for parallel processing)
        Transform data and move to new database

        Args:
            docs(:obj:`pymongo.Cursor`): documents to be processed
            l(:obj:`multiprocess.Lock`): order verbose message.
            counter(:obj:`Counter`): to keep track of processor progress across parallel processes.
        """
        bulk_write = []
        if self.max_entries == float("inf"):
            limit = 0
        else:
            limit = self.max_entries
        docs = self.from_collection.find(filter={}, projection={'_id': 0},
                                      no_cursor_timeout=True, batch_size=1000,
                                      limit=limit)
        i = 0
        async for doc in docs:
            i += 1          
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
            canon_anc_ids, canon_anc_names = await self.get_rank(anc_ids, anc_names)
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


def main():
    src = MigrateTaxon(to_database="test")
    src.index_primary('tax_id')
    # cursors = src.split_cursors(1, 'datanator', 'taxon_tree')    
    # lock = Lock()
    # counter = Counter(0)    
    # procs = [Process(target=src.process_cursors, args=(docs, lock, counter)) for docs in cursors]
    # for p in procs: p.start()
    # for p in procs: p.join()
    loop = asyncio.get_event_loop()
    loop.run_until_complete(src.process_cursors())
    

if __name__ == '__main__':
    main()