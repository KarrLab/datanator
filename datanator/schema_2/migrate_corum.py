from datanator_query_python.config import motor_client_manager, config
import simplejson as json
import asyncio
from pymongo import UpdateOne
from pymongo.errors import BulkWriteError
from pprint import pprint
import os


class MigrateCorum:

    def __init__(self, collection="corum", to_database="datanator-test",
                 from_database="datanator", max_entries=float("inf")):
        self.collection = collection
        self.from_database = from_database
        self.to_database = to_database
        self.from_collection = motor_client_manager.client.get_database(from_database)[collection]
        self.to_collection = motor_client_manager.client.get_database(to_database)[collection]
        self.max_entries = max_entries

    async def index_primary(self, _key, background=True):
        """Index key (single key ascending)

        Args:
            _key(:obj:`str`): Name of key to be indexed
        """
        yield self.to_collection.create_index(_key, background=background)

    async def process_cursor(self, skip=0):
        """Transform data and move to new database

        Args:
            docs(:obj:`pymongo.Cursor`): documents to be processed
        """
        bulk_write = []
        query = {}
        if self.max_entries == float('inf'):
            limit = 0
        else:
            limit = self.max_entries
        docs = self.from_collection.find(filter=query, projection={'_id': 0},
                                        no_cursor_timeout=True, batch_size=10,
                                        skip=skip, limit=limit)
        i = 0
        async for doc in docs:
            i += 1
            if i == self.max_entries:
                break
            if i != 0 and i % 50 == 0:
                print("Processing file {}".format(i + skip))
                try:
                    self.to_collection.bulk_write(bulk_write)
                    bulk_write = []
                except BulkWriteError as bwe:
                    pprint(bwe.details)
                    bulk_write = []
            doc.pop("complex_id")
            doc["ncbi_taxonomy_id"] = doc["SWISSPROT_organism_NCBI_ID"]
            doc.pop("SWISSPROT_organism_NCBI_ID")
            doc["schema_version"] = "2"
            bulk_write.append(UpdateOne({'ComplexID': doc.get("ComplexID")}, {'$set': json.loads(json.dumps(doc, ignore_nan=True))}, upsert=True))
        if len(bulk_write) != 0:
            try:
                self.to_collection.bulk_write(bulk_write)
            except BulkWriteError as bwe:
                pprint(bwe.details)
            finally:
                print("Done.")  


def main():
    loop = asyncio.get_event_loop()
    src = MigrateCorum()
    src.index_primary('ComplexID')
    loop.run_until_complete(src.process_cursor(skip=0))

if __name__ == '__main__':
    main()
