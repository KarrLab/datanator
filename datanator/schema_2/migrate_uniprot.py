from datanator_query_python.util import mongo_util
from datanator_query_python.config import config
import simplejson as json
from pymongo import UpdateOne
from pymongo.errors import BulkWriteError
from pprint import pprint


class MigrateUniprot(mongo_util.MongoUtil):

    def __init__(self, readPreference='primary',
                username=config.SchemaMigration.USERNAME,
                password=config.SchemaMigration.PASSWORD,
                collection='uniprot',
                to_database='datanator-test',
                from_database='datanator',
                max_entries=float('inf')):
        super().__init__(readPreference=readPreference, username=username,
                        password=password)
        self.collection = collection
        self.to_database = to_database
        self.from_collection = self.client.get_database(from_database)[collection]
        self.to_collection = self.client.get_database(to_database)[collection]
        self.max_entries = max_entries

    def index_primary(self, _key, background=True):
        """Index key (single key ascending)

        Args:
            _key(:obj:`str`): Name of key to be indexed
        """
        self.to_collection.create_index(_key, background=background)
    
    def get_multi_cursor(self, count):
        """Manually return multiple cursors.

        Args:
            count(:obj:`int`): Number of cursors wanted.

        Return:
            (:obj:`list` of :obj:`pymongo.Cursor`)
        """
        query = {}
        total_docs = self.from_collection.count_documents(query)
        step_size = int((total_docs % count + total_docs) / count)
        result = []
        for i in range(0, count):
            result.append(self.from_collection.find(filter=query, skip=i * step_size, limit=step_size))
        return result

    def get_rank(self, ids, names):
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
            doc = self.from_collection.find_one(filter=query, projection=projection)
            rank = doc.get('rank', None)
            if rank in roi:
                canon_anc_id.append(_id)
                canon_anc_name.append(name)
        return canon_anc_id, canon_anc_name
    
    def process_cursor(self, skip=0):
        """Process mongodb cursor
        Transform data and move to new database

        Args:
            docs(:obj:`pymongo.Cursor`): documents to be processed
        """
        bulk_write = []
        query = {"modifications": {"$exists": True}}
        docs = self.from_collection.find(filter=query, projection={'_id': 0},
                                      no_cursor_timeout=True, batch_size=500,
                                      skip=skip)
        for i, doc in enumerate(docs):
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
            uniprot_id = doc.get('uniprot_id')
            doc["add_id"] = [{"name_space": "gene_name_alt", "value": doc.get("gene_name_alt")},
                             {"name_space": "gene_name_orf", "value": doc.get("gene_name_orf")},
                             {"name_space": "gene_name_oln", "value": doc.get("gene_name_oln")}]
            doc.pop('gene_name_alt', None)
            doc.pop('gene_name_orf', None)
            doc.pop('gene_name_oln', None)
            doc['schema_version'] = "2"
            doc['canon_anc_names'] = []  # place holder
            doc['canon_anc_ids'] = [] # place holder
            modifications = doc.get('modifications')
            if modifications is not None:
                bw = []
                for mod in modifications:
                    mod['uniprot_id'] = uniprot_id
                    mod['schema_version'] = "2"
                    reference = doc['reference']
                    doc['reference'] = {"namespace": "doi", "value": reference["doi"]}
                    bw.append(json.loads(json.dumps(mod, ignore_nan=True)))  
                self.client.get_database(self.to_database)['protein_modifications'].insert_many(bw)
            doc.pop('modifications', None)
            bulk_write.append(UpdateOne({'uniprot_id': uniprot_id}, {'$set': json.loads(json.dumps(doc, ignore_nan=True))}, upsert=True))
        if len(bulk_write) != 0:
            try:
                self.to_collection.bulk_write(bulk_write)
            except BulkWriteError as bwe:
                pprint(bwe.details)
            finally:
                print("Done.")   


def main():
    src = MigrateUniprot(max_entries=101)
    src.index_primary('uniprot_id')
    src.process_cursor()

if __name__ == '__main__':
    main()