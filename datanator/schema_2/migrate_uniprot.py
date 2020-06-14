from datanator_query_python.util import mongo_util
from datanator_query_python.config import config
from datanator_query_python.query import query_taxon_tree


class MigrateUniprot(mongo_util.MongoUtil):

    def __init__(self, readPreference='primary', collection='uniprot',
                max_entries=float('inf'), username=config.Config.USERNAME,
                password=config.Config.PASSWORD):
        super().__init__(readPreference=readPreference, username=username,
                        password=password)
        self.collection = collection
        self.from_collection = self.client.get_database('datanator')[collection]
        self.to_collection = self.client.get_database('datanator-test')[collection]
        self.max_entries = max_entries

    def index_primary(self, _key):
        """Index key (single key ascending)

        Args:
            _key(:obj:`str`): Name of key to be indexed
        """
        self.to_collection.create_index(_key)
    
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
    
    def process_cursor(self, docs):
        """Process mongodb cursor (for parallel processing)
        Transform data and move to new database

        Args:
            docs(:obj:`pymongo.Cursor`): documents to be processed
        """
        for i, doc in enumerate(docs):
            if i == self.max_entries:
                break
            doc["add_id"] = [{"name_space": "gene_name_alt", "value": doc.get("gene_name_alt")},
                             {"name_space": "gene_name_orf", "value": doc.get("gene_name_orf")},
                             {"name_space": "gene_name_oln", "value": doc.get("gene_name_oln")}]
            doc.pop('gene_name_alt', None)
            doc.pop('gene_name_orf', None)
            doc.pop('gene_name_oln', None)
        pass