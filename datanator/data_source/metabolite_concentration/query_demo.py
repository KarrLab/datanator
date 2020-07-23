from datanator_query_python.util import mongo_util
from datanator_query_python.config import config


class QueryDemo(mongo_util.MongoUtil):
    def __init__(self, MongoDB=None,
                       db=None,
                       collection_str=None,
                       password=None,
                       username=None,
                       max_entries=20):
        super().__init__(MongoDB=MongoDB,
                         db=db,
                         username=username,
                         password=password)
        self.collection = self.db_obj[collection_str]
        self.max_entries = max_entries

    def get_canon_ancestors(self, tax_id):
        """Getting canon ancestor information by taxon ID

        Args:
            (:obj:`int`): Taxon ID of organism.

        Return:
            (:obj:`list` of :obj:`Obj`)
        """
        query = {"tax_id": tax_id}
        projection = {"canon_anc_ids": 1, "canon_anc_names": 1,
                      "_id": 0}
        doc = self.collection.find_one(filter=query,
                                       projection=projection)
        result = []
        if doc is None:
            return result
        for _id, name in zip(doc["canon_anc_ids"], doc["canon_anc_names"]):
            obj = {"ncbi_taxonomy_id": _id,
                   "name": name}
            result.append(obj)
        return result

    def demo_find(self, tax_id):
        """Find organism with canon ancestor tax_id

        Args:
            (:obj:`int`): Ancestor ID.

        Return:
            (:obj:`list`)
        """
        query = {"canon_anc_ids": tax_id}
        projection = {"canon_anc_ids": 1, "canon_anc_names": 1,
                      "_id": 0}
        docs = self.collection.find(filter=query,
                                    projection=projection)  
        result = []  
        if docs is None:
            return result
        for i, doc in enumerate(docs):  
            if i == self.max_entries:
                break
            result.append(doc)
        return result
