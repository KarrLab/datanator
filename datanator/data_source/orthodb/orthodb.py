from datanator_query_python.util import mongo_util
from pymongo.collation import Collation, CollationStrength
from datanator.util import file_util
import csv


class OrthoDB(mongo_util.MongoUtil):
    def __init__(self,
                 MongoDB=None,
                 db=None,
                 des_col=None,
                 username=None,
                 password=None,
                 max_entries=float('inf'),
                 verbose=True):
        super().__init__(MongoDB=MongoDB,
                         db=db,
                         username=username,
                         password=password)
        self.max_entries = max_entries
        self.db = db
        self.collection = self.db_obj[des_col]
        self.verbose = verbose
        self.taxon = self.client["datanator-test"]["taxon_tree"]
        self.file_manager = file_util.FileUtil()
        self.collation = Collation(locale='en', strength=CollationStrength.SECONDARY)

    def pairwise_name_group(self, 
                            url,
                            batch_size=100):
        """Parse file in https://v101.orthodb.org/download/odb10v1_OGs.tab.gz
        into MongoDB.

        Args:
            (:obj:`str`): URL of the file.
            (:obj:`int`): Number of docs to be inserted at once.
        """
        with open(url) as f:
            x = csv.DictReader(f, fieldnames=["_id", "taxon", "name"],
                               delimiter="\t")
            tmp = []
            for i, row in enumerate(x):
                if i == self.max_entries:
                    break
                elif i % batch_size == 0 and i != 0:
                    self.collection.insert_many(tmp)
                    tmp = []
                else:
                    obj = {"orthodb_id": row["_id"],
                           "orthodb_name": row["name"]}
                    tmp.append(obj)
            if len(tmp) != 0:
                self.collection.insert_many(tmp)