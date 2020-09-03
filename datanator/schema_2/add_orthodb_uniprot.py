from datanator.util import x_ref
from datanator_query_python.config import config


class AddOrtho(x_ref.XRef):
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
                         password=password,
                         des_col=des_col,
                         verbose=verbose,
                         max_entries=max_entries)
        self.collection = self.db_obj[des_col]
        self.max_entries = max_entries
        self.verbose = verbose

    def add_ortho(self, skip=0):
        """Add OrthoDB to existing uniprot entries.

        Args:
            (:obj:`int`, optional): Skipping for x number of records.
        """
        docs = self.collection.find({},
                                    projection={"_id": 0, "uniprot_id": 1},
                                    skip=skip, 
                                    no_cursor_timeout=True)
        count = self.collection.count_documents({})
        for i, doc in enumerate(docs):
            if i == self.max_entries:
                print("Done!")
                break
            if self.verbose and i % 500 == 0:
                print("Processing doc {} out of {} ...".format(i+skip, count))
            uniprot_id = doc["uniprot_id"]
            obj, _ = self.uniprot_id_to_orthodb(uniprot_id)
            self.collection.update_one({"uniprot_id": uniprot_id},
                                       {"$set": {"orthodb_id": obj["orthodb_id"],
                                                 "orthodb_name": obj["orthodb_name"]}})
        print("Done!")

def main():
    conf = config.DatanatorAdmin()
    src = AddOrtho(MongoDB=conf.SERVER,
                    db="datanator-test",
                    des_col="uniprot",
                    username=conf.USERNAME,
                    password=conf.PASSWORD,
                    verbose=True)
    src.add_ortho(skip=95000)


if __name__ == "__main__":
    main()