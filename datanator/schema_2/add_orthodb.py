from datanator.util import x_ref
from datanator_query_python.config import config
import requests
import csv


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
        self.orthodb = self.client["datanator"]["orthodb"]
        self.max_entries = max_entries
        self.verbose = verbose

    def add_ortho(self, skip=0):
        """Add OrthoDB to existing uniprot entries.

        Args:
            (:obj:`int`, optional): Skipping for x number of records.
        """
        query = {"orthodb_id": {"$exists": False}}
        docs = self.collection.find(query,
                                    projection={"_id": 0, "uniprot_id": 1},
                                    skip=skip, 
                                    no_cursor_timeout=True)
        count = self.collection.count_documents(query)
        for i, doc in enumerate(docs):
            if i == self.max_entries:
                print("Done!")
                break
            if self.verbose and i % 500 == 0:
                print("Processing doc {} out of {} ...".format(i+skip, count))
            uniprot_id = doc.get("uniprot_id")
            if uniprot_id is None:
                continue
            obj, _ = self.uniprot_id_to_orthodb(uniprot_id)
            self.collection.update_one({"uniprot_id": uniprot_id},
                                       {"$set": {"orthodb_id": obj["orthodb_id"],
                                                 "orthodb_name": obj["orthodb_name"]}})
        print("Done!")

    def add_x_ref_uniprot(self, 
                url,
                batch_size=100,
                skip=0):
        """Add orthodb gene id to uniprot_id

        Args:
            (:obj:`str`): URL of the file.
            (:obj:`int`, optional): Number of docs to be inserted at once.
            (:obj:`int`, optional): First number of rows to skip.
        """
        with open(url) as f:
            x = csv.reader(f,
                           delimiter="\t")
            count = 0
            for i, row in enumerate(x):
                if i == self.max_entries:
                    break
                elif i < skip or row is None:
                    continue
                orthodb_gene = row[0]
                tmp = requests.get("https://www.orthodb.org/search?query={}&limit=1".format(orthodb_gene)).json()["data"]
                if len(tmp) == 0:
                    continue
                else:
                    orthodb_id = tmp[0]               
                refseq = row[2]
                uniprot_id = row[4] 
                doc = self.orthodb.find_one({"orthodb_id": orthodb_id})
                if doc is None:
                    r = requests.get("https://dev.orthodb.org/group?id={}".format(orthodb_id)).json()
                    orthodb_name = r["data"]["name"]
                    self.orthodb.insert_one({"orthodb_id": orthodb_id,
                                             "orthodb_name": orthodb_name})
                else:
                    orthodb_name = doc["orthodb_name"]
                name = row[7]
                if uniprot_id != "":
                    count += 1
                    if self.verbose and count % 10 == 0:
                        print("Process row {} with orthodb gene ID {}, doc count {}...".format(i, orthodb_gene, count)) 
                    self.collection.update_one({"uniprot_id": uniprot_id},
                                                {"$set": {"orthodb_id": orthodb_id,
                                                          "orthodb_name": orthodb_name},
                                                "$addToSet": {"add_id": {"$each": [{"namespace": "orthodb_gene",
                                                                                    "value": orthodb_gene},
                                                                                    {"namespace": "refseq",
                                                                                    "value": refseq}]}}},
                                                upsert=False)

def main():
    conf = config.DatanatorAdmin()

    # # add to uniprot collection
    # db = "datanator-test"
    # des_col = "uniprot"
    # src = AddOrtho(MongoDB=conf.SERVER,
    #                 db=db,
    #                 des_col=des_col,
    #                 username=conf.USERNAME,
    #                 password=conf.PASSWORD,
    #                 verbose=True)
    # src.add_ortho(skip=422500)

    # add to uniprot collection
    db = "datanator-test"
    des_col = "uniprot"
    src = AddOrtho(MongoDB=conf.SERVER,
                    db=db,
                    des_col=des_col,
                    username=conf.USERNAME,
                    password=conf.PASSWORD,
                    verbose=True)
    src.add_x_ref_uniprot(skip=9000)

    # # add to rna_halflife_new collection
    # db = "datanator"
    # des_col = "rna_halflife_new"
    # src = AddOrtho(MongoDB=conf.SERVER,
    #                 db=db,
    #                 des_col=des_col,
    #                 username=conf.USERNAME,
    #                 password=conf.PASSWORD,
    #                 verbose=True)
    # src.add_ortho(skip=0)


if __name__ == "__main__":
    main()