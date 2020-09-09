from datanator.util import x_ref
from datanator_query_python.config import config
import requests
import pandas as pd
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
            dic = {}
            count = 0   # number of orthodb genes
            uniprot_doc = 0  # number of uniprot doc processed
            for i, row in enumerate(x):
                if i == self.max_entries:
                    break
                elif i < skip or row is None:
                    continue
                orthodb_gene = row[0]
                _id = row[1]
                name = row[2]
                if self.verbose and i % 500 == 0:
                    print("Process row {} with orthodb gene ID {} ...".format(i, orthodb_gene)) 
                if dic.get(orthodb_gene) is None:               # new orthodb gene
                    count += 1
                    dic[orthodb_gene] = {name: [_id]}
                else:
                    if dic[orthodb_gene].get(name) is None:     # new namespace
                        dic[orthodb_gene][name] = [_id]
                    else:
                        dic[orthodb_gene][name].append(_id)
                
                if count % batch_size == 0 and self.verbose:
                    for key, val in dic.items():
                        if val.get("UniProt") is None:
                            continue
                        else:
                            add_ids = []                            
                            uniprot_doc += 1
                            for k, v in val.items():
                                if k != "UniProt":
                                    for l in v:
                                        add_ids.append({"namespace": k,
                                                        "value": l})
                                else:
                                    add_ids.append({"namespace": "orthodb_gene",
                                                    "value": key})
                                    uniprot_id = v[0]
                            self.collection.update_one({"uniprot_id": uniprot_id},
                                                       {"$addToSet": {"add_id": {"$each": add_ids}}},
                                                       upsert=False)
                            if uniprot_doc % 100 == 0 and self.verbose:
                                print("     Processing uniprot doc {}... with uniprot_id {}".format(uniprot_doc, uniprot_id))
                    dic = {}
                else:
                    continue
                    

                # tmp = requests.get("https://www.orthodb.org/search?query={}&limit=1".format(orthodb_gene)).json()["data"]
                # if len(tmp) == 0:
                #     continue
                # else:
                #     orthodb_id = tmp[0]               
 
                # doc = self.orthodb.find_one({"orthodb_id": orthodb_id})
                # if doc is None:
                #     r = requests.get("https://dev.orthodb.org/group?id={}".format(orthodb_id)).json()
                #     orthodb_name = r["data"]["name"]
                #     self.orthodb.insert_one({"orthodb_id": orthodb_id,
                #                              "orthodb_name": orthodb_name})
                # else:
                #     orthodb_name = doc["orthodb_name"]


    def display_tab(self, _file, 
                    skip=0, batch_size=100):
        """Display rows of tab files.

        Args:
            _file (:obj:`str`): Location of tab file.
        """
        with open(_file) as f:
            x = csv.reader(f,
                           delimiter="\t")
            for i, row in enumerate(x):
                if i == self.max_entries + skip:
                    break
                elif i < skip or row is None:
                    continue
                print(row)

    def parse_og2_genes(self, _file, skip=0, batch_size=500):
        """https://v101.orthodb.org/download/odb10v1_OG2genes.tab.gz
        Only stores top level group {2, 2157, 2759, 10239}

        Args:
            _file (:obj:`str`): Location of tab file.
        """
        with open(_file) as f:
            x = csv.reader(f,
                           delimiter="\t")
            count = 0
            batch = []
            batch_count = 0
            self.collection.create_index("orthodb_gene")
            for i, row in enumerate(x):
                if i == self.max_entries + skip:
                    break
                elif i < skip or row is None:
                    continue
                og = row[0]
                gene = row[1]
                if self.verbose and i % 100 == 0:
                    print("Processing row {} with gene name {} and og {} ...".format(i, gene, og))

                if og.endswith("at2") or og.endswith("at2157") or og.endswith("at2759") or og.endswith("at10239"):
                    count += 1
                    # if self.verbose and count % 100 == 0:
                    #     print("     Updating doc {} ...".format(count))
                    batch.append({"orthodb_gene": gene,
                                  "top_level_group": og})
                    # self.collection.update_one({"orthodb_gene": gene},
                    #                             {"$set": {"top_level_group": og}},
                    #                             upsert=True)
                    if self.verbose and count % batch_size == 0:
                        batch_count += 1
                        print("     Inserting batch {} ...".format(batch_count))
                        self.collection.insert_many(batch)
                        batch = []
                else:
                    continue

        if len(batch) != 0:
            self.collection.insert_many(batch)
            print("Done.")


    def add_uniprot(self, 
                url,
                skip=0):
        """Add uniprot_id to docs in orthodb_gene collection

        Args:
            url(:obj:`str`): URL of the file.
            skip(:obj:`int`, optional): First number of rows to skip.
        """
        with open(url) as f:
            x = csv.reader(f,
                           delimiter="\t")
            for i, row in enumerate(x):
                if i == self.max_entries:
                    break
                elif i < skip or row is None:
                    continue
                orthodb_gene = row[0]
                _id = row[1]
                name = row[2]
                if self.verbose and i % 500 == 0:
                    print("Process row {}...".format(i))
                if name != "UniProt": 
                    continue
                else:
                    self.collection.update_one({"orthodb_gene": orthodb_gene},
                                               {"$set": {"uniprot_id": _id}},
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

    # # add x ref to uniprot collection
    # db = "datanator-test"
    # des_col = "uniprot"
    # src = AddOrtho(MongoDB=conf.SERVER,
    #                 db=db,
    #                 des_col=des_col,
    #                 username=conf.USERNAME,
    #                 password=conf.PASSWORD,
    #                 verbose=True)
    # src.add_x_ref_uniprot('./docs/orthodb/odb10v1_gene_xrefs.tab',
    #                       skip=66500)

    # # add group-gene pairing in new collection.
    # db = "datanator"
    # des_col = "orthodb_gene"
    # src = AddOrtho(MongoDB=conf.SERVER,
    #                 db=db,
    #                 des_col=des_col,
    #                 username=conf.USERNAME,
    #                 password=conf.PASSWORD,
    #                 verbose=True)
    # src.parse_og2_genes('./docs/orthodb/odb10v1_OG2genes.tab',
    #                       skip=4418100)

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

    # add uniprot_id to orthodb_gene collection
    db = "datanator"
    des_col = "orthodb_gene"
    src = AddOrtho(MongoDB=conf.SERVER,
                    db=db,
                    des_col=des_col,
                    username=conf.USERNAME,
                    password=conf.PASSWORD,
                    verbose=True)
    src.add_uniprot('./docs/orthodb/odb10v1_gene_xrefs.tab',
                    skip=0)


if __name__ == "__main__":
    main()