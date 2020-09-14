from datanator.util import x_ref
from datanator_query_python.config import config
import requests
import pandas as pd
import csv
from pymongo import UpdateOne
import json


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
        self.uniprot = self.client["datanator-test"]["uniprot"]
        self.max_entries = max_entries
        self.verbose = verbose

    def add_ortho(self, skip=0,
                  batch_size=100):
        """Add OrthoDB to existing uniprot entries, specifically
        orthodb_id and orthodb_name

        Args:
            skip(:obj:`int`, optional): Skipping for x number of records.
            batch_size(:obj:`int`, optional): Bulk write size
        """
        query = {"orthodb_id": {"$exists": False}}
        docs = self.collection.find(query,
                                    projection={"_id": 1, "uniprot_id": 1},
                                    skip=skip, 
                                    no_cursor_timeout=True)
        count = self.collection.count_documents(query)
        bulk = []
        for i, doc in enumerate(docs):
            if i == self.max_entries:
                print("Done!")
                break
            if self.verbose and i % 100 == 0:
                print("Processing doc {} out of {} ...".format(i+skip, count))
            uniprot_id = doc.get("uniprot_id")
            if uniprot_id is None:
                continue
            obj, _ = self.uniprot_id_to_orthodb(uniprot_id)
            bulk.append(UpdateOne({"_id": doc["_id"]},
                                    {"$set": {"orthodb_id": obj["orthodb_id"],
                                                "orthodb_name": obj["orthodb_name"]}}))
            if len(bulk) == batch_size:
                print("     Bulk updating...")
                self.collection.bulk_write(bulk)
                bulk = []
            else:
                continue
        if len(bulk) != 0:
            self.collection.bulk_write(bulk)
        print("Done!")

    def add_x_ref_uniprot(self, 
                url,
                batch_size=100,
                skip=0):
        """Add orthodb gene id to add_id field in uniprot collection

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
                    bulk = []
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
                            bulk.append(UpdateOne({"uniprot_id": uniprot_id},
                                                    {"$addToSet": {"add_id": {"$each": add_ids}}},
                                                    upsert=False))
                            if uniprot_doc % 100 == 0 and self.verbose:
                                print("     Processing uniprot doc {}... with uniprot_id {}".format(uniprot_doc, uniprot_id))
                    if bulk != []:
                        self.collection.bulk_write(bulk)
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
                skip=0,
                batch_size=100):
        """Add uniprot_id to docs in orthodb_gene collection

        Args:
            url(:obj:`str`): URL of the file.
            skip(:obj:`int`, optional): First number of rows to skip.
        """
        bulk = []
        bulk_uniprot = []
        with open(url) as f:
            batch_count = 0
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
                    bulk.append(UpdateOne({"orthodb_gene": orthodb_gene},
                                            {"$set": {"uniprot_id": _id}},
                                            upsert=False))
                    bulk_uniprot.append(UpdateOne({"uniprot_id": _id},
                                                    {"$addToSet": {"add_id": {"namespace": "orthodb_gene",
                                                                              "value": orthodb_gene}}},
                                                    upsert=False))
                    if len(bulk) == batch_size:
                        batch_count += 1
                        print("     Bulk inserting batch {}".format(batch_count))
                        self.collection.bulk_write(bulk)
                        self.uniprot.bulk_write(bulk_uniprot)
                        bulk = []
                        bulk_uniprot = []
                    else:
                        continue
            if len(bulk) != 0:
                self.collection.bulk_write(bulk)
                print("Done.")

    def add_x_ref_rna_halflife(self,
                               skip=0,
                               batch_size=100):
        """Add orthodb_id and orthodb_name to rna_halflives_new collection.

        Args:
            skip (:obj:`int`, optional): First number of docs to  skip. Defaults to 0.
            batch_size (:obj:`int`, optional): Bulk write size. Defaults to 100.
        """
        # con_0 = {"orthodb_id": None}
        con_1 = {"orthodb_id": {"$exists": False}}
        con_2 = {"uniprot_id": {"$ne": None}}
        # comb = {"$or": [con_0, con_1]}
        query = {"$and": [con_1, con_2]}
        docs = self.collection.find(filter=query,
                                    projection={"halflives": 0},
                                    skip=skip,
                                    no_cursor_timeout=True)
        count = self.collection.count_documents(query)
        bulk = []
        for i, doc in enumerate(docs):
            if i == self.max_entries:
                break
            if self.verbose and i % 50 == 0:
                print("Processing doc {} out of {}...".format(i + skip, count))
            uniprot_id = doc["uniprot_id"]
            obj = self.uniprot.find_one({"uniprot_id": uniprot_id},
                                        projection={"orthodb_id": 1,
                                                    "orthodb_name": 1,
                                                    "ncbi_taxonomy_id": 1})
            if obj is not None:
                bulk.append(UpdateOne({"_id": obj["_id"]},
                                      {"$set": {"orthodb_id": obj["orthodb_id"],
                                                "orthodb_name": obj["orthodb_name"]}}))
                if len(bulk) >= batch_size:
                    print("     Bulk insertion...")
                    self.collection.bulk_write(bulk)
                    bulk = []
                continue            
            name = doc["protein_names"][0]
            level = self.uniprot_taxon(uniprot_id)
            if level == []:
                continue
            elif level[1] not in [2, 2157, 2759]:
                level = 10239
            else:
                level = level[1]
            orthodb_id = self.name_level(name, level=level)
            if orthodb_id == "":  # didn't find anything using protein name, use KEGG instead
                orthodb_id = self.name_level(doc["ko_number"], level=level)
                if orthodb_id == "": # Still cannot find anything
                    print("     Cannot find Orthodb group ID for uniprot_id {}".format(uniprot_id))
                    bulk.append(UpdateOne({"_id": doc["_id"]},
                                        {"$set": {"orthodb_id": None,
                                                    "orthodb_name": None}},
                                        upsert=False))
                    continue
            tmp = self.orthodb.find_one({"orthodb_id": orthodb_id})
            if tmp is not None:
                orthodb_name = tmp["orthodb_name"]
            else:   # orthodb_id not in collection orthodb
                orthodb_name = requests.get("https://dev.orthodb.org/group?id={}".format(orthodb_id)).json()["data"]["name"]
                self.orthodb.insert_one({"orthodb_id": orthodb_id,
                                        "orthodb_name": orthodb_name})
            bulk.append(UpdateOne({"_id": doc["_id"]},
                                  {"$set": {"orthodb_id": orthodb_id,
                                            "orthodb_name": orthodb_name}},
                                   upsert=False))
            if len(bulk) >= batch_size:
                print("     Bulk insertion...")
                self.collection.bulk_write(bulk)
                bulk = []
            else:
                continue
        if len(bulk) != 0:
            self.collection.bulk_write(bulk)
        print("Done.")

    def add_x_ref_rna_mod(self,
                          skip=0,
                          batch_size=10):
        """Add orthodb_id and orthodb_name to rna_modification collection.

        Args:
            skip (:obj:`int`, optional): First number of docs to  skip. Defaults to 0.
            batch_size (:obj:`int`, optional): Bulk write size. Defaults to 100.
        """
        con_0 = {"orthodb_id": {"$exists": False}}
        con_1 = {"kegg_orthology_id": {"$exists": True}}
        query = {"$and": [con_0, con_1]}
        docs = self.collection.find(filter=query,
                                    projection={"kegg_orthology_id": 1})
        levels = [2, 2157, 2759, 10239]
        bulk = []
        count = self.collection.count_documents(query)
        for i, doc in enumerate(docs):
            ids = []
            names = []
            if i == self.max_entries:
                break
            if self.verbose and i % 10 == 0:
                print("Processing doc {} out of {}...".format(i, count))
            ko = doc["kegg_orthology_id"]
            for l in levels:
                j = self.name_level(ko, level=l)
                if j == "":
                    ids.append(""),
                    names.append("")
                else:
                    n = self.orthodb.find_one({"orthodb_id": j})
                    ids.append(j)
                    names.append(n)
            bulk.append(UpdateOne({"_id": doc["_id"]},
                                  {"$set": {"orthodb_id": ids,
                                            "orthodb_name": names}}))
            if len(bulk) >= batch_size:
                print("     Bulk writing ...")
                self.collection.bulk_write(bulk)
                bulk = []
        if len(bulk) != 0:
            self.collection.bulk_write(bulk)
        print("Done.")

    def ko_passthrough(self,
                          batch_size=100):
        """Use KEGG and phylogenetic information to find orthodb group ID.

        Args:
            batch_size(:obj:`int`): Size of each bulk write
        """
        con_0 = {"orthodb_id": ""}
        con_1 = {"orthodb_id": None}
        comb = {"$or": [con_0, con_1]}
        comb_nor = {"$nor": [con_0, con_1]}
        con_2 = {"ko_number": {"$ne": None}}
        query = {"$and": [comb, con_2]}
        kos = self.collection.distinct("ko_number",
                                        filter=query)
        ancs = [2, 2157, 2759, 10239]
        for i, ko in enumerate(kos):
            if i == self.max_entries:
                break
            if self.verbose and i % 50 == 0:
                print("Processing ko {}, {} out of {}...".format(ko, i, len(kos)))
            for anc in ancs:
                if self.verbose and i % 50 == 0:
                    print("     Processing species with ancestor {}...".format(anc))
                con_3 = {"canon_anc_ids": anc}
                con_4 = {"ko_number": ko}
                q = {"$and": [comb_nor, con_3, con_4]}
                doc = self.collection.find_one(filter=q,
                                               projection={"orthodb_id": 1,
                                                           "orthodb_name": 1,
                                                           "gene_name": 1,
                                                           "protein_name": 1})
                if doc is None:
                    continue
                else:
                    orthodb_name = doc["orthodb_name"]
                    orthodb_id = doc["orthodb_id"]
                q = {"$and": [comb, con_3, con_4]}
                self.collection.update_many(q,
                                            {"$set": {"orthodb_id": orthodb_id,
                                                      "orthodb_name": orthodb_name}})

    def final_passthrough(self,
                          skip=0,
                          batch_size=100):
        """Use gene_name or protein_name to identify the remaining protein
        records without OrthoDB group ID, needs orthodb API.

        Args:
            batch_size(:obj:`int`): Bulk write size.
        """
        con_0 = {"ko_number": None}
        con_1 = {"orthodb_id": {"$in": [None, ""]}}
        query = {"$and": [con_0, con_1]}
        docs = self.collection.find(filter=query,
                                    projection={"protein_name": 1,
                                                "gene_name": 1,
                                                "canon_anc_ids": 1,
                                                "uniprot_id": 1},
                                    skip=skip)
        count = self.collection.count_documents(query)        
        bulk = []
        with open("./docs/orthodb/id_cache.json") as f:
            # try:
            cache = json.load(f)
            # except json.decoder.JSONDecodeError:
            #     cache = {}
        for i, doc in enumerate(docs):
            if i == self.max_entries:
                break
            if self.verbose and i % 10 == 0:
                print("Processing doc {} out of {}...".format(i, count))
                with open("./docs/orthodb/id_cache.json", "w+") as x:
                    json.dump(cache, x, indent=4)
            if len(doc.get("canon_anc_ids", [])) <= 1:
                continue
            dic_name = doc["protein_name"] if doc["gene_name"] is None else doc["gene_name"]
            tmp_0 = doc["canon_anc_ids"][1]
            if tmp_0 not in [2, 2157, 2759]:
                anc = 10239
            else:
                anc = tmp_0
            tmp_1 = cache.get("{}at{}_id".format(dic_name, anc))
            if tmp_1 is None:
                orthodb_id = self.name_level(dic_name, level=anc)
                cache["{}at{}_id".format(dic_name, anc)] = orthodb_id
                tmp_2 = self.collection.find_one({"orthodb_id": orthodb_id})
                orthodb_name = "" if tmp_2 is None else tmp_2["orthodb_name"]
                cache["{}at{}_name".format(dic_name, anc)] = orthodb_name
            else:
                orthodb_id = tmp_1
                orthodb_name = cache["{}at{}_name".format(dic_name, anc)]

            bulk.append(UpdateOne({"_id": doc["_id"]},
                                    {"$set": {"orthodb_id": orthodb_id,
                                            "orthodb_name": orthodb_name}}))
            if len(bulk) >= batch_size:
                print("     Bulk writing ...")
                self.collection.bulk_write(bulk)
                bulk = []
        if len(bulk) != 0:
            self.collection.bulk_write(bulk)
        print("Done!")



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
    # src.add_ortho(skip=0)

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
    #                       skip=30517000)

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

    # # add uniprot_id to orthodb_gene collection and add orthodb_gene info to uniprot collection
    # db = "datanator"
    # des_col = "orthodb_gene"
    # src = AddOrtho(MongoDB=conf.SERVER,
    #                 db=db,
    #                 des_col=des_col,
    #                 username=conf.USERNAME,
    #                 password=conf.PASSWORD,
    #                 verbose=True)
    # src.add_uniprot('./docs/orthodb/odb10v1_gene_xrefs.tab',
    #                 skip=64249000)

    db = "datanator"
    des_col = "rna_halflife_new"
    src = AddOrtho(MongoDB=conf.SERVER,
                    db=db,
                    des_col=des_col,
                    username=conf.USERNAME,
                    password=conf.PASSWORD,
                    verbose=True)
    src.add_x_ref_rna_halflife()

    # db = "datanator"
    # des_col = "rna_modification"
    # src = AddOrtho(MongoDB=conf.SERVER,
    #                 db=db,
    #                 des_col=des_col,
    #                 username=conf.USERNAME,
    #                 password=conf.PASSWORD,
    #                 verbose=True)
    # src.add_x_ref_rna_mod()

    # db = "datanator-test"
    # des_col = "uniprot"
    # src = AddOrtho(MongoDB=conf.SERVER,
    #                 db=db,
    #                 des_col=des_col,
    #                 username=conf.USERNAME,
    #                 password=conf.PASSWORD,
    #                 verbose=True)
    # src.ko_passthrough()

    # db = "datanator-test"
    # des_col = "uniprot"
    # src = AddOrtho(MongoDB=conf.SERVER,
    #                 db=db,
    #                 des_col=des_col,
    #                 username=conf.USERNAME,
    #                 password=conf.PASSWORD,
    #                 verbose=True)
    # src.final_passthrough()    


if __name__ == "__main__":
    main()
