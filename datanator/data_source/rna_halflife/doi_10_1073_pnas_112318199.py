from datanator.util import rna_halflife_util, x_ref
from datanator_query_python.config import config
from pymongo import UpdateOne
import chardet
from pymongo.collation import Collation, CollationStrength
import csv


class Halflife(rna_halflife_util.RnaHLUtil):

    def __init__(self, cache_dir=None, server=None, src_db=None, protein_col=None,
                authDB=None, readPreference=None, username=None, password=None,
                verbose=None, max_entries=None, des_db=None, rna_col=None):
        """Init
        
        Args:
            cache_dir (:obj:`str`, optional): Cache directory for logs. Defaults to None.
            server (:obj:`str`, optional): MongoDB server address. Defaults to None.
            db (:obj:`str`, optional): Database where initial uniprot collection resides. Defaults to None.
            collection_str (:obj:`str`, optional): name of collection. Defaults to None.
            authDB (:obj:`str`, optional): MongoDB authentication database. Defaults to None.
            readPreference (:obj:`str`, optional): MongoDB read preference. Defaults to None.
            username (:obj:`str`, optional): MongoDB username. Defaults to None.
            password (:obj:`str`, optional): MongoDB password. Defaults to None.
            verbose (:obj:`bool`, optional): Wheter to display verbose messages. Defaults to None.
            max_entries (:obj:`int`, optional): Number of records to be processed. Defaults to None.
            uniprot_col_db (:obj:`int`, optional): Database to which new uniprot records will be inserted. Defaults to None.
        """
        super().__init__(server=server, username=username, password=password, src_db=src_db,
        des_db=des_db, protein_col=protein_col, rna_col=rna_col, authDB=authDB, readPreference=readPreference,
        max_entries=max_entries, verbose=verbose, cache_dir=cache_dir)
        self.xref = x_ref.XRef(MongoDB=server,
                               db=src_db,
                               username=username,
                               password=password)
        self.uniprot_col = self.client["datanator-test"]["uniprot"]
        self.collation = Collation('en', strength=CollationStrength.SECONDARY)
        self.max_entries = max_entries
        self.verbose = verbose

    def fill_rna_halflife(self, 
                url,
                batch_size=100,
                skip=2):
        """Fill rna_halflife collection

        Args:
            (:obj:`str`): URL of the file.
            (:obj:`int`, optional): Number of docs to be inserted at once.
            (:obj:`int`, optional): First number of rows to skip.
        """
        with open(url, newline='') as f:
            # result = chardet.detect(f.read())
            # dialect = csv.Sniffer().sniff(f.read(1024))
            # f.seek(0)
            x = csv.reader(f,
                           delimiter=",")
            tax_id = 511145   # parent 83333
            tax_name = "Escherichia coli str. K-12 substr. MG1655"
            batch = []
            for i, row in enumerate(x):
                halflives = []
                if i == self.max_entries:
                    break
                if i < skip:
                    continue
                if self.verbose and i % 50 == 0:
                    print("Process row {} with gene symbol {} ...".format(i, gene_symbol))
                oln = row[0].lower()
                gene_symbol = row[1]
                gene_name = row[2]
                con_0 = {"add_id.value": {"$in": [gene_symbol, oln]}}
                con_1 = {"gene_name": gene_symbol}
                con_2 = {"ncbi_taxonomy_id": {"$in": [tax_id, 83333]}}
                query = {"$and": [{"$or": [con_0, con_1]}, con_2]}
                uniprot_obj = self.uniprot_col.find_one(filter=query)
                if uniprot_obj is not None:
                    uniprot_id = uniprot_obj["uniprot_id"]
                    orthodb_id = uniprot_obj["orthodb_id"]
                    orthodb_name = uniprot_obj["orthodb_name"]
                else:
                    print("     Row {} contains entry not found in UniProt collection ...".format(i + 1))
                    x = self.xref.gene_tax_to_uniprot(gene_symbol, 83333)        
                    if x == {}:
                        x = self.xref.gene_tax_to_uniprot(oln, 83333)
                        if x == {}: 
                            print("         Row {} contains entry not found in UniProt website ...".format(i + 1))
                            continue
                    uniprot_id = x["uniprot_id"]
                    orthodb_id = x["orthodb_id"]
                    orthodb_name = x["orthodb_name"]
                halflife_lb = float(row[3]) * 60 if row[3] != "" else None # unit in seconds
                halflife_m9 = float(row[4]) * 60 if row[4] != "" else None
                if halflife_lb is not None:
                    obj_0 = {"gene_name": gene_symbol,
                             "gene_description": gene_name,
                             "species": tax_name,
                             "ncbi_taxonomy_id": tax_id,
                             "ordered_locus_name": oln,
                             "halflife": halflife_lb,
                             "unit": "s",
                             "reference": {"doi": "10.1073/pnas.112318199"},
                             "growth_medium": "LB + 0.2% glucose medium at 30ºC"}
                    halflives.append(obj_0)
                if halflife_m9 is not None:
                    obj_1 = {"gene_name": gene_symbol,
                             "gene_description": gene_name,
                             "species": tax_name,
                             "ordered_locus_name": oln,
                             "ncbi_taxonomy_id": tax_id,
                             "halflife": halflife_m9,
                             "unit": "s",
                             "reference": {"doi": "10.1073/pnas.112318199"},
                             "growth_medium": "M9 + 0.2% glucose medium at 30ºC"}
                    halflives.append(obj_1)
                if halflives == []:
                    continue
                else:
                    batch.append(UpdateOne({"uniprot_id": uniprot_id},
                                           {"$set": {"orthodb_id": orthodb_id,
                                                     "orthodb_name": orthodb_name},
                                            "$addToSet": {"halflives": {"$each": halflives},
                                                          "protein_names": gene_name}}, upsert=True))
                
                if len(batch) >= batch_size:
                    print("     Bulk updating...")
                    self.rna_hl_collection.bulk_write(batch)
                    batch = []
            if len(batch) != 0:
                print("     Bulk updating...")
                self.rna_hl_collection.bulk_write(batch)
            print("Done.")  


def main():
    conf = config.DatanatorAdmin()         

    # add to RNA halflife
    des_db = 'datanator'
    src_db = 'datanator-test'
    protein_col = 'uniprot'
    rna_col = 'rna_halflife_new'
    username = conf.USERNAME
    password = conf.PASSWORD
    MongoDB = conf.SERVER
    src = Halflife(server=MongoDB, src_db=src_db,
        protein_col=protein_col, authDB='admin', readPreference='nearest',
        username=username, password=password, verbose=True,
        des_db=des_db, rna_col=rna_col)
    src.fill_rna_halflife("./docs/rna_halflife/3181Table5.csv", skip=9)


if __name__ == "__main__":
    main()