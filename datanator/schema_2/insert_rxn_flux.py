from pathlib import Path
from datanator_query_python.util import mongo_util
import pandas as pd
from bioservices import *
import copy
from pymongo.collation import Collation, CollationStrength
import math


class InsRxn(mongo_util.MongoUtil):
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
        self.kegg = KEGG()
        self.taxon = self.client["datanator-test"]["taxon_tree"]
        self.collation = Collation(locale='en', strength=CollationStrength.SECONDARY)
        self.NON_KEGG = ["unknown", "transport", "unknow"]

    def parse_docs(self,
                   _dir="./docs/"):
        """Parse database with xls files from
        http://www.cecafdb.org/, stored in ./docs/

        Args:
            _dir(:obj:`str`): Directory in which xls files are stored.
        """
        kegg = {}
        obj = {
            "entity": {
                "type": "reaction",
                "name": "",
                "identifiers": [],
                "schema_version": "2.0"
            },
            "source": [],
            "genotype": {
                "ncbi_taxonomy_id": "",
                "name": "",
                "canon_ancestors": []
            },
            "schema_version": "2.0",
            "identifier": {},
            "values": {"reaction_flux": []},
            "environment": {}
        }
        paths = Path(_dir).glob('**/*.xls')
        for i, path in enumerate(paths):
            if self.verbose:
                print(path)
            carbons = set()
            case = 0
            case_obj = {}  # {"case_0": {}, "case_1": {}, ...} common elements of a case object
            final_obj = {} # {"case_0": [{rxn_obs_0}, {rxn_obs_1}, {rxn_obs_2}, ...], "case_1": [], ...}
            if i < 2:
                continue
            if i == self.max_entries:
                break
            df = pd.read_excel(str(path),
                              sheet_name="sheet one",
                              skiprows=[0,1,2,3,7,15,16,
                                        17,18,19,20,21])
            if self.verbose:
                print(df.shape[0], df.shape[1])
            for i, row in df.iterrows():
                if not isinstance(row[0], str) and math.isnan(row[0]):
                    continue
                elif row[0] == "Experiment Name (ID)":
                    obj["source"].append({"namespace": "publication",
                                           "value": row[1]})
                elif row[0] == "Coordinator":
                    obj["source"].append({"namespace": "publication_url",
                                           "value": row[1]})
                elif row[0] == "Strains":
                    for species in row[1:]:
                        if not isinstance(species, str) and math.isnan(species):
                            break
                        else:
                            case_obj["case_{}".format(case)] = copy.deepcopy(obj)                            
                            case_obj["case_{}".format(case)]["genotype"] = self.build_genotype_obj(species)
                            final_obj["case_{}".format(case)] = []
                            case += 1  # number of columns (cases)
                elif row[0] == "culture medium":
                    for j, (k, o) in enumerate(case_obj.items()):
                        o["environment"]["media"] = row[j+1]
                elif row[0] == "carbon source":
                    for j, (k, o) in enumerate(case_obj.items()):
                        o["environment"]["carbon_source"] = row[j+1]
                        if row[j+1].strip() == "D-Glucose":
                            carbons.add("glucose")
                        carbons.add(row[j+1].strip())
                elif row[0] == "growth rate":
                    for j, (k, o) in enumerate(case_obj.items()):
                        o["environment"]["growth_rate"] = row[j+1]
                elif row[0] == "Case-specific description":
                    for j, (k, o) in enumerate(case_obj.items()):
                        o["environment"]["condition"] = row[j+1]
                elif "<==>" in row[0]:
                    kegg_rxn_id = row[1]
                    if self.verbose:
                        print(kegg_rxn_id)
                    if kegg_rxn_id.strip() in self.NON_KEGG:  # non kegg description rather than rn:RXXXXXX
                        print("Condition 0: {}".format(kegg_rxn_id.strip()))
                        s, p = self.parse_rxn_str(row[0])
                        s_n = s
                        p_n = p
                    elif kegg.get(kegg_rxn_id) is None:
                        print("Condition 1: {}".format(kegg_rxn_id.strip()))
                        s, s_n, p, p_n = self.get_kegg_rxn(kegg_rxn_id)
                        tmp = {"substrates": s,
                               "products": p,
                               "s_n": s_n,
                               "p_n": p_n}
                        kegg[kegg_rxn_id] = tmp
                    else:
                        s = kegg[kegg_rxn_id]["substrates"]
                        p = kegg[kegg_rxn_id]["products"]
                        s_n = kegg[kegg_rxn_id]["s_n"]
                        p_n = kegg[kegg_rxn_id]["p_n"]
                    for j, (k, o) in enumerate(case_obj.items()):
                        tmp = copy.deepcopy(o)
                        tmp["entity"]["identifiers"].append({"namespace": "kegg_reaction_id",
                                                             "value": kegg_rxn_id})
                        tmp["identifier"] = {"namespace": "reaction_eqn_ik",
                                             "substrates": s,
                                             "products": p}
                        tmp["entity"]["identifiers"].append({"namespace": "reaction_eqn_name",
                                                             "substrates": s_n,
                                                             "products": p_n})
                        if row[j+3] == "X":
                            val = 0
                        else:
                            val = row[j+3]
                        tmp["values"]["reaction_flux"].append({"value": val,
                                                               "units": "mmol*(g^-1)*(h^-1)"})
                        final_obj[k].append(tmp)
                elif row[0].lower() in carbons:
                    for j, (_, obs) in enumerate(final_obj.items()):
                        rate = row[j+1]
                        for o in obs:
                            if "WQZGKKKJIJFFOK-DVKNGEFBSA-N" in o["identifier"]["substrates"]: #glucose
                                o["values"]["reaction_flux"][0]["value"] = rate
                            else:
                                o["values"]["reaction_flux"][0]["value"] *= rate / 100
                else:
                    continue
            for _, v in final_obj.items():
                self.collection.insert_many(v)

    def get_kegg_rxn(self, _id):
        """Use bioservice to request kegg reaction information.

        Args:
            _id(:obj:`str`): Kegg reaction id.

        Return:
            (:obj:`tuple` of :obj:`list`): substrates and products in inchikey lists
        """
        agg = self.kegg.get(_id)
        obj = self.kegg.parse(agg)
        eqn = obj["EQUATION"]
        defi = obj["DEFINITION"]
        substrates = []
        products = []
        for i, (c, d) in enumerate(zip(eqn.split(" <=> "), defi.split(" <=> "))):
            if i == 0:  # substrates string
                # kegg compound id to ChEBI id
                tmp = [KEGG().parse(KEGG().get(x))["DBLINKS"]["ChEBI"] for x in c.split(" + ")]
                # chebi id to inchikey
                for x in tmp:
                    try:
                        substrates.append(ChEBI().getCompleteEntity("CHEBI:"+x[0:5]).inchiKey)
                    except AttributeError:
                        substrates.append(ChEBI().getCompleteEntity("CHEBI:"+x[0:5]).smiles)
                substrates_name = [x for x in d.split(" + ")]
            elif i == 1: # products string
                tmp = [KEGG().parse(KEGG().get(x))["DBLINKS"]["ChEBI"] for x in c.split(" + ")]
                # chebi id to inchikey
                for x in tmp:
                    try:
                        products.append(ChEBI().getCompleteEntity("CHEBI:"+x[0:5]).inchiKey)
                    except AttributeError:
                        products.append(ChEBI().getCompleteEntity("CHEBI:"+x[0:5]).smiles)
                products_name = [x for x in d.split(" + ")]
        return substrates, substrates_name, products, products_name

    def parse_rxn_str(self, rxn):
        """Parse reaction string given in documents.

        Args:
            (:obj:`str`): Reaction equation.

        Return:
            (:obj:`tuple` of :obj:`lists`)
        """
        substrates = []
        products = []
        for i, c in enumerate(rxn.split(" <==> ")):
            if i == 0:
                substrates = [x for x in c.split(" + ")]
            elif i == 1:
                products = [x for x in c.split(" + ")]
        return substrates, products

    def build_genotype_obj(self, name):
        """Given species name, build genotype object.

        Args: 
            (:obj:`str`): Name of the species.

        Return:
            (:obj:`Obj`)
        """
        con_0 = {"tax_name": name}
        con_1 = {"name_txt": name}
        query = {"$or": [con_0, con_1]}
        obj = {"canon_ancestors": []}
        doc = self.taxon.find_one(filter=query,
                                  projection={"_id": 0,
                                              "tax_name": 1,
                                              "tax_id": 1,
                                              "canon_anc_ids": 1,
                                              "canon_anc_names": 1},
                                  collation=self.collation)
        if doc is not None:
            obj["ncbi_taxonomy_id"] = doc["tax_id"]
            obj["name"] = doc["tax_name"]
            for _id, name in zip(doc["canon_anc_ids"], doc["canon_anc_names"]):
                obj["canon_ancestors"].append({"ncbi_taxonomy_id": _id,
                                               "name": name})
            return obj
        else:
            name_sans_strain = name.split(" ")
            name = "{} {}".format(name_sans_strain[0], name_sans_strain[1])
            con_0 = {"tax_name": name}
            con_1 = {"name_txt": name}
            query = {"$or": [con_0, con_1]}
            doc = self.taxon.find_one(filter=query,
                                        projection={"_id": 0,
                                                    "tax_name": 1,
                                                    "tax_id": 1,
                                                    "canon_anc_ids": 1,
                                                    "canon_anc_names": 1},
                                        collation=self.collation)
            if doc is not None:
                obj["ncbi_taxonomy_id"] = doc["tax_id"]
                obj["name"] = doc["tax_name"]
                for _id, name in zip(doc["canon_anc_ids"], doc["canon_anc_names"]):
                    obj["canon_ancestors"].append({"ncbi_taxonomy_id": _id,
                                                "name": name})
                return obj
            else:
                return obj