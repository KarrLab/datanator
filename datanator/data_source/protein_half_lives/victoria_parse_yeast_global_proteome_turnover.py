import pandas as pd
import bson
from pymongo import MongoClient
from datanator_query_python.util import mongo_util
from datanator_query_python.config import config

class parse_yeast_global_proteome_turnover(mongo_util.MongoUtil):
    def __init__(self,
                 MongoDB=None,
                 db=None,
                 collection=None,
                 username=None,
                 password=None,
                 authSource = 'admin',
                 readPreference = 'nearest'):
        super(parse_yeast_global_proteome_turnover,self).__init__(MongoDB=MongoDB, db=db,
                                                                  username = username,
                                                                  password = password,
                                                                  authSource = authSource,
                                                                  readPreference=readPreference)
        self.collection = collection
        self.entity_col = self.db_obj["entity"]
        self.identifier_col = self.db_obj["identifier"]
        self.obs_col = self.db_obj["observation"]
        
    def build_entity(self, data, i):
        """Build entity object from obj.
        Go into entity collection

        Args:
            data (:obj:`Obj`): source object.
            i (:obj: `int`): index (row labels) of object in dataframe.
        Return:
            (:obj:`Obj`), e.g.
            {
                "entity": {
                    "type": "protein",
                    "name": "Type 2A phosphatase-associated protein 42",
                    "identifiers": [{"namespace": "uniprot_id",
                                     "value": "Q04372"},
                                    {"namespace": "ensembl",
                                     "value": "YMR028W"},
                                    {"namespace": "gene_name",
                                     "value": "TAP42"}
                                   ],
                    "synonyms": []
                }
            }        
        """
        entity = {}
        entity["type"] = "protein"
        names = data.iloc[i,8].split(";")
        entity["name"] = names[0]
        entity["identifiers"] = []
        entity["identifiers"].append({"namespace": "uniprot_id",
                                      "value": data.iloc[i,0]})
        if data.iloc[i,1]!=None:
            ensembl_ids = data.iloc[i,1].split(";")
            for ensembl in ensembl_ids:
                entity["identifiers"].append({"namespace": "ensembl",
                                              "value": ensembl})
        if data.iloc[i,7]!=None:
            gene_names = data.iloc[i,7].split(";")
            for gene in gene_names:
                entity["identifiers"].append({"namespace": "gene_name",
                                              "value": gene})   
        entity["synonyms"] = []
        if len(names)>1:
            entity["synonyms"] = names[1:]
        return entity

    def build_obs(self, data, i, species):
        """Build observation objects from obj.
        Go into observations collection.
        Args:
            data (:obj:`Obj`): source object.
            i (:obj: `int`): index (row labels) of object in dataframe.
            species (:obj:`str`): species of yeast.
        Return:
            obj(:obj:`Obj`)
            {
                "entity": {
                    "type": "protein",
                    "name": "Type 2A phosphatase-associated protein 42",
                    "identifiers": [{"namespace": "uniprot_id",
                                     "value": "Q04372"}]
                },
                "genotype":{
                    "taxon": {}
                },
                "environment:{
                    "media":
                },
                "values": [],
                "source": {}, ...
            }
        """
        entity = {}
        entity["type"] = "protein"
        names = data.iloc[i,8].split(";")
        entity["name"] = names[0]
        entity["identifiers"] = []
        entity["identifiers"].append({"namespace": "uniprot_id",
                                      "value": data.iloc[i,0]})
        values_p = []
        if data.iloc[i,4]!="n.d.":
            values_p.append({"type": "Half-life",
                             "value": str(float(data.iloc[i,4])*60),
                             "units": "s"})
        else:
            if ">=" in data.iloc[i,5]:
                values_p.append({"type": "Half-life",
                                 "value": "greater than or equal to "+str(float(data.iloc[i,5][3:])*3600),
                                 "units": "s"})
        if data.iloc[i,2]!="n.d.":
            values_p.append({"type": "Degradation rates",
                             "value": str(float(data.iloc[i,2])/60),
                             "units": "s^(-1)"})
        values_p.append({"type": "R^2 (quality of curve fitting)",
                         "value": data.iloc[i,3]})
        values_p.append({"type": "Cross validation of slope",
                         "value": data.iloc[i,6]})
        environment = {}
        query={}
        if species=="Saccharomyces cerevisiae BY4742":
            genotype = {"taxon": {"ncbi_taxonomy_id": 559292,
                                  "name": species}}
            query = {"tax_id": 559292}
            environment["media"] = "5 ml synthetic medium, 30 mg/l heavy [13C6/15N2] L-lysine, 6.7 g/l yeast nitrogen base, 2 g/l dropout mix, all amino acids except lysine and 2% glucose"
            environment["temperature"] = bson.Int64(30)
        elif species=="Schizosaccharomyces pombe MKSP201":
            genotype = {"taxon": {"ncbi_taxonomy_id": 4896,
                                  "name": species}}
            query = {"tax_id": 4896}
            environment["media"] = "Edinburgh minimal medium supplemented with 75 mg/l leucine, histidine, uracil, and adenine, heavy [13C6/15N2] L-lysine"
        genotype["taxon"]["canon_ancestors"] = []
        projection = {"_id":0,"canon_anc_ids":1,"canon_anc_names":1}
        doc = self.client["datanator-test"]["taxon_tree"].find_one(filter = query, projection = projection)
        if doc!=None:
            for j in range(len(doc["canon_anc_names"])):
                d = {}
                d["ncbi_taxonomy_id"] = doc["canon_anc_ids"][j]
                d["name"] = doc["canon_anc_names"][j]
                genotype["taxon"]["canon_ancestors"].append(d)
        source = [{"namespace":"doi","value":"10.1016/j.celrep.2014.10.065"}]
        ob_p = {"entity":entity,
                "genotype":genotype,
                "environment":environment,
                "values":values_p,
                "source":source,
                "schema_version":"2.0"}
        return ob_p
    
    def build_obs_multi_ids(self, data, i, species):
        """Build observation object from data with more than one uniprot_id per row
        Go into observations collection

        Args:
            data (:obj:`Obj`): source object.
            i (:obj: `int`): index (row labels) of object in dataframe.
            species (:obj:`str`): species of yeast.
        Return:
            obj(:obj:`Obj`)
            {
                "entity": {
                    "type": "protein",
                    "name": "Cytoplasmic protein",
                    "identifiers": [{}... {}]
                },
                "genotype":{
                    "taxon": {}
                },
                "values": [],
                "source": {}, ...
            } 
        """
        uniprot_ids = data.iloc[i,0].split(";")
        for uniprot in uniprot_ids:
            query = {"identifiers.value":uniprot}
            projection = {"_id":0, "name":1}
            doc = self.client["datanator-demo"]["entity"].find_one(filter=query,projection=projection)
            if doc!=None:
                entity = {}
                entity["type"] = "protein"
                entity["name"] = doc["name"]
                entity["identifiers"] = []
                entity["identifiers"].append({"namespace": "uniprot_id",
                                              "value": uniprot})
                values_p = []
                if data.iloc[i,4]!="n.d.":
                    values_p.append({"type": "Half-life",
                                     "value": str(float(data.iloc[i,4])*60),
                                     "units": "s"})
                else:
                    if ">=" in data.iloc[i,5]:
                        values_p.append({"type": "Half-life",
                                         "value": "greater than or equal to "+str(float(data.iloc[i,5][3:])*3600),
                                         "units": "s"})
                if data.iloc[i,2]!="n.d.":
                    values_p.append({"type": "Degradation rates",
                                     "value": str(float(data.iloc[i,2])/60),
                                     "units": "s^(-1)"})
                values_p.append({"type": "R^2 (quality of curve fitting)",
                                 "value": data.iloc[i,3]})
                values_p.append({"type": "Cross validation of slope",
                                 "value": data.iloc[i,6]})
                environment = {}
                query={}
                if species=="Saccharomyces cerevisiae BY4742":
                    genotype = {"taxon": {"ncbi_taxonomy_id": 559292,
                                          "name": species}}
                    query = {"tax_id": 559292}
                    environment["media"] = "5 ml synthetic medium, 30 mg/l heavy [13C6/15N2] L-lysine, 6.7 g/l yeast nitrogen base, 2 g/l dropout mix, all amino acids except lysine and 2% glucose"
                    environment["temperature"] = bson.Int64(30)
                elif species=="Schizosaccharomyces pombe MKSP201":
                    genotype = {"taxon": {"ncbi_taxonomy_id": 4896,
                                          "name": species}}
                    query = {"tax_id": 4896}
                    environment["media"] = "Edinburgh minimal medium supplemented with 75 mg/l leucine, histidine, uracil, and adenine, heavy [13C6/15N2] L-lysine"
                genotype["taxon"]["canon_ancestors"] = []
                projection = {"_id":0,"canon_anc_ids":1,"canon_anc_names":1}
                taxon_doc = self.client["datanator-test"]["taxon_tree"].find_one(filter = query, projection = projection)
                if taxon_doc!=None:
                    for j in range(len(taxon_doc["canon_anc_names"])):
                        d = {}
                        d["ncbi_taxonomy_id"] = taxon_doc["canon_anc_ids"][j]
                        d["name"] = taxon_doc["canon_anc_names"][j]
                        genotype["taxon"]["canon_ancestors"].append(d)
                source = [{"namespace":"doi","value":"10.1016/j.celrep.2014.10.065"}]
                ob_p = {"entity":entity,
                        "genotype":genotype,
                        "environment":environment,
                        "values":values_p,
                        "source":source,
                        "schema_version":"2.0"}
                query = {"$and":[{"namespace":"uniprot_id"},
                                     {"value":entity["identifiers"][0]["value"]}]}
                self.identifier_col.update_one(query,
                                               {"$set": {"namespace": "uniprot_id",
                                                         "value": entity["identifiers"][0]["value"]}},
                                               upsert=True)
                #update observation collection
                con_1 = {"source":{"$elemMatch":{"namespace":"doi","value":"10.1016/j.celrep.2014.10.065"}}}
                con_2 = {"identifier":{"namespace":"uniprot_id","value": uniprot}}
                query = {"$and": [con_1,con_2]}                                    
                self.obs_col.update_one(query,
                                    {"$set": {"entity": ob_p["entity"],
                                              "genotype": ob_p["genotype"],
                                              "environment": ob_p["environment"],
                                              "schema_version": "2.0",
                                              "identifier":{"namespace":"uniprot_id","value": uniprot}},
                                     "$addToSet": {"values": {"$each": ob_p["values"]},
                                                   "source": {"$each": ob_p["source"]}}},
                                     upsert=True)
    
    def process_docs(self):
        files = ['mmc2.txt','mmc3.txt']
        for file in files:
            table = pd.read_csv(file,
                                delimiter='\t',
                                dtype={"UNIPROT":str,
                                       "ENSG":str,
                                       "Degradation rates (min-1)":str,
                                       "R2 (quality of curve fitting)":str,
                                       "t1/2 (min)":str,
                                       "t1/2 (hours)":str,
                                       "Cv of slope":str,
                                       "Gene name":str,
                                       "Protein name":str})
            table = table.where(pd.notnull(table), None)
            for i in range(len(table)):
                uniprot_ids = table.iloc[i,0].split(";")
                names = table.iloc[i,8]
                species = "Saccharomyces cerevisiae BY4742"
                if file == "mmc3.txt":
                    species = "Schizosaccharomyces pombe MKSP201"
                if len(uniprot_ids)>1:
                    self.build_obs_multi_ids(table,i,species)
                elif len(uniprot_ids)==1 and names!=None:
                    #update entity collection
                    entity = self.build_entity(table,i)
                    query = {"identifiers":{"$elemMatch":{"namespace":"uniprot_id",
                                                          "value":entity["identifiers"][0]["value"]}}}
                    self.entity_col.update_one(query,
                                               {"$set": {
                                                     "type": entity["type"],
                                                     "name": entity["name"],                      
                                                     "schema_version": "2.0"},
                                                "$addToSet": {
                                                    "identifiers": {"$each": entity["identifiers"]},
                                                    "synonyms": {"$each": entity["synonyms"]}}},
                                               upsert=True)
                    obs = self.build_obs(table,i,species)

                    #update identifier collection
                    query = {"$and":[{"namespace":"uniprot_id"},
                                     {"value":entity["identifiers"][0]["value"]}]}
                    self.identifier_col.update_one(query,
                                                   {"$set": {"namespace": "uniprot_id",
                                                             "value": entity["identifiers"][0]["value"]}},
                                                   upsert=True)

                    #update observation collection
                    con_1 = {"source":{"$elemMatch":{"namespace":"doi","value":"10.1016/j.celrep.2014.10.065"}}}
                    con_2 = {"identifier":{"namespace":"uniprot_id","value":entity["identifiers"][0]["value"]}}
                    query = {"$and": [con_1,con_2]}                                    
                    self.obs_col.update_one(query,
                                        {"$set": {"entity": obs["entity"],
                                                  "genotype": obs["genotype"],
                                                  "environment": obs["environment"],
                                                  "schema_version": "2.0",
                                                  "identifier":{"namespace":"uniprot_id","value": entity["identifiers"][0]["value"]}},
                                         "$addToSet": {"values": {"$each": obs["values"]},
                                                       "source": {"$each": obs["source"]}}},
                                         upsert=True)
def main():
    conf=config.Victoria()
    conf_main = config.Config()
    username = conf.USERNAME
    password = conf.PASSWORD
    MongoDB = conf_main.SERVER

    src = parse_yeast_global_proteome_turnover(MongoDB = MongoDB,
                                               username=username,
                                               password=password,
                                               collection = "observation",
                                               db = "datanator-demo")
    
    src.process_docs()

if __name__== '__main__':
    main()
