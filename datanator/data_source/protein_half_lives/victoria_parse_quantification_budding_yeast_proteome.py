import pandas as pd
from pymongo import MongoClient
from datanator_query_python.util import mongo_util
from datanator_query_python.config import config

class parse_quantification_budding_yeast_proteome(mongo_util.MongoUtil):
    def __init__(self,
                 MongoDB=None,
                 db=None,
                 collection=None,
                 username=None,
                 password=None,
                 authSource = 'admin'):
        super(parse_quantification_budding_yeast_proteome,self).__init__(MongoDB=MongoDB, db=db,
                                                                         username = username,
                                                                         password = password,
                                                                         authSource = authSource)
        self.collection = collection
        self.entity_col = self.db_obj["entity"]
        self.identifier_col = self.db_obj["identifier"]
        self.obs_col = self.db_obj["observation"]
        
    def build_entity(self, data, i, doc):
        """Build entity object from obj.
        Go into entity collection
        Args:
            data (:obj:`Obj`): source object.
            i (:obj: `int`): index (row labels) of object in dataframe.
            doc (:obj: `Obj`): MongoDB document containing identifiers
        Return:
            (:obj:`Obj`), e.g.
            {
                "entity": {
                    "type": "protein",
                    "name": "Cytoplasmic protein",
                    "identifiers": [{}... {}]
                }
            }        
        """
        entity = {}
        entity["type"] = "protein"
        entity["name"] = doc["protein_name"]
        entity["identifiers"] = []
        entity["identifiers"].append({"namespace": "uniprot_id",
                                      "value": doc["uniprot_id"]})
        entity["identifiers"].append({"namespace": "gene_name_oln",
                                      "value": data.iloc[i,0]})
        if data.iloc[i,1]!=None:
            entity["identifiers"].append({"namespace": "gene_name",
                                          "value": data.iloc[i,1]})
        return entity

    def build_obs(self, data, i, doc):
        """Build observation objects from obj.
        Go into observations collection.
        Args:
            data (:obj:`Obj`): source object.
            i (:obj: `int`): index (row labels) of object in dataframe
            doc (:obj: `Obj`): MongoDB document containing identifier information.
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
        entity = {}
        entity["type"] = "protein"
        entity["name"] = doc["protein_name"]
        entity["identifiers"] = []
        entity["identifiers"].append({"namespace": "uniprot_id",
                                      "value": doc["uniprot_id"]})
        values_p = []
        if str(data.iloc[i,2])=="300":
            values_p.append({"type": "Half-life",
                             "value": str(float(data.iloc[i,2])*60),
                             "description": "non-exponential profile",
                             "units": "s"})
        else:
            values_p.append({"type": "Half-life",
                             "value": str(float(data.iloc[i,2])*60),
                             "units": "s"})

        genotype = {"taxon": {"ncbi_taxonomy_id":4932,
                              "name": "Saccharomyces cerevisiae"}}
        genotype["taxon"]["canon_ancestors"] = []
        query = {"tax_id":4932}
        projection = {"_id":0,"canon_anc_ids":1,"canon_anc_names":1}
        doc = self.client["datanator-test"]["taxon_tree"].find_one(filter = query, projection = projection)
        if doc!=None:
            for j in range(len(doc["canon_anc_names"])):
                d = {}
                d["ncbi_taxonomy_id"] = doc["canon_anc_ids"][j]
                d["name"] = doc["canon_anc_names"][j]
                genotype["taxon"]["canon_ancestors"].append(d)
        genotype["growthPhase"] = "log phase"
        environment = {}
        environment["media"] = "1.7 ml yeast extract/peptone and dextrose medium, 35 μg/ml cycloheximide"
        source = [{"namespace":"doi","value":"10.1073/pnas.0605420103"}]
        ob_p = {"entity":entity,
                "genotype":genotype,
                "environment":environment,
                "values":values_p,
                "source":source,
                "schema_version":"2.0"}
        return ob_p
        
    def process_docs(self):
        data = pd.read_csv('SuppDataSet.txt',delimiter='\t')
        data = data.where(pd.notnull(data), None)
        for i in range(len(data)):
            query = {"gene_name_oln":data.iloc[i,0]}
            projection = {"_id":0,"uniprot_id":1,"protein_name":1}
            doc = self.client["datanator"]["uniprot"].find_one(filter=query,projection=projection)
            if doc!=None:
                #update entity collection
                entity = self.build_entity(data, i, doc)
                query = {"identifiers":{"$elemMatch":{"namespace":"uniprot_id",
                                                      "value": doc["uniprot_id"]}}}
                self.entity_col.update_one(query,
                                           {"$set": {
                                                 "type": entity["type"],
                                                 "name": entity["name"],                      
                                                 "schema_version": "2.0"},
                                            "$addToSet": {
                                                "identifiers": {"$each": entity["identifiers"]}}},
                                           upsert=True)
                obs = self.build_obs(data, i, doc)
                #update identifier collection
                query = {"$and":[{"namespace":"uniprot_id"},
                                 {"value":doc["uniprot_id"]}]}
                self.identifier_col.update_one(query,
                                            {"$set": {"namespace": "uniprot_id",
                                                      "value": doc["uniprot_id"],
                                                      "description": entity["name"]}},upsert=True)
                #update observation collection
                con_1 = {"source":{"$elemMatch":{"namespace":"doi","value":"10.1073/pnas.0605420103"}}}
                con_2 = {"identifier":{"namespace":"uniprot_id","value":doc["uniprot_id"]}}
                query = {"$and":[con_1,con_2]}
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

    src = parse_quantification_budding_yeast_proteome(MongoDB = MongoDB,
                                                      username=username,
                                                      password=password,
                                                      collection = "observation",
                                                      db = "datanator-demo")
    src.process_docs()
    
if __name__== '__main__':
    main()
