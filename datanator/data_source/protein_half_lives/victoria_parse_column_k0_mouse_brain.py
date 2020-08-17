import pandas as pd
from pymongo import MongoClient
from datanator_query_python.util import mongo_util
from datanator_query_python.config import config

class parse_column_k0_mouse_brain(mongo_util.MongoUtil):
    def __init__(self,
                 MongoDB=None,
                 db=None,
                 collection=None,
                 username=None,
                 password=None,
                 authSource = 'admin',
                 readPreference = 'nearest'):
        super(parse_column_k0_mouse_brain,self).__init__(MongoDB=MongoDB, db=db,
                                                         username = username,
                                                         password = password,
                                                         authSource = authSource,
                                                         readPreference=readPreference)
        self.collection = collection
        self.entity_col = self.db_obj["entity"]
        self.identifier_col = self.db_obj["identifier"]
        self.obs_col = self.db_obj["observation"]
      
    def build_entity_peptide(self, data, i):
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
                    "name": "Cytoplasmic protein",
                    "identifiers": [{}... {}]
                }
            }        
        """
        entity = {}
        entity["type"] = "protein"    
        entity["name"] = data.iloc[i,0]
        entity["identifiers"] = []
        entity["identifiers"].append({"namespace": "uniprot_id",
                                      "value": data.iloc[i,1]})
        return entity
    
    def build_entity_protein(self, data, i):
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
                    "name": "Cytoplasmic protein",
                    "identifiers": [{}... {}]
                }
            }        
        """
        entity = {}
        entity["type"] = "protein"
        entity["identifiers"] = []
        entity["identifiers"].append({"namespace": "uniprot_id",
                                      "value": data.iloc[i,0]})
        query = {"$and":[{"namespace":"uniprot_id"},
                         {"value":entity["identifiers"][0]["value"]}]}
        projection = {"_id":0,"description":1}
        doc = self.client["datanator-demo"]["identifier"].find_one(filter = query, projection = projection)
        if doc!=None:
            entity["name"] = doc["description"]
        return entity
    
    def build_obs_peptide(self, data, i, description):
        """Build observation objects from obj.
        Go into observations collection.
        Args:
            data (:obj:`Obj`): source object.
            i (:obj: `int`): index (row labels) of object in dataframe
            description (:obj: `str`): description of proteins analyzed, e.g. "brain_peptide","liver_peptide","blood_peptide"
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
        entity["name"] = data.iloc[i,0]
        entity["identifiers"] = []
        entity["identifiers"].append({"namespace": "uniprot_id",
                                      "value": data.iloc[i,1]})

        values_p = []
        values_p.append({"type": "k0 (turnover rate)",
                         "value": data.iloc[i,2],
                         "description": description
                         })

        genotype = {"taxon": {"ncbi_taxonomy_id":10090,
                              "name":"Mus musculus"}}
        genotype["taxon"]["canon_ancestors"] = []
        query = {"tax_id":10090}
        projection = {"_id":0,"canon_anc_ids":1,"canon_anc_names":1}
        doc = self.client["datanator-test"]["taxon_tree"].find_one(filter = query, projection = projection)
        if doc!=None:
            for j in range(len(doc["canon_anc_names"])):
                d = {}
                d["ncbi_taxonomy_id"] = doc["canon_anc_ids"][j]
                d["name"] = doc["canon_anc_names"][j]
                genotype["taxon"]["canon_ancestors"].append(d)
                
        source = [{"namespace":"doi","value":"10.1073/pnas.1006551107"}]
        
        ob_p = {"entity":entity,
                "genotype":genotype,
                "values":values_p,
                "source":source,
                "schema_version":"2.0"}
        return ob_p

    def build_obs_protein(self, data, i, description):
        """Build observation objects from obj.
        Go into observations collection.
        Args:
            data (:obj:`Obj`): source object.
            i (:obj: `int`): index (row labels) of object in dataframe
            description (:obj: `str`): description of proteins analyzed, e.g. "brain_protein","liver_protein", "blood_protein"
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
        entity["identifiers"] = []
        entity["identifiers"].append({"namespace": "uniprot_id",
                                      "value": data.iloc[i,0]})
        query = {"$and":[{"namespace":"uniprot_id"},
                         {"value":entity["identifiers"][0]["value"]}]}
        projection = {"_id":0,"description":1}
        doc = self.client["datanator-demo"]["identifier"].find_one(filter = query, projection = projection)
        if doc!=None:
            entity["name"] = doc["description"]
            
        values_p = []
        
        values_p.append({"type": "k0 (turnover rate)",
                         "value": data.iloc[i,1],
                         "description": description
                         })
        genotype = {"taxon": {"ncbi_taxonomy_id":10090,
                              "name":"Mus musculus"}}
        genotype["taxon"]["canon_ancestors"] = []
        query = {"tax_id":10090}
        projection = {"_id":0,"canon_anc_ids":1,"canon_anc_names":1}
        doc = self.client["datanator-test"]["taxon_tree"].find_one(filter = query, projection = projection)
        if doc!=None:
            for j in range(len(doc["canon_anc_names"])):
                d = {}
                d["ncbi_taxonomy_id"] = doc["canon_anc_ids"][j]
                d["name"] = doc["canon_anc_names"][j]
                genotype["taxon"]["canon_ancestors"].append(d)
                
        source = [{"namespace":"doi","value":"10.1073/pnas.1006551107"}]
        
        ob_p = {"entity":entity,
                "genotype":genotype,
                "values":values_p,
                "source":source,
                "schema_version":"2.0"}
        return ob_p
    
    def process_docs_peptide(self):
        peptide_files = ['brain_peptide.txt','liver_peptide.txt','blood_peptide.txt']
        for file in peptide_files:
            data = pd.read_csv(file,delimiter="\t",dtype={"Protein Name":str,"Uniprot Accession #":str,"k0":str})
            data = data.where(pd.notnull(data), None)
            for i in range(len(data)):
                #update entity collection
                entity = self.build_entity_peptide(data,i)
                query = {"identifiers":{"$elemMatch":{"namespace":"uniprot_id",
                                                      "value":entity["identifiers"][0]["value"]}}}
                self.entity_col.update_one(query,
                                           {"$set": {
                                                 "type": entity["type"],
                                                 "name": entity["name"],                      
                                                 "schema_version": "2.0"},
                                            "$addToSet": {
                                                "identifiers": {"$each": entity["identifiers"]}}},
                                           upsert=True)
                
                obs = self.build_obs_peptide(data,i,file[:file.find(".")].replace("_"," "))
                #update identifier collection
                query = {"$and":[{"namespace":"uniprot_id"},
                                 {"value":entity["identifiers"][0]["value"]}]}
                self.identifier_col.update_one(query,
                                            {"$set": {"namespace": "uniprot_id",
                                                      "value": entity["identifiers"][0]["value"],
                                                      "description": entity["name"]}},upsert=True)

                #update observation collection
                con_1 = {"source":{"$elemMatch":{"namespace":"doi","value":"10.1073/pnas.1006551107"}}}
                con_2 = {"identifier":{"namespace":"uniprot_id","value":entity["identifiers"][0]["value"]}}
                query = {"$and": [con_1,con_2]}
                self.obs_col.update_one(query,
                                    {"$set": {"entity": obs["entity"],
                                              "genotype": obs["genotype"],
                                              "schema_version": "2.0",
                                              "identifier":{"namespace":"uniprot_id","value": entity["identifiers"][0]["value"]}},
                                     "$addToSet": {"values": {"$each": obs["values"]},
                                                   "source": {"$each": obs["source"]}}},
                                     upsert=True)
    def process_docs_protein(self):
        protein_files = ['brain_protein.txt','liver_protein.txt','blood_protein.txt']
        for file in protein_files:
            data = pd.read_csv(file,delimiter="\t",dtype={"Uniprot Accession #":str,"k0":str})
            data = data.where(pd.notnull(data), None)
            for i in range(len(data)):
                #update entity collection
                entity = self.build_entity_protein(data,i)
                query = {"identifiers":{"$elemMatch":{"namespace":"uniprot_id",
                                                      "value":entity["identifiers"][0]["value"]}}}
                self.entity_col.update_one(query,
                                           {"$set": {
                                                 "type": entity["type"],
                                                 "name": entity["name"],                      
                                                 "schema_version": "2.0"},
                                            "$addToSet": {
                                                "identifiers": {"$each": entity["identifiers"]}}},
                                           upsert=True)
                
                obs = self.build_obs_protein(data,i,file[:file.find(".")].replace("_"," "))
                #update identifier collection
                query = {"$and":[{"namespace":"uniprot_id"},
                                 {"value":entity["identifiers"][0]["value"]}]}
                self.identifier_col.update_one(query,
                                            {"$set": {"namespace": "uniprot_id",
                                                      "value": entity["identifiers"][0]["value"],
                                                      "description":entity["name"]}},upsert=True)

                #update observation collection
                con_1 = {"source":{"$elemMatch":{"namespace":"doi","value":"10.1073/pnas.1006551107"}}}
                con_2 = {"identifier":{"namespace":"uniprot_id","value":entity["identifiers"][0]["value"]}}
                query = {"$and": [con_1,con_2]}
                                                         
                self.obs_col.update_one(query,
                                    {"$set": {"entity": obs["entity"],
                                              "genotype": obs["genotype"],
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

    src = parse_column_k0_mouse_brain(MongoDB = MongoDB,
                                      username=username,
                                      password=password,
                                      collection = "observation",
                                      db = "datanator-demo")
    
    #src.process_docs_peptide()
    src.process_docs_protein()
if __name__== '__main__':
    main()
