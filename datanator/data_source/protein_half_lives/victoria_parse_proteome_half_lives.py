import pandas as pd
from pymongo import MongoClient
from datanator_query_python.util import mongo_util
from datanator_query_python.config import config

class parse_proteome_half_life(mongo_util.MongoUtil):
    def __init__(self,
                 MongoDB=None,
                 db=None,
                 collection=None,
                 username=None,
                 password=None,
                 authSource = 'admin',
                 readPreference = 'nearest'):
        super(parse_proteome_half_life,self).__init__(MongoDB=MongoDB, db=db,
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
            i (:obj: `int`): index (row labels) of object in dataframe
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
        entity["name"] = str(data.iloc[i,1]).replace("?","-")
        entity["identifiers"] = []
        if i==39:
            entity["identifiers"].append({"namespace": "gene_symbol",
                                          "value":"K-"+"\u03B1"+"-1"})
        else:
            entity["identifiers"].append({"namespace": "gene_symbol",
                                          "value":str(data.iloc[i,0])})
        return entity
    
    def build_obs(self, data, i, table):
        """Build observation objects from obj. (Tables S4, S6-S9) 
        Go into observations collection.
        Args:
            data (:obj:`Obj`): source object.
            i (:obj: `int`): index (row labels) of object in dataframe
            table (:obj: `int`): table number
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
        if table == 4:
            entity["name"] = str(data.iloc[i,1]).replace("?","-")
            if i==39:
                entity["identifiers"].append({"namespace": "gene_symbol",
                                              "value":"K-"+"\u03B1"+"-1"})
            else:
                entity["identifiers"].append({"namespace": "gene_symbol",
                                              "value":str(data.iloc[i,0])})
        else:
            entity["identifiers"].append({"namespace": "gene_symbol",
                                          "value":str(data.iloc[i,0])})
            query = {"identifiers":{"$elemMatch":{"value":entity["identifiers"][0]["value"]}}}
            projection = {"_id":0,"name":1}
            doc = self.client["datanator-demo"]["entity"].find_one(filter = query, projection = projection)
            if doc!=None:
                entity["name"]=doc["name"]

        values_p = []
        if table == 4:
            values_p.append({"type":"Half-life",
                             "value":str(float(data.iloc[i,2])*3600),
                             "uncertainty":float(data.iloc[i,3])*3600,
                             "units":"s"})
        else:
            v = ""
            if ">" in str(data.iloc[i,1]):
                v = "greater than " + str(float(data.iloc[i,1][1:])*3600)
            else:
                v = str(float(data.iloc[i,1])*3600)
            if data.iloc[i,2]!=None:
                values_p.append({"type":"Half-life",
                                 "value":v,
                                 "uncertainty":float(data.iloc[i,2])*3600,
                                 "units":"s"})
            else:
                values_p.append({"type":"Half-life",
                                 "value":v,
                                 "units":"s"})
        
        genotype = {"taxon": {"ncbi_taxonomy_id":int(9606),
                              "name":"Homo sapiens"}}
        genotype["taxon"]["canon_ancestors"] = []
        query = {"tax_id":int(9606)}
        projection = {"_id":0,"canon_anc_ids":1,"canon_anc_names":1}
        doc = self.client["datanator-test"]["taxon_tree"].find_one(filter=query, projection=projection)
        if doc!=None:
            for j in range(len(doc["canon_anc_names"])):
                d = {}
                d["ncbi_taxonomy_id"] = doc["canon_anc_ids"][j]
                d["name"] = doc["canon_anc_names"][j]
                genotype["taxon"]["canon_ancestors"].append(d)
        
        genotype["organ"]="lung"
        genotype["cellType"]="non‐small cell lung carcinoma"
        genotype["cellLine"]="H1299"
        
        environment = {}
        if table == 4:
            environment["condition"] = "normal growth condition"
            environment["media"] = "2 ml, without Cyclohexamide (CHX), RPMI 1640, 0.05% Penicillin‐Streptomycin antibiotics, 10% Fetal Calf Serum (FCS), L‐Glutamine, no riboflavin and phenol red"
        elif table == 6:
            environment["condition"] = "after paclitaxel addition"
            environment["media"] = "2 ml, 1 μM paclitaxel, RPMI 1640, 0.05% Penicillin‐Streptomycin antibiotics, 10% Fetal Calf Serum (FCS), L‐Glutamine, no riboflavin and phenol red"
        elif table == 7:
            environment["condition"] = "during serum starvation"
            environment["media"] = "5 ml, methionine/cysteine‐free RPMI 1640, 10% dialyzed FBS, 2mM L‐glutamine"
        elif table == 8:
            environment["condition"] = "after cisplatin addition"
            environment["media"] = "2 ml, 80 μM  cisplatin, RPMI 1640, 0.05% Penicillin‐Streptomycin antibiotics, 10% Fetal Calf Serum (FCS), L‐Glutamine, no riboflavin and phenol red"
        elif table == 9:
            environment["condition"] = "after transcription inhibitor actinomycin‐D addition"
            environment["media"] = "2 ml, 1 ug/ml actinomycin‐D diluted from stock, RPMI 1640, 0.05% Penicillin‐Streptomycin antibiotics, 10% Fetal Calf Serum (FCS), L‐Glutamine, no riboflavin and phenol red"
            
        source = [{"namespace":"doi","value":"10.1126/science.1199784"}]
        
        ob_p = {"entity":entity,
                "genotype":genotype,
                "environment":environment,
                "values":values_p,
                "source":source,
                "schema_version":"2.0"}
        print(str(i)+"\t"+entity["identifiers"][0]["value"])
        return ob_p

    def build_obs_s5(self, data, i, cpt):
        """Build observation objects from obj. (Tables S5) 
        Go into observations collection.
        Args:
            data (:obj:`Obj`): source object.
            i (:obj: `int`): index (row labels) of object in dataframe
            cpt (:obj: `bool`): false before CPT addition, true after CPT addition
            
        Return:
            obj(:obj:`Obj`)
            {
                "identifier":{"namespace": "gene_symbol",
                              "value": "BAG1"
                },
                "entity": {
                    "type": "protein",
                    "name": "BAG1",
                    "identifiers": [
                                    {"namespace": "gene_symbol",
                                     "value: "BCL2-associated athanogene"
                                    }
                                   ]
                },
                "genotype":{
                    "taxon": {},
                    "organ": "lung",
                    "cellType": "non‐small cell lung carcinoma",
                    "cellLine": "H1299"
                },
                "environment":{
                    "condition": "after CPT addition",
                    "media": "2 ml, 10 μM camptothecin, RPMI 1640, 0.05% Penicillin‐Streptomycin antibiotics, 10% Fetal Calf Serum (FCS), L‐Glutamine, no riboflavin and phenol red"
                },
                "values": [
                    {"type": "Half-life",
                     "value": "34200.0",
                     "uncertainty": "5400",
                     "units": "s"
                    }
                ]
                "source": [
                            {"namespace": "doi",
                             "value": "10.1126/science.1199784"
                            }
                ]
            }
        """
        
        entity = {} 
        entity["type"] = "protein"
        
        entity["identifiers"] = []
        if i==14:
            entity["identifiers"].append({"namespace": "gene_symbol",
                                          "value":"K-"+"\u03B1"+"-1"})
        else:
            entity["identifiers"].append({"namespace": "gene_symbol",
                                          "value":str(data.iloc[i,0])})
        query = {"identifiers":{"$elemMatch":{"value":entity["identifiers"][0]["value"]}}}
        projection = {"_id":0,"name":1}
        doc = self.client["datanator-demo"]["entity"].find_one(filter = query, projection = projection)
        if doc!=None:
            entity["name"]=doc["name"]
            
        values_p = []
        if (cpt!=True):
            values_p.append({"type":"Half-life",
                             "value":str(float(data.iloc[i,1])*3600),
                             "uncertainty":float(data.iloc[i,2])*3600,
                             "units":"s"})
        else:
            v=""
            if ">" in str(data.iloc[i,3]):
                v = "greater than " + str(float(data.iloc[i,3][1:])*3600)
            else:
                v = str(float(data.iloc[i,3])*3600)
            if data.iloc[i,4]!=None:
                values_p.append({"type":"Half-life",
                                 "value":v,
                                 "uncertainty":float(data.iloc[i,4])*3600,
                                 "units":"s"})
            else:
                values_p.append({"type":"Half-life",
                                 "value":v,
                                 "units":"s"})
        
        genotype = {"taxon": {"ncbi_taxonomy_id":int(9606),
                              "name":"Homo sapiens"}}
        genotype["taxon"]["canon_ancestors"] = []
        query = {"tax_id":int(9606)}
        projection = {"_id":0,"canon_anc_ids":1,"canon_anc_names":1}
        doc = self.client["datanator-test"]["taxon_tree"].find_one(filter = query, projection = projection)
        if doc!=None:
            for j in range(len(doc["canon_anc_names"])):
                d = {}
                d["ncbi_taxonomy_id"] = doc["canon_anc_ids"][j]
                d["name"] = doc["canon_anc_names"][j]
                genotype["taxon"]["canon_ancestors"].append(d)
        
        genotype["organ"]="lung"
        genotype["cellType"]="non‐small cell lung carcinoma"
        genotype["cellLine"]="H1299"
        
        environment = {}
        if (cpt):
            environment["condition"] = "after CPT addition"
            environment["media"] = "2 ml, 10 μM camptothecin (CPT), RPMI 1640, 0.05% Penicillin‐Streptomycin antibiotics, 10% Fetal Calf Serum (FCS), L‐Glutamine, no riboflavin and phenol red"
        else:
            environment["condition"] = "before CPT addition"
            environment["media"] = "2 ml, without camptothecin (CPT), RPMI 1640, 0.05% Penicillin‐Streptomycin antibiotics, 10% Fetal Calf Serum (FCS), L‐Glutamine, no riboflavin and phenol red"

        source = [{"namespace":"doi","value":"10.1126/science.1199784"}]
        
        ob_p = {"entity":entity,
                "genotype":genotype,
                "environment":environment,
                "values":values_p,
                "source":source,
                "schema_version":"2.0"}
        print(str(i)+"\t"+entity["identifiers"][0]["value"])
        return ob_p
    
    def process_docs(self):
            #read Table S4
            data = pd.read_csv('table_S4.txt',delimiter="\t")
            data = data.where(pd.notnull(data), None)
            con_1 = {"source":{"$elemMatch":{"namespace":"doi","value":"10.1126/science.1199784"}}}
            for i in range(len(data)):
                #update entity collection
                entity = self.build_entity(data,i)
                query = {"identifiers":{"$elemMatch":{"namespace":"gene_symbol","value":entity["identifiers"][0]["value"]}}}
                self.entity_col.update_one(query,
                                               {"$set": {
                                                     "type": entity["type"],
                                                     "name": entity["name"],                      
                                                     "schema_version": "2.0"},
                                                "$addToSet": {
                                                    "identifiers": {"$each": entity["identifiers"]}}},
                                               upsert=True)
                
                obs_s4 = self.build_obs(data,i,4)
                #update identifier collection
                query = {"$and":[{"namespace":"gene_symbol"},{"value":entity["identifiers"][0]["value"]}]}
                self.identifier_col.update_one(query,
                                                {"$set": {"namespace": "gene_symbol",
                                                          "value": entity["identifiers"][0]["value"],
                                                          "description":str(entity["name"])}},upsert=True)
                
                #update observation collection
                con_2 = {"identifier":{"namespace":"gene_symbol","value":entity["identifiers"][0]["value"]}}
                query = {"$and":[con_1,con_2]}
                self.obs_col.update_one(query,
                                        {"$set": {"entity": obs_s4["entity"],
                                                  "genotype": obs_s4["genotype"],
                                                  "environment":obs_s4["environment"],
                                                  "schema_version": "2.0",
                                                  "identifier":{"namespace":"gene_symbol","value": entity["identifiers"][0]["value"]}},
                                         "$addToSet": {"values": {"$each": obs_s4["values"]},
                                                       "source": {"$each": obs_s4["source"]}}},
                                         upsert=True)
               
            #read Table S5
            data_s5 = pd.read_csv('table_S5.txt',delimiter="\t")
            data_s5 = data_s5.where(pd.notnull(data_s5), None)
            for i in range(len(data_s5)):
                obs_s5 = self.build_obs_s5(data_s5,i,False)
                con_2 = {"identifier":{"namespace":"gene_symbol","value":obs_s5["entity"]["identifiers"][0]["value"]}}
                con_3 = {"environment":{"$elemMatch":{"condition":obs_s5["environment"]["condition"]}}}
                query = {"$and":[con_1,con_2,con_3]}
                self.obs_col.update_one(query,
                                        {"$set": {"entity": obs_s5["entity"],
                                                  "genotype": obs_s5["genotype"],
                                                  "environment":obs_s5["environment"],
                                                  "schema_version": "2.0",
                                                  "identifier":{"namespace":"gene_symbol","value": obs_s5["entity"]["identifiers"][0]["value"]}},
                                         "$addToSet": {"values": {"$each": obs_s5["values"]},
                                                       "source": {"$each": obs_s5["source"]}}},
                                         upsert=True)
                
                obs_s5 = self.build_obs_s5(data_s5,i,True)
                con_3 = {"environment":{"$elemMatch":{"condition":obs_s5["environment"]["condition"]}}}
                query = {"$and":[con_1,con_2,con_3]}
                self.obs_col.update_one(query,
                                        {"$set": {"entity": obs_s5["entity"],
                                                  "genotype": obs_s5["genotype"],
                                                  "environment":obs_s5["environment"],
                                                  "schema_version": "2.0",
                                                  "identifier":{"namespace":"gene_symbol","value": obs_s5["entity"]["identifiers"][0]["value"]}},
                                         "$addToSet": {"values": {"$each": obs_s5["values"]},
                                                       "source": {"$each": obs_s5["source"]}}},
                                         upsert=True)

            #read Tables S6-S9
            for k in range(6,10):
                data = pd.read_csv('table_S'+str(k)+'.txt',delimiter="\t")
                data = data.where(pd.notnull(data), None)
                for i in range(len(data)):
                    obs = self.build_obs(data,i,k)
                    con_2 = {"identifier":{"namespace":"gene_symbol","value":obs["entity"]["identifiers"][0]["value"]}}
                    con_3 = {"environment":{"$elemMatch":{"condition":obs["environment"]["condition"]}}}
                    query = {"$and":[con_1,con_2,con_3]}
                    self.obs_col.update_one(query,
                                            {"$set": {"entity": obs["entity"],
                                                      "genotype": obs["genotype"],
                                                      "environment":obs["environment"],
                                                      "schema_version": "2.0",
                                                      "identifier":{"namespace":"gene_symbol","value": obs["entity"]["identifiers"][0]["value"]}},
                                             "$addToSet": {"values": {"$each": obs["values"]},
                                                           "source": {"$each": obs["source"]}}},
                                             upsert=True)
            
def main():
    conf=config.Victoria()
    conf_main = config.Config()
    username = conf.USERNAME
    password = conf.PASSWORD
    MongoDB = conf_main.SERVER

    src = parse_proteome_half_life(MongoDB = MongoDB,
                                   username=username,
                                   password=password,
                                   collection = "observation",
                                   db = "datanator-demo")
    
    src.process_docs()
if __name__== '__main__':
    main()



