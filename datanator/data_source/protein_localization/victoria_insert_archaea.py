import pandas as pd
import numpy as np
from pymongo import MongoClient
from datanator_query_python.util import mongo_util
from datanator_query_python.config import config

class insert_archaea(mongo_util.MongoUtil):
    def __init__(self,
                 MongoDB=None,
                 db=None,
                 collection=None,
                 max_entries=float('inf'),
                 username=None,
                 password=None,
                 authSource = 'admin',
                 readPreference = 'nearest'):
        super().__init__(MongoDB=MongoDB, db=db,
                         username = username,
                         password = password,
                         authSource = authSource,
                         readPreference=readPreference)
        self.max_entries = max_entries
        self.collection = collection
        self.entity_col = self.db_obj["entity"]
        self.identifier_col = self.db_obj["identifier"]
        self.obs_col = self.db_obj["observation"]
        
    def build_entity(self, data, i):
        """Build entity object from data.
        Go into entity collection

        Args:
            data (:obj:`Obj`): source object.
            i (:obj: `int`): index (row labels) of object in dataframe
        Return:
            (:obj:`Obj`), e.g.
            {
                "entity": {
                    "type": "protein",
                    "name": "GNAT family N-acetyltransferase",
                    "identifiers": [{}... {}]
                }
            }        
        """
        entity = {}
        entity["type"] = "protein"
        entity["name"]=str(data.iloc[i,0])[str(data.iloc[i,0]).rfind("|")+2:]
        
        entity["identifiers"]=[]
        entity["identifiers"].append({"namespace": "Seq_ID",
                                      "value": str(data.iloc[i,0])[str(data.iloc[i,0]).find("W"):str(data.iloc[i,0]).rfind("|")]})
        return entity
        

    def build_obs(self, data, i):
        """Build observation objects from obj.
        Go into observations collection.
        Args:
            data (:obj:`Obj`): source object.
            i (:obj: `int`): index (row labels) of object in dataframe
        Return:
            obj(:obj:`Obj`)
            {
                "entity": {
                    "type": "protein",
                    "name": "GNAT family N-acetyltransferase",
                    "identifiers": [{}... {}]
                },
                "values": [],
                "source": {}, ...
            }
        """
        
        entity={}
        entity["type"]="protein"
        entity["name"]=str(data.iloc[i,0])[str(data.iloc[i,0]).rfind("|")+2:]
        
        entity["identifiers"]=[]
        entity["identifiers"].append({"namespace": "Seq_ID",
                                      "value": str(data.iloc[i,0])[str(data.iloc[i,0]).find("W"):str(data.iloc[i,0]).rfind("|")]})
        print(str(data.iloc[i,0])[str(data.iloc[i,0]).find("W"):])

        values_p=[]
        if data.iloc[i,27]==None:
            values_p.append({"type":"localization_score","value":float(data.iloc[i,1]),"description":"Cytoplasmic Membrane"})
            values_p.append({"type":"localization_score","value":float(data.iloc[i,2]),"description":"Cell Wall"})
            values_p.append({"type":"localization_score","value":float(data.iloc[i,3]),"description":"Extracellular"})                 
            values_p.append({"type":"localization_score","value":float(data.iloc[i,4]),"description":"Cytoplasmic"})
            values_p.append({"type":"localization","value":str(data.iloc[i,5]),"description":"predicted"})

        else:
            values_p.append({"type":"localization_score","value":float(data.iloc[i,27]),"description":"Cytoplasmic Membrane"})
            values_p.append({"type":"localization_score","value":float(data.iloc[i,28]),"description":"Cell Wall"})
            values_p.append({"type":"localization_score","value":float(data.iloc[i,29]),"description":"Extracellular"})                 
            values_p.append({"type":"localization_score","value":float(data.iloc[i,30]),"description":"Cytoplasmic"})
            values_p.append({"type":"localization","value":str(data.iloc[i,31]),"description":"predicted"})
        
        genotype={}
        genotype["cellType"]="Archaea"
        source=[{"namespace":"PSORTb","value": "Version 3.0"}]
        ob_p ={"entity":entity,
               "genotype":genotype,
               "values":values_p,
               "source":source,
               "schema_version":"2.0"}
        return ob_p
        
                         
    def process_docs(self, skip=0):
        #read .tab file in chunks
        f=pd.read_csv('Computed-Archaea-PSORTdb-3.00.tab',delimiter="\t",chunksize=100000,low_memory=False)
        chunks = []
        for chunk in f:
            chunks.append(chunk)
        data = pd.concat(chunks)

        for i in range(len(data)):
            #update entity collection
            con_0={"identifiers":{"$elemMatch":{"namespace":"Seq_ID","value":str(data.iloc[i,0])[str(data.iloc[i,0]).find("W"):str(data.iloc[i,0]).rfind("|")]}}}
            entity = self.build_entity(data,i)
            self.entity_col.update_one(con_0,{"$set": {"type": entity["type"],
                                                         "name": entity["name"],                      
                                                         "schema_version": "2.0"},
                                                "$addToSet": {"identifiers": {"$each": entity["identifiers"]}}},
                                               upsert=True)  
            #update identifier collection
            con_2 = {"namespace":"Seq_ID"}
            con_3 = {"value":str(data.iloc[i,0])[str(data.iloc[i,0]).find("W"):str(data.iloc[i,0]).rfind("|")]}
            query={"$and":[con_2,con_3]}
            obs = self.build_obs(data,i)
            self.identifier_col.update_one(query,
                                            {"$set": {"namespace": (obs["entity"]["identifiers"][0])["namespace"],
                                              "value": (obs["entity"]["identifiers"][0])["value"]}},upsert=True)

            #update observation collection
            con_4 = {"identifier":{"namespace":"Seq_ID","value":str(data.iloc[i,0])[str(data.iloc[i,0]).find("W"):str(data.iloc[i,0]).rfind("|")]}}
            con_5 = {"source":{"$elemMatch":{"namespace":"PSORTb","value":"Version 3.0"}}}
            query={"$and":[con_4,con_5]}
            self.obs_col.update_one(query,
                                        {"$set": {"entity": obs["entity"],
                                                  "genotype":obs["genotype"],
                                                  "schema_version": "2.0",
                                                  "identifier": obs["entity"]["identifiers"][0]},
                                         "$addToSet": {"values": {"$each": obs["values"]},
                                                       "source": {"$each": obs["source"]}}},
                                         upsert=True)
        
def main():
    conf=config.Victoria()
    conf_main = config.Config()
    username = conf.USERNAME
    password = conf.PASSWORD
    MongoDB = conf_main.SERVER
    src = insert_archaea(MongoDB = MongoDB,
                                       username=username,
                                       password=password,
                                       collection = "observation",
                                       db = "datanator-demo")
    src.process_docs()
    
if __name__== '__main__':
    main()









        
