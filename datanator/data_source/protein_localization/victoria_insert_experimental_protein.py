import pandas as pd
from pymongo import MongoClient
from datanator_query_python.util import mongo_util
from datanator_query_python.config import config

class insert_experimental_protein(mongo_util.MongoUtil):
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
        entity["name"]=str(data.iloc[i,6])
        entity["identifiers"]=[]
        entity["identifiers"].append({"namespace": "uniprot_id",
                "value": str(data.iloc[i,0])})
        entity["identifiers"].append({"namespace": "Refseq_Accession",
                "value":str(data.iloc[i,1])})
        entity["identifiers"].append({"namespace": "Other_Accession",
                "value":str(data.iloc[i,2])})
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
                    "name": "Cytoplasmic protein",
                    "identifiers": [{}... {}]
                },
                "value": [],
                "source": {}, ...
            }
        """
        
        entity={}
        entity["type"]="protein"
        entity["name"]=str(data.iloc[i,6])
        print(entity["name"])
        entity["identifiers"]=[]
        entity["identifiers"].append({"namespace": "uniprot_id",
                "value": str(data.iloc[i,0])})
        entity["identifiers"].append({"namespace":"Refseq_Accession",
                                          "value":str(data.iloc[i,1])})
        entity["identifiers"].append({"namespace": "Other_Accession",
                "value":str(data.iloc[i,2])})
        
        values_p=[]
        exp=str(data.iloc[i,3]).split(",")
        for x in exp:
            values_p.append({"type":"localization","value":x,"description":"experimental"})
        values_p.append({"type":"localization","value":str(data.iloc[i,4]),"description":"secondary"})
        values_p.append({"type":"MultipleSCL","value":str(data.iloc[i,5])})
                         
        genotype = {"taxon": {"ncbi_taxonomy_id":int(data.iloc[i,9]),
                              "name":str(data.iloc[i,10])}}
        genotype["taxon"]["canon_ancestors"]=[]
        query={"tax_id":int(data.iloc[i,9])}
        projection={"_id":0,"canon_anc_ids":1,"canon_anc_names":1}
        doc = self.client["datanator-test"]["taxon_tree"].find_one(filter=query, projection=projection)
        if doc!=None:
            for j in range(len(doc["canon_anc_names"])):
                d={}
                d["ncbi_taxonomy_id"]=doc["canon_anc_ids"][j]
                d["name"]=doc["canon_anc_names"][j]
                genotype["taxon"]["canon_ancestors"].append(d)
        genotype["cellType"]=str(data.iloc[i,13])
        source=[{"namespace":"ePSORTdb","value": "Version "+str(data.iloc[i,17])}]
        ob_p ={"entity":entity,
               "genotype":genotype,
               "values":values_p,
               "source":source,
               "schema_version":"2.0"}
        
        return ob_p
        
                         
    def process_docs(self, skip=0):
        #read .tsv file
        data=pd.read_csv('Experimental-PSORTdb-v4.00.tsv',delimiter="\t")
        data = data.where(pd.notnull(data), None)
        
        #update entity collection
        
        for i in range(len(data)):
            
            if data.iloc[i,0]!="None":
                con_0={"identifiers":{"$elemMatch":{"namespace":"uniprot_id","value":str(data.iloc[i,0])}}}
            elif data.iloc[i,1]!="None":
                con_0 = {"identifiers":{"$elemMatch":{"namespace": "Refseq_Accession",
                "value":str(data.iloc[i,1])}}}
            elif data.iloc[i,2]!="None":
                con_0 = {"identifiers":{"$elemMatch":{"namespace":"Other_Accession","value":str(data.iloc[i,2])}}}
              
            entity = self.build_entity(data,i)
             
            self.entity_col.update_one(con_0,
                                           {"$set": {"type": entity["type"],
                                                     "name": entity["name"],                      
                                                     "schema_version": "2.0"},
                                            "$addToSet": {"identifiers": {"$each": entity["identifiers"]}}},
                                           upsert=True)
            
            
            obs = self.build_obs(data,i)
            
            if obs["entity"]["identifiers"][0]["value"]!="None":
                id_name="uniprot_id"
                value=obs["entity"]["identifiers"][0]["value"]
            else:
                if obs["entity"]["identifiers"][1]["value"]!="None":
                    id_name="Refseq_Accession"
                    value=obs["entity"]["identifiers"][1]["value"]
                elif obs["entity"]["identifiers"][2]["value"]!="None":
                    id_name="Other_Accession"
                    value=obs["entity"]["identifiers"][2]["value"]
                else:
                    id_name="uniprot_id"
                    value="None"
    
                
            query={"$and":[{"namespace":id_name},{"value":value}]}
            con_1={"identifier":{"value":value}}
            #update identifier collection
            self.identifier_col.update_one(query,
                                                {"$set": {"namespace": id_name,
                                                  "value": value,
                                                  "description":str(entity["name"])}},upsert=True)
            
            #update observation collection
            con_2 = {"source":{"$elemMatch":{"namespace":"ePSORTdb","value":"Version "+str(data.iloc[i,17])}}}
            query={"$and":[con_1,con_2]}
            self.obs_col.update_one(query,
                                        {"$set": {"entity": obs["entity"],
                                                  "genotype": obs["genotype"],
                                                  "schema_version": "2.0",
                                                  "identifier":{"namespace":id_name,"value":value}},
                                         "$addToSet": {"values": {"$each": obs["values"]},
                                                       "source": {"$each": obs["source"]}}},
                                         upsert=True)
                       
def main():
    conf=config.Victoria()
    conf_main = config.Config()
    username = conf.USERNAME
    password = conf.PASSWORD
    MongoDB = conf_main.SERVER
    src = insert_experimental_protein(MongoDB = MongoDB,
                                       username=username,
                                       password=password,
                                       collection = "observation",
                                       db = "datanator-demo")
    src.process_docs()
    
if __name__== '__main__':
    main()









        
