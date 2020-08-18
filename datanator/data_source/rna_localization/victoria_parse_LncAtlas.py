import pandas as pd
from datetime import datetime
from pymongo import MongoClient
from datanator_query_python.util import mongo_util
from datanator_query_python.config import config

class parse_LncAtlas(mongo_util.MongoUtil):
    def __init__(self,
                 MongoDB=None,
                 db=None,
                 collection=None,
                 username=None,
                 password=None,
                 authSource = 'admin'):
        super(parse_LncAtlas,self).__init__(MongoDB=MongoDB, db=db,
                                            username = username,
                                            password = password,
                                            authSource = authSource)
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
                    "type": "RNA",
                    "name": "TSPAN6",
                    "identifiers": [{}... {}]
                }
            }        
        """
        entity = {}
        entity["type"] = "RNA"
        entity["name"] = data.iloc[i,4]
        entity["identifiers"] = []
        query = {"$and":[{"gene_name":data.iloc[i,4]},{"ncbi_taxonomy_id":9606}]}
        projection = {"_id":0,"uniprot_id":1}
        doc = self.client["datanator"]["uniprot"].find_one(filter = query, projection = projection)
        if doc!=None:
            entity["identifiers"].append({"namespace":"uniprot_id","value":doc["uniprot_id"]})
        else:
            entity["identifiers"].append({"namespace":"ensembl","value":data.iloc[i,0]})
        entity["identifiers"].append({"namespace":"gene_name","value":data.iloc[i,4]})
        entity["related"] = []
        if data.iloc[i,5]=="coding":
            entity["related"].append({"namespace":"coding_type",
                                      "value":"protein_coding"})
        else:
            entity["related"].append({"namespace":"coding_type",
                                      "value":"non-coding"})
        print(str(i)+" "+entity["name"])
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
                    "type": "RNA",
                    "name": "TSPAN6",
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
        entity["type"] = "RNA"
        entity["name"] = data.iloc[i,4]
        entity["identifiers"] = []
        query = {"$and":[{"gene_name":data.iloc[i,4]},{"ncbi_taxonomy_id":9606}]}
        projection = {"_id":0,"uniprot_id":1}
        doc = self.client["datanator"]["uniprot"].find_one(filter = query, projection = projection)
        if doc!=None:
            entity["identifiers"].append({"namespace":"uniprot_id","value":doc["uniprot_id"]})
        else:
            entity["identifiers"].append({"namespace":"ensembl","value":data.iloc[i,0]})
        entity["identifiers"].append({"namespace":"gene_name","value":data.iloc[i,4]})

        values_p = []
        value_type = "expression_value"
        if data.iloc[i,2]=="ratio2":
            value_type = "Cytoplasmic–Nuclear Relative Concentration Index"
        elif data.iloc[i,2]=="ratioKc":
            value_type = "Chromatin-Subnuclear Relative Concentration Index"
        elif data.iloc[i,2]=="ratioKin":
            value_type = "Insoluble-Subcytoplasmic Relative Concentration Index"
        elif data.iloc[i,2]=="ratioKmem":
            value_type = "Membrane-Subcytoplasmic Relative Concentration Index"
        elif data.iloc[i,2]=="ratioKno":
            value_type = "Nucleolus-Subnuclear Relative Concentration Index"
        elif data.iloc[i,2]=="ratioKnp":
            value_type = "Nucleoplasm-Subnuclear Relative Concentration Index"
        if "ratio" in data.iloc[i,2]:
            values_p.append({"type":value_type,
                             "value":data.iloc[i,3]})
        else:
            description = data.iloc[i,2]
            if description=="cytosolsub":
                description = "subcytoplasmic"
            elif description=="nucleussub":
                description = "subnuclear"
            values_p.append({"type":value_type,
                             "value":data.iloc[i,3],
                             "units":"FPKM (fragments per kilobase per million mapped)",
                             "description":description})
            
        genotype = {"taxon": {"ncbi_taxonomy_id":9606,
                              "name": "Homo sapiens"}}
        genotype["taxon"]["canon_ancestors"] = []
        query = {"tax_id":9606}
        projection = {"_id":0,"canon_anc_ids":1,"canon_anc_names":1}
        doc = self.client["datanator-test"]["taxon_tree"].find_one(filter = query, projection = projection)
        if doc!=None:
            for j in range(len(doc["canon_anc_names"])):
                d = {}
                d["ncbi_taxonomy_id"] = doc["canon_anc_ids"][j]
                d["name"] = doc["canon_anc_names"][j]
                genotype["taxon"]["canon_ancestors"].append(d)
        genotype["cellLine"] = data.iloc[i,1]
        source = []
        source.append({"namespace":"doi","value":"10.1261/rna.060814.117"})
        source.append({"namespace":"LncAtlas","value":datetime.now().isoformat(timespec='minutes')})
        ob_p = {"entity":entity,
                "genotype":genotype,
                "values":values_p,
                "source":source,
                "schema_version":"2.0"}
        print(str(i)+" "+entity["name"])
        return ob_p

    def process_docs(self):
        data = pd.read_csv('lncAtlas.csv',dtype={"ENSEMBL ID":str,
                                                 "Data Source":str,
                                                 "Data Type":str,
                                                 "Value":str,
                                                 "Gene Name":str,
                                                 "Coding Type":str})
        data = data.where(pd.notnull(data), None)
        for i in range(len(data)):
            #update identifier collection
            query_id = {"$and":[{"gene_name":data.iloc[i,4]},{"ncbi_taxonomy_id":9606}]}
            query_identifier = {}
            query_entity = {}
            query_cellLine = {"genotype.cellLine":data.iloc[i,1]}
            projection = {"_id":0,"uniprot_id":1}
            doc = self.client["datanator"]["uniprot"].find_one(filter = query_id, projection = projection)
            entity = self.build_entity(data, i)
            obs = self.build_obs(data, i)
            cond_source = {"source":{"$elemMatch":{"namespace":"doi","value":"10.1261/rna.060814.117"}}}
            if doc!=None:
                #uniprot_id is in MongoDB's uniprot collection
                query_identifier = {"$and":[{"namespace":"uniprot_id"},
                                 {"value":doc["uniprot_id"]}]}
                #update identifier collection
                self.identifier_col.update_one(query_identifier,
                                               {"$set": {"namespace": "uniprot_id",
                                                "value": doc["uniprot_id"],
                                                "description": entity["name"]}},upsert=True)
                query_entity = {"identifiers":{"$elemMatch":{"namespace":"uniprot_id",
                                                             "value": doc["uniprot_id"]}}}
                #update observation collection
                query_identifier = {"identifier":{"namespace":"uniprot_id","value":doc["uniprot_id"]}}
                query_obs = {"$and":[cond_source,query_identifier,query_cellLine]}
                self.obs_col.update_one(query_obs,
                                        {"$pull":{"source":{"namespace":"LncAtlas"}}})
                self.obs_col.update_one(query_obs,
                                        {"$set": {"entity": obs["entity"],
                                                  "genotype": obs["genotype"],
                                                  "schema_version": "2.0",
                                                  "identifier":{"namespace":"uniprot_id","value": doc["uniprot_id"]}},
                                         "$addToSet": {"values": {"$each": obs["values"]},
                                                       "source": {"$each": obs["source"]}}},
                                         upsert=True)
            else:
                #can't find uniprot_id in MongoDB's uniprot collection
                query_identifier = {"$and":[{"namespace":"ensembl"},
                                            {"value":data.iloc[i,0]}]}
                #update identifier collection
                self.identifier_col.update_one(query_identifier,
                                               {"$set": {"namespace": "ensembl",
                                                "value": data.iloc[i,0],
                                                "description": entity["name"]}},upsert=True)
                query_entity = {"identifiers":{"$elemMatch":{"namespace":"ensembl",
                                                             "value":data.iloc[i,0]}}}
                #update observation collection
                query_identifier = {"identifier":{"namespace":"ensembl","value":data.iloc[i,0]}}
                query_obs = {"$and":[cond_source,query_identifier,query_cellLine]}
                self.obs_col.update_one(query_obs,
                                        {"$pull":{"source":{"namespace":"LncAtlas"}}})
                self.obs_col.update_one(query_obs,
                                        {"$set": {"entity": obs["entity"],
                                                  "genotype": obs["genotype"],
                                                  "schema_version": "2.0",
                                                  "identifier":{"namespace":"ensembl","value": data.iloc[i,0]}},
                                         "$addToSet": {"values": {"$each": obs["values"]},
                                                       "source": {"$each": obs["source"]}}},
                                         upsert=True)
            #update entity collection
            self.entity_col.update_one(query_entity,
                                       {"$set": {
                                             "type": entity["type"],
                                             "name": entity["name"],                      
                                             "schema_version": "2.0"},
                                        "$addToSet": {
                                            "identifiers": {"$each": entity["identifiers"]}}},
                                       upsert=True)
def main():
    conf=config.Victoria()
    conf_main = config.Config()
    username = conf.USERNAME
    password = conf.PASSWORD
    MongoDB = conf_main.SERVER

    src = parse_LncAtlas(MongoDB = MongoDB,
                         username=username,
                         password=password,
                         collection = "observation",
                         db = "datanator-demo")
    src.process_docs()
       
if __name__== '__main__':
    main()
