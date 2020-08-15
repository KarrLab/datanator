import pandas as pd
import bson
from pymongo import MongoClient
from datanator_query_python.util import mongo_util
from datanator_query_python.config import config

class parse_mouse_global_quantification(mongo_util.MongoUtil):
    def __init__(self,
                 MongoDB=None,
                 db=None,
                 collection=None,
                 username=None,
                 password=None,
                 authSource = 'admin',
                 readPreference = 'nearest'):
        super(parse_mouse_global_quantification,self).__init__(MongoDB=MongoDB, db=db,
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
                    "name": "Histone H2A",
                    "identifiers": [{"namespace": "uniprot_id",
                                     "value": "A2AB79"},
                                    {"namespace": "gene_name",
                                     "value": "H2aw"},
                                    {"namespace": "MGI",
                                     "value": "2448458"},
                                    {"namespace": "RefSeq",
                                     "value": "NP_835736.1"
                                    }
                                   ]
                }
            }        
        """
        entity = {}
        entity["type"] = "protein"
        entity["name"] = data.iloc[i,2]
        entity["identifiers"] = []
        entity["identifiers"].append({"namespace":"uniprot_id","value":data.iloc[i,0]})
        if data.iloc[i,1]!=None:
            gene_names = str(data.iloc[i,1]).split(" ")
            for gene in gene_names:
                entity["identifiers"].append({"namespace":"gene_name","value":gene})
        if data.iloc[i,3]!=None:
            ensembl_ids = str(data.iloc[i,3][:len(str(data.iloc[i,3]))-1]).split(";")
            for ensembl in ensembl_ids:
                entity["identifiers"].append({"namespace":"ensembl","value":ensembl})
        if data.iloc[i,4]!=None:
            entity["identifiers"].append({"namespace":"MGI","value":str(data.iloc[i,4]).replace(";","")})
        if data.iloc[i,5]!=None:
            ref_seq_ids = str(data.iloc[i,5][:len(str(data.iloc[i,5]))-1]).split(";")
            for ref_seq in ref_seq_ids:
                entity["identifiers"].append({"namespace":"RefSeq","value":ref_seq})
        return entity
    
    def build_obs(self, data, i, doc):
        """Build observation objects from obj.
        Go into observations collection.
        Args:
            data (:obj:`Obj`): source object.
            i (:obj: `int`): index (row labels) of object in dataframe
            doc (:obj: `Obj`): MongoDB document containing identifiers 
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
        entity["name"] = doc["name"]
        entity["identifiers"] = []
        for identifier in doc["identifiers"]:
            if identifier["namespace"] == "uniprot_id":
                entity["identifiers"].append({"namespace":"uniprot_id","value":identifier["value"]})
        for identifier in doc["identifiers"]:
            if identifier["namespace"]=="gene_name":
                entity["identifiers"].append({"namespace":"gene_name","value":identifier["value"]})

        values_p = []
        values_p.append({"type": "protein length",
                         "value": data.iloc[i,9],
                         "units": "amino acids"})
        
        values_p.append({"type": "protein molecular weight",
                         "value": data.iloc[i,10],
                         "units": "kDa"})
    
        values_p.append({"type": "protein copy number",
                         "value": data.iloc[i,11],
                         "description": "experiment",
                         "units":"molecules/cell"})
        
        if data.iloc[i,12]!=None:
            values_p.append({"type": "protein copy number",
                             "value": data.iloc[i,12],
                             "description": "replicate",
                             "units":"molecules/cell"})

        values_p.append({"type": "protein copy number",
                         "value": data.iloc[i,13],
                         "description": "average",
                         "units":"molecules/cell"})

        if data.iloc[i,14]!=None:
            values_p.append({"type": "mRNA copy number",
                             "value": data.iloc[i,14],
                             "description": "experiment",
                             "units":"molecules/cell"})

        if data.iloc[i,15]!=None:
            values_p.append({"type": "mRNA copy number",
                             "value": data.iloc[i,15],
                             "description": "replicate",
                             "units":"molecules/cell"})
            
        if data.iloc[i,16]!=None:
            values_p.append({"type": "mRNA copy number",
                             "value": data.iloc[i,16],
                             "description": "average",
                             "units":"molecules/cell"})

        values_p.append({"type": "protein half-life",
                         "value": data.iloc[i,17]*3600,
                         "description": "experiment",
                         "units": "s"})

        if data.iloc[i,18]!=None:
            values_p.append({"type": "protein half-life",
                             "value": float(data.iloc[i,18])*3600,
                             "description": "replicate",
                             "units": "s"})

        values_p.append({"type": "protein half-life",
                         "value": data.iloc[i,19]*3600,
                         "description": "average",
                         "units": "s"})

        if data.iloc[i,20]!=None:
            values_p.append({"type": "mRNA half-life",
                             "value": float(data.iloc[i,20])*3600,
                             "description": "experiment",
                             "units": "s"})

        if data.iloc[i,21]!=None:
            values_p.append({"type": "mRNA half-life",
                             "value": float(data.iloc[i,20])*3600,
                             "description": "replicate",
                             "units": "s"})

        if data.iloc[i,22]!=None:
            values_p.append({"type": "mRNA half-life",
                             "value": float(data.iloc[i,22])*3600,
                             "description": "average",
                             "units": "s"})
            
        if data.iloc[i,23]!=None:
            values_p.append({"type": "transcription rate (vsr)",
                             "value": data.iloc[i,23],
                             "description": "experiment",
                             "units": "molecules/(cell*h)"})

        if data.iloc[i,24]!=None:
            values_p.append({"type": "transcription rate (vsr)",
                             "value": data.iloc[i,24],
                             "description": "replicate",
                             "units": "molecules/(cell*h)"})

        if data.iloc[i,25]!=None:
            values_p.append({"type": "transcription rate (vsr)",
                             "value": data.iloc[i,25],
                             "description": "average",
                             "units": "molecules/(cell*h)"})

        if data.iloc[i,26]!=None:
            values_p.append({"type": "translation rate constant (ksp)",
                             "value": data.iloc[i,26],
                             "description": "experiment",
                             "units": "molecules/(mRNA*h)"})
            
        if data.iloc[i,27]!=None:
            values_p.append({"type": "translation rate constant (ksp)",
                             "value": data.iloc[i,27],
                             "description": "replicate",
                             "units": "molecules/(mRNA*h)"})

        if data.iloc[i,28]!=None:
            values_p.append({"type": "translation rate constant (ksp)",
                             "value": data.iloc[i,28],
                             "description": "average",
                             "units": "molecules/(mRNA*h)"})

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
        genotype["cellType"] = "fibroblast"
        genotype["cellLine"] = "NIH3T3"

        environment = {}
        environment["media"] = "Dulbecco’s Modified Eagle’s Medium (DMEM) Glutamax lacking arginine and lysine, supplemented with 10% dialyzed fetal bovine serum"
        environment["temperature"] = bson.Int64(37)
        
        source = [{"namespace":"doi","value": "10.1038/nature10098"}]
        ob_p = {"entity":entity,
            "genotype":genotype,
            "environment":environment,
            "values":values_p,
            "source":source,
            "schema_version":"2.0"}
        return ob_p
    
    def process_docs(self):
        id_file= pd.read_csv('mouse_uniprot_mapped.tab',delimiter='\t')
        id_file = id_file.where(pd.notnull(id_file), None)
        for i in range(len(id_file)):
            #update entity collection
            entity = self.build_entity(id_file,i)
            query = {"identifiers":{"$elemMatch":{"namespace":"uniprot_id",
                                              "value":id_file.iloc[i,0]}}}
            self.entity_col.update_one(query,
                                           {"$set": {
                                                 "type": entity["type"],
                                                 "name": entity["name"],                      
                                                 "schema_version": "2.0"},
                                            "$addToSet": {
                                                "identifiers": {"$each": entity["identifiers"]}}},
                                           upsert=True)
            #update identifier collection
            query = {"$and":[{"namespace":"uniprot_id"},{"value":id_file.iloc[i,0]}]}
            self.identifier_col.update_one(query,
                                            {"$set": {"namespace": "uniprot_id",
                                              "value": id_file.iloc[i,0],
                                              "description":entity["name"]}},upsert=True)
        data = pd.read_csv('table3_global_quantification_mouse.txt',
                           delimiter='\t',
                           dtype={'Protein length [amino acids]':str,
                                  'Protein molecular weight [kDa]':str,
                                  'Protein copy number experiment [molecules/cell]':str,
                                  'Protein copy number replicate [molecules/cell]':str,
                                  'Protein copy number average [molecules/cell]':str,
                                  'mRNA copy number experiment [molecules/cell]':str,
                                  'mRNA copy number replicate [molecules/cell]':str,
                                  'mRNA copy number average [molecules/cell]':str,
                                  'Protein half-life experiment [h]':str,
                                  'Protein half-life replicate [h]':str,
                                  'Protein half-life average [h]':str,
                                  'mRNA half-life experiment [h]':str,
                                  'mRNA half-life replicate [h]':str,
                                  'mRNA half-life average [h]':str,
                                  'transcription rate (vsr) experiment [molecules/(cell*h)]':str,
                                  'transcription rate (vsr) replicate [molecules/(cell*h)]':str,
                                  'transcription rate (vsr) average [molecules/(cell*h)]':str,
                                  'translation rate constant (ksp) experiment [molecules/(mRNA*h)]':str,
                                  'translation rate constant (ksp) replicate [molecules/(mRNA*h)]':str,
                                  'translation rate constant (ksp) average [molecules/(mRNA*h)]':str})
        data = data.where(pd.notnull(data), None)
        for i in range(len(data)):
            if data.iloc[i,4]!=None:
                uniprot_ids = str(data.iloc[i,4]).split(";")
                for uniprot in uniprot_ids:
                    query = {"identifiers.value": uniprot}
                    projection = {"_id":0, "name":1, "identifiers":1}
                    doc = self.client["datanator-demo"]["entity"].find_one(filter = query, projection = projection)
                    if doc!=None:
                        obs = self.build_obs(data, i, doc)
                        con_1 = {"source":{"$elemMatch":{"namespace":"doi","value":"10.1038/nature10098"}}}
                        con_2 = {"identifier":{"namespace":"uniprot_id","value":obs["entity"]["identifiers"][0]["value"]}}
                        query = {"$and": [con_1,con_2]}
                        self.obs_col.update_one(query,
                                            {"$set": {"entity": obs["entity"],
                                                      "genotype": obs["genotype"],
                                                      "environment": obs["environment"],
                                                      "schema_version": "2.0",
                                                      "identifier":{"namespace":"uniprot_id","value": obs["entity"]["identifiers"][0]["value"]}},
                                             "$addToSet": {"values": {"$each": obs["values"]},
                                                           "source": {"$each": obs["source"]}}},
                                             upsert=True)
        
def main():
    conf=config.Victoria()
    conf_main = config.Config()
    username = conf.USERNAME
    password = conf.PASSWORD
    MongoDB = conf_main.SERVER

    src = parse_mouse_global_quantification(MongoDB = MongoDB,
                                            username=username,
                                            password=password,
                                            collection = "observation",
                                            db = "datanator-demo")
    
    src.process_docs()
if __name__== '__main__':
    main()



        
