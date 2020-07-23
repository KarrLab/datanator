import pandas as pd
import urllib.request
import json
import numpy as np
from pymongo import MongoClient
from datanator_query_python.util import mongo_util
from datanator_query_python.config import config

class insert_protein_atlas(mongo_util.MongoUtil):
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
        
    def build_entity(self, data):
        """Build entity object from obj.
        Go into entity collection

        Args:
            data (:obj:`Obj`): source object.
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
        entity["type"]="protein"
        entity["name"]=data["Gene description"]
        entity["synonyms"]=data["Gene synonym"]
        entity["identifiers"]=[]
        if data["Uniprot"]==[]:
            entity["identifiers"].append({"namespace": "uniprot_id",
                "value":"None"})
        else:
            for uniprot_id in data["Uniprot"]:
                entity["identifiers"].append({"namespace": "uniprot_id",
                    "value":uniprot_id})
        entity["identifiers"].append({"namespace":"ensembl",
                    "value":data["Ensembl"]})
        entity["identifiers"].append({"namespace":"gene_name",
                    "value":data["Gene"]})
        entity["related"]=[]
        
        arr = ["protein_class","biological_process","molecular_function","disease_involvement"]
        for a in arr:
            if data[a.replace("_"," ").capitalize()]==None:
                entity["related"].append({"namespace":a,"value":None})
            else:
                for x in data[a.replace("_"," ").capitalize()]:
                    entity["related"].append({"namespace":a,"value":x})

        print(entity["name"])
        return entity
    
    def build_obs(self, data):
        """Build observation objects from obj.
        Go into observations collection.
        Args:
            data (:obj:`Obj`): source object.
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
        entity["name"]=data["Gene description"]
        entity["identifiers"]=[]
        entity["identifiers"].append({"namespace":"gene_name",
                    "value":data["Gene"]})
        
        values_p=[]
        values_p.append({"type":"protein_evidence","value":data["Evidence"]})
        protein_evidence=["HPA","UniProt","NeXtProt","MS"]
        for e in protein_evidence:
            values_p.append({"type":e+"_evidence","value":data[e+" evidence"]})
        
        
        rna_data=["tissue","cell_line","brain_regional","blood_cell","blood_lineage","cancer","mouse_brain_regional","pig_brain_regional"]
        for rna in rna_data:
            values_p.append({"type":"RNA_"+rna+"_specificity","value":data["RNA "+rna.replace("_"," ")+" specificity"]})
            values_p.append({"type":"RNA_"+rna+"_distribution","value":data["RNA "+rna.replace("_"," ")+" distribution"]})
            values_p.append({"type":"RNA_"+rna+"_specificity_score","value":data["RNA "+rna.replace("_"," ")+" specificity score"]})
        
        for i in range(5):
            if data["RNA "+rna_data[i].replace("_"," ")+" specific NX"]!=None:
                for x in data["RNA "+rna_data[i].replace("_"," ")+" specific NX"].keys():
                    values_p.append({"type":"RNA_"+rna_data[i]+"_specific_NX","value":data["RNA "+rna_data[i].replace("_"," ")+" specific NX"][x],"description":x.replace("\\","")})
            else:
                values_p.append({"type":"RNA_"+rna_data[i]+"_specific_NX","value":None})
        
        if data["RNA cancer specific FPKM"]!=None:
            for x in data["RNA cancer specific FPKM"]:
                values_p.append({"type":"RNA_cancer_specific_FPKM","value":data["RNA cancer specific FPKM"][x],"description":x})
        else:
            values_p.append({"type":"RNA_cancer_specific_FPKM","value":None})
        brain=["mouse","pig"]
        for b in brain:
            if data["RNA "+b+" brain regional specific pTPM"]!=None:
                for x in data["RNA "+b+" brain regional specific pTPM"]:
                    values_p.append({"type":"RNA_"+b+"_brain_regional_specific_pTPM","value":data["RNA "+b+" brain regional specific pTPM"][x],"description":x})
            else:
                values_p.append({"type":"RNA_"+b+"_brain_regional_specific_pTPM","value":None})
                
        
        for antibody in data["Antibody RRID"].keys():
            values_p.append({"type":"antibody","value":antibody,"description":"Human Protein Atlas Number"})
            values_p.append({"type":"antibody","value":data["Antibody RRID"][antibody],"description":"RRID"})
          
        values_p.append({"type":"reliability (IH)","value":data["Reliability (IH)"]})
        values_p.append({"type":"reliability (Mouse Brain)","value":data["Reliability (Mouse Brain)"]})
        values_p.append({"type":"reliability (IF)","value":data["Reliability (IF)"]})
        values_p.append({"type":"secretome location","value":data["Secretome location"]})
        values_p.append({"type":"blood_concentration - Conc._blood_IM [pg/L]","value":data["Blood concentration - Conc. blood IM [pg/L]"]})
        values_p.append({"type":"blood_concentration - Conc._blood_MS [pg/L]","value":data["Blood concentration - Conc. blood MS [pg/L]"]})
        
        if data["Subcellular main location"]!=None:
            for location in data["Subcellular main location"]:
                values_p.append({"type":"subcellular location","value":location,"description":"main"})
        if data["Subcellular additional location"]!=None:
            for location in data["Subcellular additional location"]:
                values_p.append({"type":"subcellular location","value":location,"description":"additional"})
        if (data["Subcellular main location"]==None) and (data["Subcellular additional location"]==None):
            values_p.append({"type":"subcellular location","value":None})
            
        cancers=["breast","cervical","colorectal","endometrial","head and neck","liver","lung","ovarian","pancreatic","prostate","renal","stomach","testis","thyroid","urothelial"]
        disease=["glioma","melanoma"]
        for cancer in cancers:
            if data["Pathology prognostics - "+cancer.capitalize()+" cancer"]==None:
                values_p.append({"type":"prognostic type","value":None,"description":cancer+" cancer"})
                values_p.append({"type":"is_prognostic","value":None,"description":cancer+" cancer"})
                values_p.append({"type":"p_val","value":None,"description":cancer+" cancer"})
            else:
                values_p.append({"type":"prognostic type","value":data["Pathology prognostics - "+cancer.capitalize()+" cancer"]["prognostic type"],"description":cancer+" cancer"})
                values_p.append({"type":"is_prognostic","value":str(data["Pathology prognostics - "+cancer.capitalize()+" cancer"]["is_prognostic"]),"description":cancer+" cancer"})
                values_p.append({"type":"p_val","value":data["Pathology prognostics - "+cancer.capitalize()+" cancer"]["p_val"],"description":cancer+" cancer"})
        
        for d in disease:
            if data["Pathology prognostics - "+d.capitalize()]==None:
                values_p.append({"type":"prognostic type","value":None,"description":d})
                values_p.append({"type":"is_prognostic","value":None,"description":d})
                values_p.append({"type":"p_val","value":None,"description":d})
            else:
                values_p.append({"type":"prognostic type","value":data["Pathology prognostics - "+d.capitalize()]["prognostic type"],"description":d})
                values_p.append({"type":"is_prognostic","value":str(data["Pathology prognostics - "+d.capitalize()]["is_prognostic"]),"description":d})
                values_p.append({"type":"p_val","value":data["Pathology prognostics - "+d.capitalize()]["p_val"],"description":d})
         
        tissue_rna = ["adipose tissue","adrenal gland","amygdala","appendix","basal ganglia","bone marrow","breast","cerebellum","cerebral cortex","cervix, uterine","colon","corpus callosum","ductus deferens","duodenum","endometrium 1","epididymis","esophagus","fallopian tube","gallbladder","heart muscle","hippocampal formation","hypothalamus","kidney","liver","lung","lymph node","midbrain","olfactory region","ovary","pancreas","parathyroid gland","pituitary gland","placenta","pons and medulla","prostate","rectum","retina","salivary gland","seminal vesicle","skeletal muscle","skin 1","small intestine","smooth muscle","spinal cord","spleen","stomach 1","testis","thalamus","thymus","thyroid gland","tongue","tonsil","urinary bladder","vagina","B-cells","dendritic cells","granulocytes","monocytes","NK-cells","T-cells","total PBMC"]
        for tissue in tissue_rna:
            values_p.append({"type":"tissue_RNA_NX","value":data["Tissue RNA - "+tissue+" [NX]"],"description":tissue})
        
        cell_RNA = ["A-431","A549","AF22","AN3-CA","ASC diff","ASC TERT1","BEWO","BJ","BJ hTERT+","BJ hTERT+ SV40 Large T+","BJ hTERT+ SV40 Large T+ RasG12V","CACO-2","CAPAN-2","Daudi","EFO-21","fHDF/TERT166","HaCaT","HAP1","HBEC3-KT","HBF TERT88","HDLM-2","HEK 293","HEL","HeLa","Hep G2","HHSteC","HL-60","HMC-1","HSkMC","hTCEpi","hTEC/SVTERT24-B","hTERT-HME1","HUVEC TERT2","K-562","Karpas-707","LHCN-M2","MCF7","MOLT-4","NB-4","NTERA-2","PC-3","REH","RH-30","RPMI-8226","RPTEC TERT1","RT4","SCLC-21H","SH-SY5Y","SiHa","SK-BR-3","SK-MEL-30","T-47d","THP-1","TIME","U-138 MG","U-2 OS","U-2197","U-251 MG","U-266/70","U-266/84","U-698","U-87 MG","U-937","WM-115"]
        for cell in cell_RNA:
            values_p.append({"type":"cell_RNA_NX","value":data["Cell RNA - "+cell.replace("\\","")+" [NX]"],"description":cell})
        
        blood_RNA = ["basophil","classical monocyte","eosinophil","gdT-cell","intermediate monocyte","MAIT T-cell","memory B-cell","memory CD4 T-cell","memory CD8 T-cell","myeloid DC","naive B-cell","naive CD4 T-cell","naive CD8 T-cell","neutrophil","NK-cell","non-classical monocyte","plasmacytoid DC","T-reg","total PBMC"]
        for b in blood_RNA:
            values_p.append({"type":"blood_RNA_NX","value":data["Blood RNA - "+b+" [NX]"],"description":b})
        
        brain_RNA = ["amygdala","basal ganglia","cerebellum","cerebral cortex","hippocampal formation","hypothalamus","midbrain","olfactory region","pons and medulla","thalamus"]
        for brain in brain_RNA:
            values_p.append({"type":"brain_RNA_NX","value":data["Brain RNA - "+brain+" [NX]"],"description":brain})
        
        genotype = {"taxon": {"ncbi_taxonomy_id":int(9606),
                              "name":"Homo sapiens (Human)"}}
        genotype["taxon"]["canon_ancestors"]=[]
        query={"tax_id":int(9606)}
        projection={"_id":0,"canon_anc_ids":1,"canon_anc_names":1}
        doc = self.client["datanator-test"]["taxon_tree"].find_one(filter=query, projection=projection)
        if doc!=None:
            for j in range(len(doc["canon_anc_names"])):
                d={}
                d["ncbi_taxonomy_id"]=doc["canon_anc_ids"][j]
                d["name"]=doc["canon_anc_names"][j]
                genotype["taxon"]["canon_ancestors"].append(d)
        genotype["chromosome"]=data["Chromosome"]
        genotype["position"]=data["Position"]
        source=[{"namespace":"doi","value":"10.1126/science.aal3321"}]
        
        ob_p ={"entity":entity,
               "genotype":genotype,
               "values":values_p,
               "source":source,
               "schema_version":"2.0"}
        
        return ob_p
        
                    
    def process_docs(self):
            #read .tsv file
            data=pd.read_csv('subcellular_location.tsv',delimiter="\t")
            data = data.where(pd.notnull(data), None)
            proteins =[]
            for i in range(3999,7999):
                if str(data.iloc[i,0]) not in proteins:
                    proteins.append(str(data.iloc[i,0]))

            for x in proteins:
                url = urllib.request.urlopen("https://www.proteinatlas.org/"+x+".json")
                data = json.loads(url.read().decode())
                
                #update entity collection
                entity = self.build_entity(data)
                query={"identifiers":{"$elemMatch":{"namespace":"gene_name","value":str(data["Gene"])}}}
                self.entity_col.update_one(query,
                                               {"$set": {
                                                     "type": entity["type"],
                                                     "name": entity["name"],                      
                                                     "schema_version": "2.0"},
                                            "$addToSet": {
                                                "synonyms": {"$each":entity["synonyms"]},
                                                "identifiers": {"$each": entity["identifiers"]},
                                                "related": {"$each":entity["related"]}},
                                            "$unset": {
                                                "protein_class": "",
                                                "biological_process": "",
                                                "molecular_function": "",
                                                "disease_involvement": ""}},
                                           upsert=True)
            
                obs = self.build_obs(data)
                #update identifier collection
                query={"$and":[{"namespace":"gene_name"},{"value":str(data["Gene"])}]}
                self.identifier_col.update_one(query,
                                                {"$set": {"namespace": "gene_name",
                                                  "value": str(data["Gene"]),
                                                  "description":str(entity["name"])}},upsert=True)

                #update observation collection
                con_1={"identifier":{"namespace":"gene_name","value":str(data["Gene"])}}
                con_2={"source":{"$elemMatch":{"namespace":"doi","value":"10.1126/science.aal3321"}}}
                query={"$and":[con_1,con_2]}
                self.obs_col.update_one(query,
                                        {"$set": {"entity": obs["entity"],
                                                  "genotype": obs["genotype"],
                                                  "schema_version": "2.0",
                                                  "identifier":{"namespace":"gene_name","value":str(data["Gene"])}},
                                         "$addToSet": {"values": {"$each": obs["values"]},
                                                       "source": {"$each": obs["source"]}}},
                                         upsert=True)
                
                

def main():
    conf=config.Victoria()
    conf_main = config.Config()
    username = conf.USERNAME
    password = conf.PASSWORD
    MongoDB = conf_main.SERVER
    src = insert_protein_atlas(MongoDB = MongoDB,
                                       username=username,
                                       password=password,
                                       collection = "observation",
                                       db = "datanator-demo")
    src.process_docs()
if __name__== '__main__':
    main()
