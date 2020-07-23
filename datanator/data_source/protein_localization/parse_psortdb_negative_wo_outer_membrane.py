import pandas as pd
import json

class ParsePsort:
    def __init__(self, max_entries):
        self.max_entries = max_entries
        
    def parse_psortdb(self):
        """
        To parse database psortdb gram negative without outer membrane file
        and create JSON files conforming to datanator_pattern/observation_compiled.json
        
        Args:
            max_entries(:obj:'int'): number of rows to parse.
            A JSON file will be created for each of the first <max_entries> rows

        Return:
            ()
        """
        data=pd.read_csv('Computed-Gram_negative_without_outer_membrane-PSORTdb-3.00.tab',delimiter="\t",low_memory=False)
        data = data.where(pd.notnull(data), None)
        for i in range(self.max_entries):
            d={}
            #entity
            d["entity"]={}
            d["entity"]["type"]="protein"
            d["entity"]["name"]=str(data.iloc[i,0])[str(data.iloc[i,0]).rfind("|")+2:]
            d["entity"]["synonyms"]=[]
            #identifiers
            d["entity"]["identifiers"]=[]
            seq_id = {}
            seq_id["namespace"]="Seq_ID"
            seq_id["value"]=str(data.iloc[i,0])[str(data.iloc[i,0]).find("ref")+4:str(data.iloc[i,0]).rfind("|")]
            d["entity"]["identifiers"].append(seq_id)
            #localizations
            d["value"]={}
            d["value"]["PPSVM_Localization"]=data.iloc[i,1]
            d["value"]["Profile_Localization"]=data.iloc[i,3]
            d["value"]["Signal_Localization"]=data.iloc[i,5]
            d["value"]["SCL-BLASTe_Localization"]=data.iloc[i,7]
            d["value"]["CMSVM_Localization"]=data.iloc[i,9]
            d["value"]["SCL-BLAST_Localization"]=data.iloc[i,11]
            d["value"]["OMPMotif_Localization"]=data.iloc[i,13]
            d["value"]["OMSVM_Localization"]=data.iloc[i,15]
            d["value"]["Motif_Localization"]=data.iloc[i,17]
            d["value"]["CytoSVM_Localization"]=data.iloc[i,19]
            d["value"]["CWSVM_Localization"]=data.iloc[i,21]
            d["value"]["ModHMM_Localization"]=data.iloc[i,23]
            d["value"]["ECSVM_Localization"]=data.iloc[i,25]
            d["value"]["Cytoplasmic Membrane_Score"]=data.iloc[i,27]
            d["value"]["Cellwall_Score"]=data.iloc[i,28]
            d["value"]["Extracellular_Score"]=data.iloc[i,29]
            d["value"]["Cytoplasmic_Score"]=data.iloc[i,30]                    
            d["value"]["Final_Localization"]=data.iloc[i,31]
            d["value"]["Final_Localization_2"]=data.iloc[i,32]
            d["value"]["Secondary_Localization"]=data.iloc[i,34]
            d["value"]["Final_Score"]=data.iloc[i,35]

            #source
            d["source"]={}
            d["source"]["namespace"]="PSORT"
            d["source"]["value"]="Version "+str(data.iloc[i,36])
            with open("Gram_Negative_WO_Outer_Membrane/"+str(data.iloc[i,0])[str(data.iloc[i,0]).find("ref")+4:str(data.iloc[i,0]).rfind("|")]+".json","w+") as f:
                json.dump(d,f,indent=4)

p1=ParsePsort(10)
p1.parse_psortdb()
