import pandas as pd
import json
import numpy as np

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        """
        Converts the dictionary's values into a JSON serializable data type
        
        """
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NpEncoder, self).default(obj)
        
class ParsePsortExperimental:
    def __init__(self, max_entries):
        self.max_entries = max_entries

    def parse_psortdb(self):
        """
        To parse database psortdb Experimental-PSORTdb-v4.00.tsv file
        and create JSON files conforming to datanator_pattern/observation_compiled.json

        Args:
            max_entries: int
                number of rows to parse.
                A JSON file will be created for each of the tsv file's first <max_entries> rows

        Return:
            ()
        """
        data=pd.read_csv('Experimental-PSORTdb-v4.00.tsv',delimiter="\t")
        data = data.where(pd.notnull(data), None)
        for i in range(self.max_entries):
            d={}
            #entity
            d["entity"]={}
            d["entity"]["type"]="protein"
            d["entity"]["name"]=str(data.iloc[i,6]).replace(".","")
            if data.iloc[i,7] != None:
                    d["entity"]["synonyms"]=str(data.iloc[i,7]).split(",")
            else:
                d["entity"]["synonyms"]=[]
            #identifiers
            d["entity"]["identifiers"]=[]
            uniprot={}
            uniprot["name_space"]="uniprot_id"
            uniprot["value"]=data.iloc[i,0]
            ref_seq = {}
            ref_seq["name_space"]="Refseq_Accession"
            ref_seq["value"]=data.iloc[i,1]
            other_accession = {}
            other_accession["name_space"]="Other_Accession"
            other_accession["value"]=data.iloc[i,2]
            d["entity"]["identifiers"].append(uniprot)
            d["entity"]["identifiers"].append(ref_seq)
            d["entity"]["identifiers"].append(other_accession)

            #localizations
            d["value"]={}
            if data.iloc[i,3] != None:
                d["value"]["experimental_localization"] = str(data.iloc[i,3]).split(",")
            else:
                d["value"]["experimental_localization"] = []
            if data.iloc[i,4] != None:
                d["value"]["secondary_localizaton"] = str(data.iloc[i,4]).split(",")
            else:
                d["value"]["secondary_localizaton"] = []

            #genotype
            d["genotype"]={}
            d["genotype"]["taxon"]={}
            d["genotype"]["taxon"]["ncbi_taxonomy_id"]=data.iloc[i,9]
            d["genotype"]["taxon"]["name"]=data.iloc[i,10]

            #environment
            d["environment"]={}
            d["environment"]["GramStain"]=data.iloc[i,13]
            
            #source
            d["source"]={}
            d["source"]["namespace"]="ePSORTdb"
            d["source"]["value"]="Version "+str(data.iloc[i,17])

            #name is the JSON file's name
            if (data.iloc[i,0]!=None):
                name = data.iloc[i,0]   #SwissProt_ID
            else:
                name = data.iloc[i,2]   #Other_Accession
            with open("Experimental_PSortdb/"+name+".json","w+") as f:
                json.dump(d,f,cls=NpEncoder,indent=4)

p1=ParsePsortExperimental(10)
p1.parse_psortdb()
