import json
import pandas as pd  
import os

class ParseJSONSchema:
    
    def __init__(self, dataset, directory):
        self.dataset = dataset
        self.directory = directory


    def update_directory(self):
        data = pd.read_csv(self.dataset, delimiter='\t', nrows=10000)
        data = data.where(pd.notnull(data), None)

        for i in range(len(data)):
            d = {}
            
            # data to "entity"
            d["entity"] = {}
            d["entity"]["type"] = "protein"
            seq_id = str(data.iloc[i,0])
            d["entity"]["name"] = seq_id[str(data.iloc[i,0]).rfind('|')+2:]
            
            # data to "identifier" in "entity"
            identifier = []
            dict_identifier = {} # dictionary for identifier in entity
            dict_identifier["namespace"] = "SeqID"
            dict_identifier["value"] = seq_id[8:22]
            identifier.append(dict_identifier)
            d["entity"]["identifiers"] = identifier

            # data to values
            d["values"] = []
            for column_name in data.columns[1:]:
                if column_name != "SeqID" and column_name != "PSortVersion":
                    values_dict = {}
                    values_dict["type"] = column_name
                    values_dict["value"] = data[column_name].iloc[i]
                    d["values"].append(values_dict)

            # data to "identifier"
            d["identifier"] = {}
            d["identifier"]["namespace"] = "SeqID"
            d["identifier"]["value"] = seq_id[8:22]


            # source
            d["source"] = []
            dict_source = {}
            dict_source["namespace"] = "PSORTsb Gram Positive"
            dict_source["value"] = "Version 3"
            d["source"].append(dict_source)

            # environment
            d["environment"] = {"GramStain": "Gram positive"}

            # Schema Version
            d["schema_version"] = "2.0"


            # Create JSON files and place in directory
            with open(self.directory+"/{}.json".format(seq_id[8:22]), "w+") as JSONfile:
                json.dump(d, JSONfile)


def main():
    json_files = ParseJSONSchema(dataset="./datanator/docs/protein_localization/computed_gram_positive/Computed-Gram_positive-PSORTdb-3.00.tab", 
                                directory="./datanator/docs/protein_localization/computed_gram_positive/JSONSchema")

    json_files.update_directory()

if __name__ == "__main__":
    main()