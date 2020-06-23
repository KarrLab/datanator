import pandas as pd
import json
import numpy as np
from datanator_query_python.config import config
from datanator_query_python.util import mongo_util


class ParsePsortExperimental(mongo_util.MongoUtil):

    def __init__(self, max_entries=float('inf'),
                       MongoDB=None,
                       db=None,
                       collection=None,
                       username=None,
                       password=None,
                       authSource='admin',
                       readPreference='nearest'):
        super().__init__(MongoDB=MongoDB, db=db,
                        username=username,
                        password=password,
                        authSource=authSource,
                        readPreference=readPreference)
        self.max_entries = max_entries
        self.collection = collection
​
    def parse_psortdb(self):
        """
        To parse database psortdb Experimental-PSORTdb-v4.00.tsv file
​
        Args:
            max_entries: int
                number of rows to parse.
                A JSON file will be created for each of the tsv file's first <max_entries> rows
​
        Return:
            ()
        """
        collection = self.db_obj[self.collection]
        # data=pd.read_csv('Experimental-PSORTdb-v4.00.tsv',delimiter="\t")
        # data = data.fillna("None")
        # header = list(data.columns.values)
        # for i in range(self.max_entries):
        #     d={}
        #     for j in range(len(header)):
        #         if isinstance(data.iloc[i,j],int):
        #             data.iloc[i,j]=int(data.iloc[i,j])
        #         else:
        #             d[header[j]] = data.iloc[i,j]
        #     #name is the JSON file's name
        #     if (data.iloc[i,0]!="None"):
        #         name = data.iloc[i,0]   #SwissProt_ID
        #     else:
        #         name = data.iloc[i,2]   #Other_Accession
        #     with open(name+".json","w+") as f:
        #         x = json.dumps(d,f,cls=NpEncoder)
        collection.update_one({"uniprot_id": "P01234"},
                                {"$set": {"protein_name": "some_name",
                                        "another_field": "another_value"},
                                "$addToSet": {"$each": {"add_id": [{"namespace": "something",
                                                            "value": "1"}]}}})
            
            
​
def main():
    conf = config.Justin()
    conf_main = config.Config()
    username = conf.USERNAME
    password = conf.PASSWORD
    MongoDB = conf_main.SERVER
    src = ParsePsortExperimental(MongoDB=MongoDB,
                                 username=username,
                                 password=password,
                                 collection="protein_localization",
                                 db="datanator-demo")
    src.parse_psortdb()    


if __name__ == '__main__':
    main()