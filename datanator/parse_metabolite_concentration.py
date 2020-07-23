import pandas as pd
import urllib.request
import json
from datanator_query_python.config import config
from datanator_query_python.util import mongo_util

class ParseMetaboliteConcentration(mongo_util.MongoUtil):
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

    def parse_metabolite(self):
        """
        Read JSON metabolite concentration files from Github and
        insert separate documents for each metabolite into MongoDB database  

        Args:
            ()
        Return:
            ()
        """
        collection = self.db_obj[self.collection]
        metabolites = ["ATP","CTP","GMP","GTP","IMP","NAD","NADH","NADP","NADPH","TTP","UTP"]
        for i in range(len(metabolites)):       
            url = urllib.request.urlopen("https://raw.githubusercontent.com/KarrLab/datanator/tutorial/docs/metabolites/"+metabolites[i]+".json")
            data = json.loads(url.read().decode())
            collection.insert_one({"inchikey":data['inchikey']})
            for j in range(len(data['concentrations'])):
                sub_data = data['concentrations'][j]
                collection.update_one({"inchikey":data['inchikey']},{"$addToSet":{'concentrations':sub_data}})
def main():
    conf=config.Victoria()
    conf_main = config.Config()
    username = conf.USERNAME
    password = conf.PASSWORD
    MongoDB = conf_main.SERVER
    src = ParseMetaboliteConcentration(MongoDB = MongoDB,
                                       username=username,
                                       password=password,
                                       collection = "metabolite_concentration",
                                       db = "datanator-demo")
    src.parse_metabolite()

if __name__== '__main__':
    main()
        
