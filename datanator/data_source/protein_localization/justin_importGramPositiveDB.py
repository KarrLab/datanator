from datanator_query_python.config import config as q_conf
from datanator_query_python.util import mongo_util
import os
import json


class ImportMetabolite(mongo_util.MongoUtil):

    def __init__(self, 
                 MongoDB, 
                 db,
                 username,
                 password,
                 collection,
                 directory):
        super().__init__(MongoDB=MongoDB,
                         db=db,
                         username=username,
                         password=password)
        self.collection = self.db_obj[collection]
        self.directory = directory


    def update_collection(self):
        for filename in os.listdir(self.directory):
            print(filename)
            if filename.endswith(".json"):
                with open(self.directory+'/'+filename, "r") as f:
                    post = json.load(f)
                    self.collection.insert_one(post) # SHOULD BE update_one NOT insert_one

                    
def main():
    conf = q_conf.Justin()
    password = conf.PASSWORD
    username = conf.USERNAME
    MongoDB = conf.SERVER
    src = ImportMetabolite(MongoDB, 
                           "datanator-demo", 
                          username,
                          password,
                          "observation",
                          "./datanator/docs/protein_localization/computed_gram_positive/JSONSchema")
    src.update_collection() 


if __name__ == "__main__":
    main()