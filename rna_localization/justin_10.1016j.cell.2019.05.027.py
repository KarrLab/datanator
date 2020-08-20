import pandas as pd 
from datanator_query_python.util import mongo_util
from datanator_query_python.config import config 


class RNALocate(mongo_util.MongoUtil):
    def __init__(self, MongoDB, db, username, password):
        super().__init__(MongoDB=MongoDB,
                         db=db,
                         username=username,
                         password=password)
        self.identifier_collection = self.db_obj['identifier']
        self.observation_collection = self.db_obj['observation']

    def parse_rna_location(self):
        data = pd.read_excel('mmc3.xlsx', 'Gene lists and orphans')
        for i in range(len(data)):
            d = {}
            d['entity'] = {'type': 'RNA',
                           'name': data['Common_Gene'][i],
                           'identifiers': [{'namespace': 'ensembl',
                                            'value': data['Ensembl_Gene'][i]}]}
            d['identifier'] = {'namespace': 'ensembl',
                               'value': data['Ensembl_Gene'][i],
                               'description': data['Common_Gene'][i]}
            d['values'] = [{'type': col_name, 'value': data[col_name][i]} for col_name in data.columns[3:11]]
            d['source'] = [{'namespace': 'doi', 'value': '10.1016/j.cell.2019.05.027'}]
            d['schema_version'] = '2.0'

            self.observation_collection.update_one({'type': 'RNA',
                                                   'name': data['Common_Gene'][i],
                                                   'identifiers': [{'namespace': 'ensembl',
                                                                     'value': data['Ensembl_Gene'][i]}]},
                                                    {'$set': d},
                                                    upsert=True)

            self.identifier_collection.update_one({'namespace': 'ensembl', 'value': data['Ensembl_Gene'][i]},
                                                  {'$set': {'description': data['Common_Gene'][i]}},
                                                  upsert=True)
            
            print("Row {} has been added".format(str(i)))



def main():
    conf = config.Justin()
    username = conf.USERNAME
    password = conf.PASSWORD
    MongoDB = conf.SERVER
    db = 'datanator-demo'
    #url = 'https://www.cell.com/cms/10.1016/j.cell.2019.05.027/attachment/618723b6-c0fb-4138-846e-fb09eb6b2f2f/mmc3'
    src = RNALocate(MongoDB, db, username, password)
    src.parse_rna_location()

if __name__ == "__main__":
    main()