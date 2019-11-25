import pandas as pd
from datanator.util import mongo_util
import requests
import io


class Halflife(mongo_util.MongoUtil):

    def __init__(self, cache_dir=None, server=None, db=None, collection_str=None,
                authDB=None, readPreference=None, username=None, password=None,
                verbose=None, max_entries=None):
        super(Halflife, self).__init__(MongoDB=server, db=db, username=username,
                                 password=password, authSource=authDB,
                                 verbose=verbose)
        self.cache_dir = cache_dir
        self.client, self.db, self.collection = self.con_db(collection_str)
        self.url = "https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-016-3219-8/MediaObjects/12864_2016_3219_MOESM5_ESM.xlsx"
        self.max_entries = max_entries

    def download_xlsx(self):
        response = requests.get(self.url)
        if self.max_entries == float('inf'):
            nrows = None
        else:
            nrows = self.max_entries
        data = pd.read_csv(io.BytesIO(response.content),
                               delimiter='\t', encoding='utf-8', nrows=nrows)
        return data