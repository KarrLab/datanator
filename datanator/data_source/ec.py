from datanator_query_python.util import mongo_util
import requests


class EC(mongo_util.MongoUtil):

    def __init__(self, server=None, db=None, username=None, password=None, 
                 authSource='admin', readPreference='nearest', collection_str='ec',
                 verbose=True, max_entries=float('inf')):
        super().__init__(MongoDB=server, db=db, verbose=verbose, max_entries=max_entries,
                        username=username, password=password, authSource=authSource,
                        readPreference=readPreference)
        self.max_entries = max_entries
        self.collection_str = collection_str
        self.client, self.db, self.collection = self.con_db(collection_str)

    def download_dat(self):
        """Download enzyme.dat file for parsing.
        (ftp://ftp.expasy.org/databases/enzyme/enzyme.dat)
        """
        url = "ftp://ftp.expasy.org/databases/enzyme/enzyme.dat"

