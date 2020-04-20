from datanator_query_python.util import mongo_util
import pandas as pd


class BrendaRxn(mongo_util.MongoUtil):

    def __init__(self, MongoDB=None, db=None, username=None, password=None,
                verbose=False, max_entries=float('inf'), authSource='admin', readPreference='nearest',
                collection_str='brenda_reaction'):
        super().__init__(MongoDB=MongoDB, db=db, username=username, password=password,
                        verbose=verbose, max_entries=max_entries, authSource=authSource,
                        readPreference=readPreference)
        self.url = 'http://bkms-react.tu-bs.de/download/Reactions_BKMS.tar.gz'
        self.collection = self.db_obj[collection_str]

    def download_and_read(self):
        dtype = {'EC_Number': str, 'Enzyme_Name': str, 'Reaction': str,
                'BRENDA_Pathway_Name': str, 'KEGG_Pathway_Name': str,
                'MetaCyc_Pathway_Name': str, 'Reaction_ID_BRENDA': str,
                'Reaction_ID_KEGG': str, 'Reaction_ID_MetaCyc': str,
                'Reaction_ID_SABIO_RK': int, 'stoichiometry': str,
                'Commentary_KEGG': str, 'Commentary_MetaCyc': str}
        df = pd.read_csv(self.url, compression='gzip', header=0, error_bad_lines=False,
                        dtype=dtype)
        return df