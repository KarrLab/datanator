from datanator_query_python.util import mongo_util
import pandas as pd
import json
import pymongo
import io


class BrendaRxn(mongo_util.MongoUtil):

    def __init__(self, MongoDB=None, db=None, username=None, password=None,
                verbose=False, max_entries=float('inf'), authSource='admin', readPreference='nearest',
                collection_str='brenda_reaction'):
        super().__init__(MongoDB=MongoDB, db=db, username=username, password=password,
                        verbose=verbose, max_entries=max_entries, authSource=authSource,
                        readPreference=readPreference)
        self.url = 'http://bkms-react.tu-bs.de/download/Reactions_BKMS.tar.gz'
        self.collection = self.db_obj[collection_str]
        self.max_entries = max_entries

    def download_and_read(self):
        # dtype = {'EC_Number': str, 'Enzyme_Name': str, 'Reaction': str,
        #         'BRENDA_Pathway_Name': str, 'KEGG_Pathway_Name': str,
        #         'MetaCyc_Pathway_Name': str, 'Reaction_ID_BRENDA': str,
        #         'Reaction_ID_KEGG': str, 'Reaction_ID_MetaCyc': str,
        #         'Reaction_ID_SABIO_RK': int, 'stoichiometry': str,
        #         'Commentary_KEGG': str, 'Commentary_MetaCyc': str}
        doc = './docs/brenda/Reactions_BKMS.tar.gz'
        df = pd.read_csv(doc, compression='gzip', header=0, error_bad_lines=False,
                         engine='c', sep='\t',
                         lineterminator='\n', low_memory=False)
        df.columns = [x.lower() for x in ['Count', 'EC_Number', 'enzyme_name', 
                                          'Reaction', 'Reaction_ID_BRENDA', 'Reaction_ID_KEGG', 
                                          'Reaction_ID_MetaCyc', 'Reaction_ID_SABIO_RK', 'BRENDA_Pathway_Name', 
                                          'KEGG_Pathway_ID', 'KEGG_Pathway_Name', 'MetaCyc_Pathway_ID', 
                                          'MetaCyc_Pathway_Name', 'Stoichiometry', 'Missing_Substrate', 'Missing_Product', 
                                          'Commentary_KEGG', 'Commentary_MetaCyc', 'Remark']]
        return df

    def clean_up(self, df):
        """Clean up column values in dataframe, i.e. split strings
        into arrays.
        
        Args:
            df (:obj:`pandas.DataFrame`): dataframe to be processed.

        Return:
            (:obj:`pandas.DataFrame`): dataframe after being processed.
        """
        df['brenda_pathway_name'] = df['brenda_pathway_name'].str.split('; ')
        df['kegg_pathway_name'] = df['kegg_pathway_name'].str.split('; ')
        df['kegg_pathway_id'] = df['kegg_pathway_id'].str.split('; ')
        df['metacyc_pathway_id'] = df['metacyc_pathway_id'].str.split('; ')
        df['metacyc_pathway_name'] = df['metacyc_pathway_name'].str.split('; ')
        df['commentary_kegg'] = df['commentary_kegg'].str.split('; ')
        df['reaction_id_brenda'] = df['reaction_id_brenda'].str.split(',')
        df['reaction_id_kegg'] = df['reaction_id_kegg'].str.split(',')
        df['reaction_id_metacyc'] = df['reaction_id_metacyc'].str.split(',')
        df['reaction_id_sabio_rk'] = pd.to_numeric(df['reaction_id_sabio_rk'], errors='coerce')
        return df

    def parse_reaction(self, df):
        """Parse 'Reaction' column into substrates and products
        e.g. 'ATP + (R)-pantoate + beta-alanine = AMP + diphosphate + (R)-pantothenate',
        'ATP + Detyrosinated alpha-tubulin + L-Tyrosine <=> alpha-Tubulin + ADP + Orthophosphate'

        Args:
            df(:obj:`pandas.DataFrame`): DataFrame containing reaction equation column.

        Return:
            (:obj:`pandas.DataFrame`): DataFrame with parsed substrates and products columns.
        """
        df['reaction'] = df['reaction'].str.replace(' <=> ', ' = ')
        df['reactants'] = df['reaction'].str.split(' = ')
        df['substrates'] = df['reactants'].apply(lambda x: str(x[0]).split(' + ') if (isinstance(x, list)) else [x])
        df['products'] = df['reactants'].apply(lambda x: str(x[1]).split(' + ') if (isinstance(x, list)) else [x])
        result = df.drop(labels=['reaction', 'reactants'], axis=1)
        return result

    def load_df(self, df):
        df_json = json.loads(df.to_json(orient='records'))
        try:        
            self.collection.insert_many(df_json)
        except pymongo.errors.InvalidOperation as e:
            return(str(e))


import datanator.config.core

def main():
    db = 'datanator'
    collection_str = 'brenda_reactions'
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    MongoDB = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']
    manager = BrendaRxn(MongoDB=MongoDB, db=db, username=username,
                        password=password, collection_str=collection_str)
    df = manager.download_and_read()
    new_df = manager.clean_up(df)
    result = manager.parse_reaction(new_df)
    manager.load_df(result)


if __name__ == '__main__':
    main()
