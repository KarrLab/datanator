import pandas as pd
from datanator.util import mongo_util
import json
import datetime
from pymongo import ASCENDING
from pymongo.collation import Collation, CollationStrength
import datanator.config.core


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
        self.collation = Collation(locale='en', strength=CollationStrength.SECONDARY)

    def download_xlsx(self, sheet_name):
        """Download supplementary xlsx file
        
        Args:
            sheet_name (:obj:`str`): name of sheet in xlsx
        
        Returns:
            (:obj:`pandas.DataFrame`): xlsx transformed to pandas.DataFrame
        """
        if self.max_entries == float('inf'):
            nrows = None
        else:
            nrows = self.max_entries
        data = pd.read_excel(self.url, sheet_name=sheet_name, nrows=nrows)
        columns = ['gene_fragment', 'cog_class', 'ar_cog', 'cog', 'function', 'gene_name', 'half_life', 'half_life_std', 'std_over_avg']
        data.columns = columns
        data['half_life'] = data['half_life'].apply(lambda x: x*60)
        data['half_life_std'] = data['half_life_std'].apply(lambda x: x*60)
        return data

    def load_halflife(self, df, growth_medium='MeOH'):
        df_json = json.loads(df.to_json(orient='records'))
        row_count = len(df.index)
        for i, doc in enumerate(df_json):
            if i == self.max_entries:
                break
            if self.verbose and i % 100 == 0:
                print('Processing {} row {} out of {}'.format(growth_medium, i, row_count))
            doc['halflives'] = [{'halflife': doc['half_life'], 'std': doc['half_life_std'], 'std_over_avg': doc['std_over_avg'],
                                'unit': 's', 'reference': [{'doi': '10.1186/s12864-016-3219-8'}], 'growth_medium': growth_medium,
                                'gene_position': doc['gene_fragment'], 'ar_cog': doc['ar_cog'], 'cog_class': doc['cog_class'],
                                'cog': doc['cog'], 'species': 'Methanosarcina acetivorans', 'ncbi_taxonomy_id': 2214}]
            doc['modified'] = datetime.datetime.utcnow()
            del doc['half_life']
            del doc['half_life_std']
            del doc['std_over_avg']
            del doc['ar_cog']
            del doc['cog']
            del doc['cog_class']

            if doc['gene_name'] != '-':
                self.collection.update_one({'gene_name': doc['gene_name']},
                                    {'$set': doc}, upsert=True, collation=self.collation)
            elif doc['function'] != '-':
                self.collection.update_one({'function': doc['function']},
                                           {'$set': doc}, upsert=True, collation=self.collation)
            else:
                self.collection.update_one({'halflives.gene_position': doc['gene_fragment']},
                                        {'$set': doc}, upsert=True, collation=self.collation)
            if i == 0:
                self.collection.create_index([("gene_name", ASCENDING)], background=True,
                                            collation=self.collation)
                self.collection.create_index([("halflives.gene_position", ASCENDING)], background=True,
                                            collation=self.collation)
                self.collection.create_index([("function", ASCENDING)], background=True,
                                            collation=self.collation)

        self.collection.update_many({'gene_fragment':{'$exists': True}}, {'$unset': {'gene_fragment': ""}})

    def add_to_halflife(self, df, growth_medium='TMA'):
        """Add df to existing rna_halflife collection
        
        Args:
            df (:obj:`pandas.DataFrame`): dataframe to be added.
            growth_medium (:obj:`str`): medium in which the cells were grown. Defaults to TMA.
        """
        df_json = json.loads(df.to_json(orient='records'))
        row_count = len(df.index)
        for i, doc in enumerate(df_json):
            if i == self.max_entries:
                break
            if self.verbose and i % 100 == 0:
                print('Processing {} row {} out of {}'.format(growth_medium, i, row_count))
            to_add = {'halflife': doc['half_life'], 'std': doc['half_life_std'], 'std_over_avg': doc['std_over_avg'],
                    'unit': 's', 'reference': [{'doi': '10.1186/s12864-016-3219-8'}], 'growth_medium': growth_medium,
                    'gene_position': doc['gene_fragment'], 'ar_cog': doc['ar_cog'], 'cog_class': doc['cog_class'],
                    'cog': doc['cog'], 'species': 'Methanosarcina acetivorans', 'ncbi_taxonomy_id': 2214}
            if doc['gene_name'] != '-':
                self.collection.update_one({'gene_name': doc['gene_name']},
                                        {'$addToSet': {'halflives': to_add},
                                        '$set': {'modified': datetime.datetime.utcnow()}}, 
                                        upsert=True, collation=self.collation)
            elif doc['function'] != '-':
                self.collection.update_one({'function': doc['function']},
                                        {'$addToSet': {'halflives': to_add},
                                        '$set': {'modified': datetime.datetime.utcnow()}}, 
                                        upsert=True, collation=self.collation)
            else:
                query = {'halflives.gene_position': doc['gene_fragment']}
                result = self.collection.find_one(filter=query, collation=self.collation)
                if result is not None:
                    self.collection.update_one(query,
                                            {'$addToSet': {'halflives': to_add},
                                            '$set': {'modified': datetime.datetime.utcnow()}}, 
                                            upsert=True, collation=self.collation)
                else:
                    doc['halflives'] = [to_add]
                    doc['modified'] = datetime.datetime.utcnow()
                    del doc['half_life']
                    del doc['half_life_std']
                    del doc['std_over_avg']
                    del doc['ar_cog']
                    del doc['cog']
                    del doc['cog_class']
                    self.collection.update_one(query, {'$set': doc}, upsert=True, collation=self.collation)
        self.collection.update_many({'gene_fragment':{'$exists': True}}, {'$unset': {'gene_fragment': ""}})        


def main():
    src_db = 'datanator'
    collection_str = 'rna_halflife'
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    server = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']       
    src = Halflife(username=username, password=password, server=server, 
                    authDB='admin', db=src_db,
                    verbose=True, collection_str=collection_str)
    df = src.download_xlsx('MeOH')
    src.load_halflife(df)
    df = src.download_xlsx('TMA')
    src.add_to_halflife(df, growth_medium='TMA')
    df = src.download_xlsx('Acetate')
    src.add_to_halflife(df, growth_medium='Acetate')

if __name__ == '__main__':
    main()