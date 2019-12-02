import pandas as pd
from datanator.util import mongo_util
import json
import datetime
from pymongo import ASCENDING
from pymongo.collation import Collation, CollationStrength
import datanator.config.core
from datanator_query_python.query import query_uniprot
from datanator.data_source import uniprot_nosql


class Halflife(mongo_util.MongoUtil):

    def __init__(self, cache_dir=None, server=None, db=None, collection_str=None,
                authDB=None, readPreference=None, username=None, password=None,
                verbose=None, max_entries=None, uniprot_col_db=None):
        super(Halflife, self).__init__(MongoDB=server, db=db, username=username,
                                 password=password, authSource=authDB,
                                 verbose=verbose)
        self.cache_dir = cache_dir
        self.client, self.db, self.collection = self.con_db(collection_str)
        self.url = "https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-016-3219-8/MediaObjects/12864_2016_3219_MOESM5_ESM.xlsx"
        self.max_entries = max_entries
        self.collation = Collation(locale='en', strength=CollationStrength.SECONDARY)
        self.uniprot_manager = query_uniprot.QueryUniprot(username=username, password=password,
                                            server=server, database='datanator', collection_str='uniprot')
        self.uniprot_col_manager = uniprot_nosql.UniprotNoSQL(MongoDB=server, db=uniprot_col_db, max_entries=max_entries,
                                                              username=username, password=password, collection_str='uniprot')

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
                                'ordered_locus_name': doc['gene_fragment'], 'ar_cog': doc['ar_cog'], 'cog_class': doc['cog_class'],
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
                self.collection.update_one({'halflives.ordered_locus_name': doc['gene_fragment']},
                                        {'$set': doc}, upsert=True, collation=self.collation)
            if i == 0:
                self.collection.create_index([("gene_name", ASCENDING)], background=True,
                                            collation=self.collation)
                self.collection.create_index([("halflives.ordered_locus_name", ASCENDING)], background=True,
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
                    'ordered_locus_name': doc['gene_fragment'], 'ar_cog': doc['ar_cog'], 'cog_class': doc['cog_class'],
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
                query = {'halflives.ordered_locus_name': doc['gene_fragment']}
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

    def fill_protein_name(self):
        """Create and fill 'protein_name' field for documents
        'gene_name' field with values other than '-'
        """
        con_0 = {'gene_name': {'$ne': '-'}}
        con_1 = {'gene_name': {'$exists': True}}
        query = {'$and': [con_0, con_1]}
        projection = {'gene_name': 1}
        docs = self.collection.find(filter=query, projection=projection, collation=self.collation)
        for doc in docs:
            gene_name = doc['gene_name']
            q = {'gene_name': gene_name}
            result = self.uniprot_manager.collection.find_one(filter=q, collation=self.collation, projection={'protein_name': 1})
            if result is not None:
                protein_name = result['protein_name']
            else:
                protein_name = ""
            self.collection.update_one({'_id': doc['_id']},
                                       {'$set': {'protein_name': protein_name}})
    
    def fill_uniprot_by_oln(self, oln):
        """Fill uniprot collection using ordered locus name
        
        Args:
            oln (:obj:`str`): Ordered locus name
        """
        gene_name, _ = self.uniprot_manager.get_gene_protein_name_by_oln(oln)
        if gene_name == "": # no such entry in uniprot collection
            self.uniprot_col_manager.load_uniprot(query=True, msg=oln)
        else:
            return

    def fill_gene_protein_name(self):
        """Fill 'gene_name' field where 'gene_name' has value of '-' and create
        'protein_name' field
        """
        query = {'gene_name': '-'}
        projection = {'gene_name': 1, 'halflives': 1}
        gene_name = ''
        protein_name = ''
        docs = self.collection.find(filter=query, projection=projection, collation=self.collation)
        for doc in docs:
            oln = doc['halflives'][0]['ordered_locus_name']
            gene_name, protein_name = self.uniprot_manager.get_gene_protein_name_by_oln(oln)
            self.collection.update_one({'_id': doc['_id']},
                                       {'$set': {'gene_name': gene_name,
                                                 'protein_name': protein_name}})

    def uniprot_names(self, results, count):
        """Extract protein_name and gene_name from returned
        tuple of uniprot query function
        
        Args:
            results (:obj:`Iter`): pymongo cursor object.
            count (:obj:`int`): Number of documents found.

        Return:
            (:obj:`tuple` of :obj:`str`): gene_name and protein_name
        """
        if count == 0:
            return '', ''
        else:
            for result in results:
                gene_name = result['gene_name']
                protein_name = result['protein_name']
                return gene_name, protein_name

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