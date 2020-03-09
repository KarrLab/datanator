import pandas as pd
from datanator.util import rna_halflife_util, file_util
import datetime
import datanator.config.core
import datetime
from pymongo.collation import Collation, CollationStrength
import tempfile
import shutil
import tabula
import math


class Halflife(rna_halflife_util.RnaHLUtil):

    def __init__(self, cache_dir=None, server=None, src_db=None, protein_col=None,
                authDB=None, readPreference=None, username=None, password=None,
                verbose=None, max_entries=None, des_db=None, rna_col=None):
        """Init
        
        Args:
            cache_dir (:obj:`str`, optional): Cache directory for logs. Defaults to None.
            server (:obj:`str`, optional): MongoDB server address. Defaults to None.
            db (:obj:`str`, optional): Database where initial uniprot collection resides. Defaults to None.
            collection_str (:obj:`str`, optional): name of collection. Defaults to None.
            authDB (:obj:`str`, optional): MongoDB authentication database. Defaults to None.
            readPreference (:obj:`str`, optional): MongoDB read preference. Defaults to None.
            username (:obj:`str`, optional): MongoDB username. Defaults to None.
            password (:obj:`str`, optional): MongoDB password. Defaults to None.
            verbose (:obj:`bool`, optional): Wheter to display verbose messages. Defaults to None.
            max_entries (:obj:`int`, optional): Number of records to be processed. Defaults to None.
            uniprot_col_db (:obj:`int`, optional): Database to which new uniprot records will be inserted. Defaults to None.
        """
        super().__init__(server=server, username=username, password=password, src_db=src_db,
        des_db=des_db, protein_col=protein_col, rna_col=rna_col, authDB=authDB, readPreference=readPreference,
        max_entries=max_entries, verbose=verbose, cache_dir=cache_dir)
        self.collation = Collation('en', strength=CollationStrength.SECONDARY)
        self.max_entries = max_entries
        self.verbose = verbose

    def fill_rna_half_life(self, df, species):
        """Load df into rna_halflife collection
        
        Args:
            df (:obj:`pandas.DataFrame`): dataframe to be loaded
            species (:obj:`list`): species name and ncbi_id.
        """
        d_0 = math.log(2) / 0.11
        d_1 = math.log(2) / 0.51
        d_2 = math.log(2) / 0.80
        d_3 = math.log(2) / 0.38
        d_4 = math.log(2) / 0.04
        reference = [{'doi': '10.1371/journal.pone.0059059'}]
        row_count = len(df.index)
        for i, row in df.iterrows():
            if i == self.max_entries:
                break
            if i % 10 == 0 and self.verbose:
                print("Processing locus {} out {}".format(i, row_count))
            if row[1:].isnull().values.all():
                continue
            symbol = row['gene_name']
            a = self.hl_helper(row['a'])
            b = self.hl_helper(row['b'])
            c = self.hl_helper(row['c'])
            d = self.hl_helper(row['d'])
            e = self.hl_helper(row['e'])
            dic_0 = {'halflife': a, 'species': species[0], 'ncbi_taxonomy_id': species[1],
                    'unit': 's', 'reference': reference,
                    'doubling_time': {'value': d_0, 'unit': 'h'}}
            dic_1 = {'halflife': b, 'species': species[0], 'ncbi_taxonomy_id': species[1],
                    'unit': 's', 'reference': reference, 
                    'doubling_time': {'value': d_1, 'unit': 'h'}}
            dic_2 = {'halflife': c, 'species': species[0], 'ncbi_taxonomy_id': species[1],
                    'unit': 's', 'reference': reference, 
                    'doubling_time': {'value': d_2, 'unit': 'h'}}
            dic_3 = {'halflife': d, 'species': species[0], 'ncbi_taxonomy_id': species[1],
                    'unit': 's', 'reference': reference, 
                    'doubling_time': {'value': d_3, 'unit': 'h'}}
            dic_4 = {'halflife': e, 'species': species[0], 'ncbi_taxonomy_id': species[1],
                    'unit': 's', 'reference': reference, 
                    'doubling_time': {'value': d_4, 'unit': 'h'}}
            halflives = [dic_0, dic_1, dic_2, dic_3, dic_4]
            # record is guaranteed to exist in uniprot because uniprot was filled prior to this operation
            protein_name = self.uniprot_query_manager.get_protein_name_by_gn(symbol, species=[272623])
            if protein_name is None:
                continue

            self.rna_hl_collection.update_one({'gene_name': symbol},
                                        {'$addToSet': {'halflives': {'$each': halflives},
                                                        'protein_synonyms': protein_name}},
                                        collation=self.collation, upsert=True)

    def hl_helper(self, value):
        """Process halflife value read from pdf
        
        Args:
            value (:obj:`str`): value from pdf
        """
        if value == 'very stable':
            return 0.
        else:
            return float(value) * 60



def main():
    src_db = 'datanator'
    des_db = 'datanator'
    rna_col = 'rna_halflife'
    protein_col = 'uniprot'
    cache_dir = tempfile.mkdtemp()
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    server = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']
    src = Halflife(server=server, src_db=src_db,
        protein_col=protein_col, authDB='admin', readPreference='nearest',
        username=username, password=password, verbose=True, max_entries=float('inf'),
        des_db=des_db, rna_col=rna_col, cache_dir=cache_dir)
    url = 'https://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0059059.s002'
    df = tabula.read_pdf(url, pandas_options={'header': None, 'na_values': 'ND'}, pages='all')
    df.columns = ['gene_name', 'a', 'b', 'c', 'd', 'e']
    # src.fill_uniprot_with_df(df, 'gene_name', identifier_type='gene_name', species=[272623])
    src.fill_rna_half_life(df, ['Lactococcus lactis subsp. lactis Il1403', 272623])
    shutil.rmtree(cache_dir)

if __name__ == '__main__':
    main()