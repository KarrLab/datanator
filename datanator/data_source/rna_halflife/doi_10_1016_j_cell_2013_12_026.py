import pandas as pd
from datanator.util import rna_halflife_util, file_util
import datetime
import datanator.config.core
import datetime
from pymongo.collation import Collation, CollationStrength
import tempfile
import shutil


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

    def fill_rna_halflife(self, df):
        """Fill RNA Half life collection with relevant documents.

        Args:
            df (:obj:`pandas.DataFrame`): dataframe to be loaded.
        """
        row_count = len(df.index)
        for i, row in df.iterrows():
            if i == self.max_entries:
                break
            if i % 10 == 0 and self.verbose:
                print("Processing locus {} out {}".format(i, row_count))
            chrom = row['chromosome']
            systematic_name = row['sys_name']
            gene_name = row['gene_name']
            _type = row['cds']
            halflife = row['half_life'] * 60
            description = row['description']
            halflife_obj = {'chromosome': chrom,
                        'systematic_name': systematic_name,
                        'gene_name': gene_name,
                        'type': 'coding_dna_sequences',
                        'halflife': halflife,
                        'unit': 's',
                        'species': 'Saccharomyces cerevisiae S288C',
                        'ncbi_taxonomy_id': 559292,
                        'reference': [{'doi': '10.1016/j.cell.2013.12.026'}]}
            self.rna_hl_collection.update_one({'gene_name': gene_name},
                                        {'$addToSet': {'halflives': halflife_obj,
                                                       'description': description}},
                                        collation=self.collation, upsert=True)

    # def fill_doi(self):
    #     """Only need to run once b/c I forgot to include reference information.
    #     """
    #     self.rna_hl_collection.update_many({},
    #                                        {'$pull': {'halflives': {'ncbi_taxonomy_id': 4932}}})



import os
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
    # url = 'https://ars.els-cdn.com/content/image/1-s2.0-S009286741301595X-mmc1.xlsx'
    url = os.path.expanduser('~/karr_lab/datanator/docs/1-s2.0-S009286741301595X-mmc1.xlsx')
    names = ['chromosome', 'sys_name', 'gene_name', 'cds', 'half_life', 'description']
    df_s1 = src.make_df(url, 'Sheet1', names=names, usecols='A:F',
                        file_type='xlsx', file_name='1-s2.0-S009286741301595X-mmc1.xlsx')
    src.fill_rna_halflife(df_s1)
    shutil.rmtree(cache_dir)

if __name__ == '__main__':
    main()