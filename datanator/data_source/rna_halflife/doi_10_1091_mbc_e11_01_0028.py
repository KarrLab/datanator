import pandas as pd
from datanator.util import rna_halflife_util, file_util
from datanator_query_python.query import query_uniprot
import datetime
import datanator.config.core
import datetime
from pymongo.collation import Collation, CollationStrength
from pymongo.errors import WriteError
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
        self.uniprot_query = query_uniprot.QueryUniprot(username=username, password=password, server=server,
        authSource=authDB, collection_str='uniprot', readPreference=readPreference)
        self.collation = Collation('en', strength=CollationStrength.SECONDARY)
        self.max_entries = max_entries
        self.verbose = verbose

    def fill_rna_halflife(self, df, start=0):
        """Fill rna_halflife collection with information parsed
        from the publication
        
        Args:
            df (:obj:`pandas.DataFrame`): dataframe to be loaded.
            start (:obj:`int`, optional): Starting row in df. Defaults to 0.
        """
        row_count = len(df.index)
        for i, row in df.iloc[start:].iterrows():
            if i == self.max_entries:
                break
            if i % 10 == 0 and self.verbose:
                print("Processing locus {} out {}".format(i, row_count))
            systematic_name = row['systematic_name']
            halflife = row['halflife'] * 60
            unit = 's'
            species = 'Saccharomyces cerevisiae W303'
            r_squared = row['r_squared']
            ncbi_taxonomy_id = 580240
            reference = [{'doi': '10.1091/mbc.e11-01-0028'}]
            obj = {'systematic_name': systematic_name,
                   'halflife': halflife,
                   'unit': unit,
                   'species': species,
                   'ncbi_taxonomy_id': ncbi_taxonomy_id,
                   'r_squared': r_squared,
                   'reference': reference}
            try:
                self.rna_hl_collection.update_one({'halflives.systematic_name': systematic_name},
                                                {'$addToSet': {'halflives': obj}},
                                                collation=self.collation, upsert=True)
            except WriteError:
                protein_name = self.uniprot_query.get_protein_name_by_gn(systematic_name, species=[559292])
                self.rna_hl_collection.insert_one({'gene_name': systematic_name,
                                                   'protein_name': protein_name,
                                                   'halflives': [obj]})

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
    # url = 'https://www.molbiolcell.org/doi/suppl/10.1091/mbc.e11-01-0028/suppl_file/mc-e11-01-0028-s10.xls'
    url = os.path.expanduser('~/karr_lab/datanator/docs/mc-e11-01-0028-s10.xls')
    names = ['systematic_name', 'halflife', 'r_squared']
    df_s1 = src.make_df(url, 'all mRNAs (with R2>0.8) ranked', names=names, usecols='A:C',
                        file_type='xls', file_name='mc-e11-01-0028-s10.xls')
    src.fill_rna_halflife(df_s1)
    shutil.rmtree(cache_dir)

if __name__ == '__main__':
    main()