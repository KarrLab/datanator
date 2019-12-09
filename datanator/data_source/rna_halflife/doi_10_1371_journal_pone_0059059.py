import pandas as pd
from datanator.util import rna_halflife_util, file_util
import datetime
import datanator.config.core
import datetime
from pymongo.collation import Collation, CollationStrength
import tempfile
import shutil
import tabula


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
    df = tabula.read_pdf(url, pandas_options={'header': None})
    df.columns = ['gene_name', 'a', 'b', 'c', 'd', 'e']
    # src.fill_uniprot_with_df(df_s1, 'oln', species=[83333, 511145])
    shutil.rmtree(cache_dir)

if __name__ == '__main__':
    main()