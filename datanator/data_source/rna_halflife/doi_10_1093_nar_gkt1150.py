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

    def fill_rna_half_life(self, df, species):
        """Load df into rna_halflife collection
        
        Args:
            df (:obj:`pandas.DataFrame`): dataframe to be loaded
            species (:obj:`list`): species name and ncbi_id.
        """
        row_count = len(df.index)
        for i, row in df.iterrows():
            if i == self.max_entries:
                break
            if i % 10 == 0 and self.verbose:
                print("Processing locus {} out {}".format(i, row_count))
            if row[2:].isnull().values.all():
                continue
            oln = row['oln']
            symbol = row['gene_symbol']
            a = row['a'] * 60
            vc_a = row['vc_a']
            b = row['b'] * 60
            vc_b = row['vc_b']
            c = row['c'] * 60
            vc_c = row['vc_c']
            d = row['d'] * 60
            vc_d = row['vc_d']
            dic_0 = {'halflife': a, 'variation_coefficient': vc_a, 'species': species[0], 'ncbi_taxonomy_id': species[1],
                    'unit': 's', 'reference': [{'doi': '10.1093/nar/gkt1150'}], 'growth_medium': 'M9 minimal medium supplemented with glucose',
                    'ordered_locus_name': oln, 'doubling_time': {'value': 6.9, 'unit': 'h'}}
            dic_1 = {'halflife': b, 'variation_coefficient': vc_b, 'species': species[0], 'ncbi_taxonomy_id': species[1],
                    'unit': 's', 'reference': [{'doi': '10.1093/nar/gkt1150'}], 'growth_medium': 'M9 minimal medium supplemented with glucose',
                    'ordered_locus_name': oln, 'doubling_time': {'value': 3.5, 'unit': 'h'}}
            dic_2 = {'halflife': c, 'variation_coefficient': vc_c, 'species': species[0], 'ncbi_taxonomy_id': species[1],
                    'unit': 's', 'reference': [{'doi': '10.1093/nar/gkt1150'}], 'growth_medium': 'M9 minimal medium supplemented with glucose',
                    'ordered_locus_name': oln, 'doubling_time': {'value': 2.3, 'unit': 'h'}}
            dic_3 = {'halflife': d, 'variation_coefficient': vc_d, 'species': species[0], 'ncbi_taxonomy_id': species[1],
                    'unit': 's', 'reference': [{'doi': '10.1093/nar/gkt1150'}], 'growth_medium': 'M9 minimal medium supplemented with glucose',
                    'ordered_locus_name': oln, 'doubling_time': {'value': 1.7, 'unit': 'h'}}
            halflives = [dic_0, dic_1, dic_2, dic_3]
            # record is guaranteed to exist in uniprot because uniprot was filled prior to this operation
            gene_name, protein_name = self.uniprot_query_manager.get_gene_protein_name_by_oln(oln, species=[83333, 511145])
            if gene_name is None:
                gene_name = symbol

            self.rna_hl_collection.update_one({'gene_name': gene_name},
                                        {'$addToSet': {'halflives': {'$each': halflives},
                                                        'protein_synonyms': protein_name}},
                                        collation=self.collation, upsert=True)



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
    url = 'https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/nar/42/4/10.1093_nar_gkt1150/1/gkt1150_Supplementary_Data.zip?Expires=1578928721&Signature=ADjsCSaceimzGs6aJ~uG7np88TzHNooAoBabdm-6utYVIZOEwRbzTdiBp~76vM4KEHz9Nir8GNrtA3AwHwGFm0bu~aorTG4xrOChS6UgfBQiUtgr8vfbDIUno1y1nxLGCKIfQrb2Gx-SVnigum2gjcveymK995zadSNZqN~w-vz-Ii0a6fH7kvKN8m9vLWf6fdo0NXSmgnkjj9KPCuS-bmK0y4ZH5Ex0Rl4qi5uCroYmDBNOhXY23pcalbpFwB1-07tA3~756gZN4Mo9uMeSVQKl5UsHzx5amB6WvSCXS8z756YoaaMCg0FQbUCcQ46fRGdHxcvPNcrPo5IMEGmi8g__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA'
    df_s1 = src.make_df(url, 'TableS1', names=['oln', 'gene_symbol', 'a', 'vc_a', 'b', 'vc_b', 'c', 'vc_c', 'd', 'vc_d'], usecols='A,B,L:S',
                        skiprows=list(range(0, 7)), file_type='zip', file_name='nar-01935-a-2013-File011.xlsx')
    # src.fill_uniprot_with_df(df_s1, 'oln', species=[83333, 511145])
    src.fill_rna_half_life(df_s1, ['Escherichia coli str. K-12 substr. MG1655', 511145])
    shutil.rmtree(cache_dir)

if __name__ == '__main__':
    main()