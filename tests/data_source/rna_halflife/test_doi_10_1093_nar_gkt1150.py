import unittest
from datanator.data_source.rna_halflife import doi_10_1093_nar_gkt1150
import tempfile
import shutil
import json
import os
from datanator_query_python.config import config
import pandas as pd


class TestProteinAggregate(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        des_db = 'test'
        src_db = 'datanator'
        cls.protein_col = 'uniprot'
        cls.rna_col = 'rna_halflife'
        cls.cache_dir = tempfile.mkdtemp()
        conf = config.TestConfig()
        username = conf.USERNAME
        password = conf.PASSWORD
        MongoDB = conf.SERVER    
        cls.src = doi_10_1093_nar_gkt1150.Halflife(server=MongoDB, src_db=src_db,
        protein_col=cls.protein_col, authDB='admin', readPreference='nearest',
        username=username, password=password, verbose=True, max_entries=20,
        des_db=des_db, rna_col=cls.rna_col, cache_dir=cls.cache_dir)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dir)
        cls.src.uniprot_collection_manager.db.drop_collection(cls.protein_col)
        cls.src.db_obj.drop_collection(cls.rna_col)
        cls.src.uniprot_collection_manager.client.close()
        cls.src.client.close()
        cls.src.uniprot_query_manager.client.close()

    @unittest.skip('downloading of file forbidden from nonacademic IP')
    def test_fill_rna_half_life(self):
        url = 'https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/nar/42/4/10.1093_nar_gkt1150/1/gkt1150_Supplementary_Data.zip?Expires=1578928721&Signature=ADjsCSaceimzGs6aJ~uG7np88TzHNooAoBabdm-6utYVIZOEwRbzTdiBp~76vM4KEHz9Nir8GNrtA3AwHwGFm0bu~aorTG4xrOChS6UgfBQiUtgr8vfbDIUno1y1nxLGCKIfQrb2Gx-SVnigum2gjcveymK995zadSNZqN~w-vz-Ii0a6fH7kvKN8m9vLWf6fdo0NXSmgnkjj9KPCuS-bmK0y4ZH5Ex0Rl4qi5uCroYmDBNOhXY23pcalbpFwB1-07tA3~756gZN4Mo9uMeSVQKl5UsHzx5amB6WvSCXS8z756YoaaMCg0FQbUCcQ46fRGdHxcvPNcrPo5IMEGmi8g__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA'
        df_s1 = self.src.make_df(url, 'TableS1', names=['oln', 'gene_symbol', 'a', 'vc_a', 'b', 'vc_b', 'c', 'vc_c', 'd', 'vc_d'], usecols='A,B,L:S',
                            skiprows=list(range(0, 7)), file_type='zip', file_name='nar-01935-a-2013-File011.xlsx')
        self.src.fill_rna_half_life(df_s1, ['Escherichia coli str. K-12 substr. MG1655', 511145])