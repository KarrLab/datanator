import unittest
import tempfile
import shutil
from datanator.util import rna_halflife_util
from datanator_query_python.config import config
import pandas as pd


class TestRnaHlUtil(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        src_db = 'datanator'
        des_db = 'test'
        conf = config.TestConfig()
        username = conf.USERNAME
        password = conf.PASSWORD
        MongoDB = conf.SERVER
        username = username
        password = password
        cls.cache_dir = tempfile.mkdtemp()
        cls.protein_col = 'uniprot'
        cls.rna_col = 'rna_halflife'
        cls.src = rna_halflife_util.RnaHLUtil(server=MongoDB, username=username,
        password=password, src_db=src_db, des_db=des_db, protein_col=cls.protein_col,
        rna_col=cls.rna_col, readPreference='nearest', cache_dir=cls.cache_dir)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dir)
        cls.src.uniprot_collection_manager.db.drop_collection(cls.protein_col)
        cls.src.uniprot_collection_manager.client.close()
        cls.src.uniprot_query_manager.client.close()

    @unittest.skip('avoid r/w db')
    def test_fill_uniprot_by_oln(self):
        self.src.fill_uniprot_by_oln('MA0002')

    @unittest.skip('links will not work from nonacademic IPs.')
    def test_make_df(self):
        url_0 = 'https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/nar/41/1/10.1093/nar/gks1019/2/gks1019-nar-00676-a-2012-File003.xlsx?Expires=1578425844&Signature=ZRFUxLdn4-vaBt5gQci~0o56KqyR9nJj9i32ig5X6YcfqiJeV3obEq8leHGdDxx6w~KABgewiQ66HTB7gmuG~2GL-YgxPKYSjt17WrYMkc-0ibw6TMlTvWZZfvw-lPe~wvpmVfNEXnTbP7jHyNLu9jeJ6yhoXvgIyQtzA5PbEI1fyXEgeZzOKMltmITqL3g3APsPsagCTC66rwrBT23Aghh6D314uilT2DZHCc68MH2nyV~qAhFqIQiOj-7VTEKqkDPvPYvuE2KNKXdvW23gk100YV~58ozbt8ijRz5Gr5gPtE~f1Ab5l260EIbWHJNabMRleInJQqUIDPFN4C38PQ__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA'
        df_0 = self.src.make_df(url_0, 'Supplementary Table 1', usecols='B:D', skiprows=[0,1,2],
        names=['ordered_locus_name', 'half_life', 'r_squared'])
        self.assertEqual(df_0.iloc[0]['ordered_locus_name'], 'Rv0002')
        url_1 = 'https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/nar/42/4/10.1093_nar_gkt1150/1/gkt1150_Supplementary_Data.zip?Expires=1578928721&Signature=ADjsCSaceimzGs6aJ~uG7np88TzHNooAoBabdm-6utYVIZOEwRbzTdiBp~76vM4KEHz9Nir8GNrtA3AwHwGFm0bu~aorTG4xrOChS6UgfBQiUtgr8vfbDIUno1y1nxLGCKIfQrb2Gx-SVnigum2gjcveymK995zadSNZqN~w-vz-Ii0a6fH7kvKN8m9vLWf6fdo0NXSmgnkjj9KPCuS-bmK0y4ZH5Ex0Rl4qi5uCroYmDBNOhXY23pcalbpFwB1-07tA3~756gZN4Mo9uMeSVQKl5UsHzx5amB6WvSCXS8z756YoaaMCg0FQbUCcQ46fRGdHxcvPNcrPo5IMEGmi8g__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA'
        df_1 = self.src.make_df(url_1, 'TableS1', file_type='zip', file_name='nar-01935-a-2013-File011.xlsx', usecols='L:O', skiprows=list(range(0, 7)),
                                names=['a', 'b', 'c', 'd'])
        self.assertEqual(df_1.iloc[0]['a'], 5.74239011770224)

    def test_fill_uniprot_with_df(self):
        pass