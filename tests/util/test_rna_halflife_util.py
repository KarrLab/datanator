import unittest
from datanator.util import rna_halflife_util
from datanator_query_python.config import config


class TestRnaHlUtil(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        src_db = 'datanator'
        des_db = 'test'
        conf = config.TestConfig()
        username = conf.MONGO_TEST_USERNAME
        password = conf.MONGO_TEST_PASSWORD
        MongoDB = conf.SERVER
        username = username
        password = password
        cls.protein_col = 'uniprot'
        cls.rna_col = 'rna_halflife'
        cls.src = rna_halflife_util.RnaHLUtil(server=MongoDB, username=username,
        password=password, src_db=src_db, des_db=des_db, protein_col=cls.protein_col,
        rna_col=cls.rna_col, readPreference='nearest')

    @classmethod
    def tearDownClass(cls):
        cls.src.uniprot_collection_manager.db.drop_collection(cls.protein_col)
        cls.src.uniprot_collection_manager.client.close()
        cls.src.uniprot_query_manager.client.close()

    @unittest.skip('avoid r/w db')
    def test_fill_uniprot_by_oln(self):
        self.src.fill_uniprot_by_oln('MA0002')

    def test_make_df(self):
        url_0 = 'https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/nar/41/1/10.1093/nar/gks1019/2/gks1019-nar-00676-a-2012-File003.xlsx?Expires=1578425844&Signature=ZRFUxLdn4-vaBt5gQci~0o56KqyR9nJj9i32ig5X6YcfqiJeV3obEq8leHGdDxx6w~KABgewiQ66HTB7gmuG~2GL-YgxPKYSjt17WrYMkc-0ibw6TMlTvWZZfvw-lPe~wvpmVfNEXnTbP7jHyNLu9jeJ6yhoXvgIyQtzA5PbEI1fyXEgeZzOKMltmITqL3g3APsPsagCTC66rwrBT23Aghh6D314uilT2DZHCc68MH2nyV~qAhFqIQiOj-7VTEKqkDPvPYvuE2KNKXdvW23gk100YV~58ozbt8ijRz5Gr5gPtE~f1Ab5l260EIbWHJNabMRleInJQqUIDPFN4C38PQ__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA'
        df_0 = self.src.make_df(url_0, 'Supplementary Table 1', usecols='B:D', skiprows=[0,1,2],
        names=['ordered_locus_name', 'half_life', 'r_squared'])
        self.assertEqual(df_0.iloc[0]['ordered_locus_name'], 'Rv0002')
