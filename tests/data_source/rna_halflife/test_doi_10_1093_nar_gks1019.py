import unittest
from datanator.data_source.rna_halflife import doi_10_1093_nar_gks1019
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
        conf = config.TestConfig()
        username = conf.MONGO_TEST_USERNAME
        password = conf.MONGO_TEST_PASSWORD
        MongoDB = conf.SERVER    
        cls.src = doi_10_1093_nar_gks1019.Halflife(server=MongoDB, src_db=src_db,
        protein_col=cls.protein_col, authDB='admin', readPreference='nearest',
        username=username, password=password, verbose=True, max_entries=20,
        des_db=des_db, rna_col=cls.rna_col)

    @classmethod
    def tearDownClass(cls):
        # cls.src.uniprot_collection_manager.db.drop_collection(cls.protein_col)
        # cls.src.db_obj.drop_collection(cls.rna_col)
        cls.src.uniprot_collection_manager.client.close()
        cls.src.client.close()
        cls.src.uniprot_query_manager.client.close()

    @unittest.skip('avoid downloading')
    def test_fill_uniprot(self):
        url_0 = 'https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/nar/41/1/10.1093/nar/gks1019/2/gks1019-nar-00676-a-2012-File003.xlsx?Expires=1578425844&Signature=ZRFUxLdn4-vaBt5gQci~0o56KqyR9nJj9i32ig5X6YcfqiJeV3obEq8leHGdDxx6w~KABgewiQ66HTB7gmuG~2GL-YgxPKYSjt17WrYMkc-0ibw6TMlTvWZZfvw-lPe~wvpmVfNEXnTbP7jHyNLu9jeJ6yhoXvgIyQtzA5PbEI1fyXEgeZzOKMltmITqL3g3APsPsagCTC66rwrBT23Aghh6D314uilT2DZHCc68MH2nyV~qAhFqIQiOj-7VTEKqkDPvPYvuE2KNKXdvW23gk100YV~58ozbt8ijRz5Gr5gPtE~f1Ab5l260EIbWHJNabMRleInJQqUIDPFN4C38PQ__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA'
        df_0 = self.src.fill_uniprot(url_0, 'Supplementary Table 1')
        self.assertEqual(df_0.iloc[0]['ordered_locus_name'], 'Rv0002')

    def test_fill_rna_halflife(self):
        d = {'half_life': [32.3, 12.2, 13.2], 'r_sqaured': [0.9, 0.7, 0.8],
            'ordered_locus_name': ['Rv0002', 'something', 'this']}
        df_0 = pd.DataFrame(d)
        self.src.fill_rna_halflife(df_0, ['aaa', 102])