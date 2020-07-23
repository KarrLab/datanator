import unittest
import shutil
import tempfile
from datanator.data_source.brenda import reaction
from datanator_query_python.config import config
import pandas as pd


class TestBrendaRxn(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        conf = config.TestConfig()
        cls.collection_str = 'brenda_reaction'
        username = conf.USERNAME
        password = conf.PASSWORD
        MongoDB = conf.SERVER
        cls.src = reaction.BrendaRxn(MongoDB=MongoDB, db='test', collection_str=cls.collection_str,
                                     username=username, password=password, authSource='admin',
                                     max_entries=20, verbose=True)

    @classmethod
    def tearDownClass(cls):
        cls.src.db_obj.drop_collection(cls.collection_str)
        cls.src.client.close()

    # @unittest.skip('passed')
    def test_download_and_read(self):
        result = self.src.download_and_read()
        self.assertEqual(result['ec_number'][1], '6.3.2.1')

    def test_clean_up(self):
        result = self.src.download_and_read()
        exp = self.src.clean_up(result)
        self.assertEqual(exp['reaction_id_brenda'][1], ['BR101'])
        self.assertEqual(exp['reaction_id_sabio_rk'][1], 2406)

    # @unittest.skip('passed')
    def test_parse_reaction(self):
        df = pd.DataFrame({'reaction': ['ATP + (R)-pantoate + beta-alanine <=> AMP + diphosphate + (R)-pantothenate',
                                        'ATP + Detyrosinated alpha-tubulin + L-Tyrosine = alpha-Tubulin + ADP + Orthophosphate']})
        result = self.src.parse_reaction(df)
        self.assertEqual(result['products'][1][1], 'ADP')
        self.assertEqual(result['substrates'][0][1], '(R)-pantoate')

    # @unittest.skip('passed')
    def test_load_df_sim(self):
        df = pd.DataFrame({'reaction': ['ATP + (R)-pantoate + beta-alanine <=> AMP + diphosphate + (R)-pantothenate',
                                        'ATP + Detyrosinated alpha-tubulin + L-Tyrosine = alpha-Tubulin + ADP + Orthophosphate']})
        result = self.src.parse_reaction(df)
        self.src.load_df(result)

    # @unittest.skip('passed')
    def test_load_df_real(self):
        result = self.src.download_and_read()
        self.src.clean_up(result)
        x = self.src.parse_reaction(result)
        self.src.load_df(x.head(100))