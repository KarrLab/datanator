import unittest
from datanator.data_source.rna_halflife import order_by_ko
from datanator_query_python.config import config


class TestReorg(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        des_db = 'test'
        src_db = 'datanator'
        cls.src_collection = 'rna_halflife'
        cls.des_collection = 'rna_halflife_test'
        conf = config.TestConfig()
        username = conf.MONGO_TEST_USERNAME
        password = conf.MONGO_TEST_PASSWORD
        MongoDB = conf.SERVER    
        cls.src = order_by_ko.Reorg(MongoDB=MongoDB, src_db=src_db,
                 verbose=False, max_entries=20, username=username, 
                 password=password, authSource='admin', readPreference='nearest',
                 des_collection=cls.des_collection, src_collection=cls.src_collection,
                 des_db=des_db)

    @classmethod
    def tearDownClass(cls):
        cls.src.des_db.drop_collection(cls.des_collection)
        cls.src.src_client.close()
        cls.src.des_client.close()

    def test_helper(self):
        doi = '10.1016/j.cell.2013.12.026'
        _, count = self.src.helper(doi)
        self.assertGreater(count, 0)

    @unittest.skip('passed')
    def test_fill_cell(self):
        self.src.fill_cell()

    @unittest.skip('passed')
    def test_fill_mbc(self):
        self.src.fill_mbc()

    @unittest.skip('passed')
    def test_fill_nar_gks(self):
        self.src.fill_nar_gks()

    @unittest.skip('passed')
    def test_fill_nar_gkt(self):
        self.src.fill_nar_gkt()

    @unittest.skip('passed')
    def test_fill_gr_131(self):
        self.src.fill_gr_131()

    @unittest.skip('passed')
    def test_fill_gb_2012(self):
        self.src.fill_gb_2012()

    @unittest.skip('passed')
    def test_fill_s12864(self):
        self.src.fill_s12864()

    @unittest.skip('passed')
    def test_fill_journal_pone(self):
        self.src.fill_journal_pone()         