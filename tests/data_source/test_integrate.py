import unittest
from datanator.data_source import integrate
import tempfile
import shutil


class TestIntegrate(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        # '../../datanator/data_source/cache'
        cls.db = 'datanator'
        cls.MongoDB = 'mongodb://mongo:27017/'
        cls.src = integrate.Integrate(
            cls.cache_dirname, cls.MongoDB, cls.db, verbose=True, max_entries=20)
        cls.client, cls.db = cls.src.con_db()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.client.close()

    @unittest.skip('passed')
    def test_con_db(self):
        self.assertNotEqual(self.db, 'Server not available')

    @unittest.skip('passed')
    def test_retrieve_sabio(self):
        collection = self.db['sabio_rk']
        taxon = 562  # E. Coli taxon id
        cursor = self.src.retrieve_sabio(collection, taxon)
        self.assertEqual(cursor.count(with_limit_and_skip=True), 20)
        self.assertEqual(cursor[0]['kinlaw_id'], 47)
        self.assertEqual(cursor[0]['parameter'][0]['observed_value'], 0.074)

    def test_find_chebi(self):
        collection = self.db['sabio_rk']
        taxon = 562  # E. Coli taxon id
        cursor = self.src.retrieve_sabio(collection, taxon)
        substrate_chebi, product_chebi = self.src.find_sabio_chebi(cursor)
        self.assertTrue('C00019' in substrate_chebi)

