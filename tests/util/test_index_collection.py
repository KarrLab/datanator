import unittest
from datanator.util import index_collection
import tempfile
import shutil


class TestIndexCollection(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'test'
        cls.MongoDB = 'mongodb://mongo:27017/'
        cls.src = index_collection.IndexCollection(
            cls.cache_dirname, cls.MongoDB, None, cls.db, verbose=True, max_entries=20)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    @unittest.skip('passed')
    def test_index_corum(self):
        col_str = 'corum'
        self.src.index_corum(col_str)
        client, _, collection = self.src.con_db(col_str)
        self.assertEqual(len(list(collection.list_indexes())), 4)
        client.close()

    @unittest.skip('passed')
    def test_index_sabio(self):
        col_str = 'sabio_rk'
        self.src.index_sabio(col_str)
        client,_,collection = self.src.con_db(col_str)
        self.assertEqual( len(list(collection.list_indexes())), 11) # 10 + 1
        client.close()

    def test_index_uniprot(self):
        col_str = 'uniprot'
        self.src.index_uniprot(col_str)
        client,_,collection = self.src.con_db(col_str)
        self.assertEqual( len(list(collection.list_indexes())), 3) # 2 + 1
        client.close()