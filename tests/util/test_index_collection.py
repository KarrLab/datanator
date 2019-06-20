import unittest
from datanator.util import index_collection
import datanator.config.core
import tempfile
import shutil


class TestMongoUtil(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'test'
        username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
        port = datanator.config.core.get_config()['datanator']['mongodb']['port']
        replSet = datanator.config.core.get_config()['datanator']['mongodb']['replSet']
        cls.src = index_collection.IndexCollection(
            cache_dirname = cls.cache_dirname, MongoDB = MongoDB, 
            replicaSet = replSet, db = cls.db, verbose=True, max_entries=20,
            username = username, password = password)

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

    @unittest.skip('passed')
    def test_index_uniprot(self):
        col_str = 'uniprot'
        self.src.index_uniprot(col_str)
        client,_,collection = self.src.con_db(col_str)
        self.assertEqual( len(list(collection.list_indexes())), 3) # 2 + 1
        client.close()