import unittest
from datanator.util import mongo_util
import datanator.config.core
import tempfile
import shutil


class TestMongoUtil(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'datanator'
        username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
        port = datanator.config.core.get_config()['datanator']['mongodb']['port']
        replSet = datanator.config.core.get_config()['datanator']['mongodb']['replSet']
        cls.src = mongo_util.MongoUtil(
            cache_dirname = cls.cache_dirname, MongoDB = MongoDB, 
            replicaSet = replSet, db = cls.db, verbose=True, max_entries=20,
            username = username, password = password)
        cls.collection_str = 'ecmdb'


    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)


    # @unittest.skip('passed')
    def test_list_all_collections(self):
        self.assertTrue('ecmdb' in self.src.list_all_collections())


    # @unittest.skip('passed')
    def test_con_db(self):
        self.assertNotEqual(self.src.con_db(self.db), 'Server not available')

    @unittest.skip('passed')
    def test_fill_db(self):
        self.collection_obj.drop()
        self.assertEqual(self.collection_obj.find().count(), 0)
        collection_obj = self.src.fill_db(self.collection_str)
        self.assertNotEqual(collection_obj.find().count(), 0)

    # @unittest.skip('passed')
    def test_print_schema(self):
        a = self.src.print_schema('ecmdb')
        self.assertEqual(a['properties']['creation_date'], {'type': 'string'})
        self.assertEqual(a['properties']['synonyms'],  {'type': 'object', 'properties': {'synonym': {'type': 'array', 
            'items': {'type': 'string'}}}, 'required': ['synonym']})


