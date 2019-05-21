import unittest
from datanator.util import mongo_util
import tempfile
import shutil


class TestMongoUtil(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'datanator'
        cls.MongoDB = 'mongodb://mongo:27017/'
        cls.src = mongo_util.MongoUtil(
            cls.cache_dirname, cls.MongoDB, 'rs0', cls.db, verbose=True, max_entries=20)
        cls.collection_str = 'ecmdb'
        cls.client, cls.db, cls.collection_obj = cls.src.con_db(
            cls.collection_str)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        # cls.client.drop_database(cls.db)
        cls.client.close()

    @unittest.skip('passed')
    def test_list_all_collections(self):
        self.assertTrue('ecmdb' in self.src.list_all_collections())


    @unittest.skip('passed')
    def test_con_db(self):
        self.assertNotEqual(self.db, 'Server not available')

    @unittest.skip('passed')
    def test_fill_db(self):
        self.collection_obj.drop()
        self.assertEqual(self.collection_obj.find().count(), 0)
        collection_obj = self.src.fill_db(self.collection_str)
        self.assertNotEqual(collection_obj.find().count(), 0)

    @unittest.skip('passed')
    def test_print_schema(self):
        a = self.src.print_schema('ecmdb')
        self.assertEqual(a['properties']['creation_date'], {'type': 'string'})
        self.assertEqual(a['properties']['synonyms'],  {'type': 'object', 'properties': {'synonym': {'type': 'array', 
            'items': {'type': 'string'}}}, 'required': ['synonym']})

    @unittest.skip('passed')
    def test_flatten_json(self):
        dic = {
            "a": 1,
            "b": 2,
            "c": [{"d": [2, 3, 4], "e": [{"f": 1, "g": 2}]}],
            'h': [1, 2, 3]
        }
        flat = self.src.flatten_json(dic)
        self.assertEqual(flat['c_0_d_0'], 2)
        self.assertEqual(flat['c_0_e_0_f'], 1)
        self.assertEqual(flat['h_0'], 1)

