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
    def test_print_dict(self):
        example_dict = { 'key1' : 'value1',
                         'key2' : 'value2',
                         'key3' : { 'key3a': 'value3a' },
                         'key4' : { 'key4a': { 'key4aa': 'value4aa',
                                               'key4ab': 'value4ab',
                                               'key4ac': 'value4ac'},
                                    'key4b': 'value4b'}
                       }
        self.src.print_dict(example_dict)

    def test_print_schema(self):
        self.src.print_dict('ecmdb')
