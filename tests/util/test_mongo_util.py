import unittest
from datanator.util import mongo_util
import tempfile
import shutil


class TestMongoUtil(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'test'
        cls.MongoDB = 'mongodb://mongo:27017/'
        cls.src = mongo_util.MongoUtil(
            cls.cache_dirname, cls.MongoDB, None, cls.db, verbose=True, max_entries=20)
        cls.collection_str = 'ecmdb'
        cls.client, cls.db, cls.collection_obj = cls.src.con_db(
            cls.collection_str)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        # cls.client.drop_database(cls.db)
        cls.client.close()

    # @unittest.skip('passed')
    def test_con_db(self):
        self.assertNotEqual(self.db, 'Server not available')

    # @unittest.skip('passed')
    def test_fill_db(self):
        self.collection_obj.drop()
        self.assertEqual(self.collection_obj.find().count(), 0)
        collection_obj = self.src.fill_db(self.collection_str)
        self.assertNotEqual(collection_obj.find().count(), 0)

