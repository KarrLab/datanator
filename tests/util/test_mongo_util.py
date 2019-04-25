import unittest
from datanator.util import mongo_util
import tempfile
import shutil


class TestMongoUtil(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'test'
        cls.MongoDB = 'mongodb://mongo_secondary:27017/'
        cls.src = mongo_util.MongoUtil(
            cls.cache_dirname, cls.MongoDB, cls.db, verbose=True, max_entries=20)
        cls.collection_str = 'ecmdb'
        cls.client, cls.db, cls.collection_obj = cls.src.con_db(
            cls.collection_str)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.client.drop_database(cls.db)
        cls.client.close()

    def test_con_db(self):
        self.assertNotEqual(self.db, 'Server not available')

    def test_fill_db(self):
        self.collection_obj.drop()
        self.assertEqual(self.collection_obj.find().count(), 0)
        collection_obj = self.src.fill_db(self.collection_str)
        self.assertNotEqual(collection_obj.find().count(), 0)

    def test_doc_feeder(self):
        collection_obj = self.src.fill_db(self.collection_str)
        query = {'m2m_id': {
            '$in': ["M2MDB000004", "M2MDB000005", "M2MDB000006"]}}
        col = self.src.doc_feeder(self.collection_str, query=query)
        for doc in col:
            if doc['m2m_id'] == "M2MDB000005":
                self.assertEqual(doc['accession'], "ECMDB00023")
            elif doc['accession'] == "ECMDB00019":
                self.assertEqual(doc['m2m_id'], "M2MDB000004")
