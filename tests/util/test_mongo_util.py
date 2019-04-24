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
        cls.collection = 'ecmdb'
        cls.src = mongo_util.MongoUtil(
            cls.cache_dirname, cls.MongoDB, cls.db, cls.collection, verbose=True, max_entries=20)
        cls.client, cls.db, cls.collection_obj = cls.src.con_db(cls.collection)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.client.close()

    def test_con_db(self):
        self.assertNotEqual(self.db, 'Server not available')

    def test_doc_feeder(self):
    	query = {'m2m_id': {'$in': ["M2MDB000004", "M2MDB000005", "M2MDB000006"]}}
    	col = self.src.doc_feeder(self.collection_obj, query=query)
    	for doc in col:
    		if doc['m2m_id'] == "M2MDB000005":
    			self.assertEqual(doc['accession'], "ECMDB00023")
    		elif doc['accession'] == "ECMDB00019":
    			self.assertEqual(doc['m2m_id'], "M2MDB000004")
