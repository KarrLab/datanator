import unittest
from datanator.core import query_nosql
import tempfile
import shutil
import configparser
import os

class TestQueryNoSQL(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'datanator'
        parser = configparser.ConfigParser(allow_no_value=True)
        parser.read(os.path.expanduser('~/.wc/datanator.ini'))
        username = parser.get('mongodb', 'user')
        password = parser.get('mongodb', 'password')
        MongoDB = parser.get('mongodb', 'server')
        port = int(parser.get('mongodb', 'port'))
        replSet = parser.get('mongodb', 'replSet')
        cls.collection_str = 'ecmdb'
        cls.src = query_nosql.DataQuery(
            cache_dirname=cls.cache_dirname, MongoDB=MongoDB, replicaSet=None, db=cls.db,
                 verbose=True, max_entries=20, username = username, password = password)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    # @unittest.skip('skip to testing for h1_hesc')
    def test_doc_feeder(self):
        query = {'m2m_id': {
            '$in': ["M2MDB000004", "M2MDB000005", "M2MDB000006"]}}
        col = self.src.doc_feeder(self.collection_str, query=query)
        for doc in col:
            if doc['m2m_id'] == "M2MDB000005":
                self.assertEqual(doc['accession'], "ECMDB00023")
            elif doc['accession'] == "ECMDB00019":
                self.assertEqual(doc['m2m_id'], "M2MDB000004")