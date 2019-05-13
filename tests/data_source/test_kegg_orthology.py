import unittest
from datanator.data_source import kegg_orthology
import tempfile
import shutil
import pymongo


class TestKeggOrthology(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.db = 'test'
        self.MongoDB = 'mongodb://mongo:27017/'

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_con_db(self):
        src = kegg_orthology.KeggOrthology(
            self.cache_dirname, self.MongoDB, self.db, replicaSet=None, verbose=True, max_entries=20)
        client, db, collection = src.con_db('kegg_orthology')
        self.assertNotEqual(collection, 'Server not available')
        client.close()

    def test_get_json_ends(self):
        pass

    def test_parse_ko_txt(self):
        src = kegg_orthology.KeggOrthology(
            self.cache_dirname, self.MongoDB, self.db, replicaSet=None, verbose=True, max_entries=20)
        client, db, collection = src.con_db('kegg_orthology')



        client.close()



