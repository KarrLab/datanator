import unittest
from datanator.query import query_kegg_orthology
from datanator.util import file_util
import tempfile
import shutil
import configparser
import os

class TestQueryKO(unittest.TestCase):

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
        cls.MongoDB = MongoDB
        cls.username = username
        cls.password = password
        cls.src = query_kegg_orthology.QueryKO(server=cls.MongoDB, database=cls.db,
                 verbose=True, max_entries=20, username = cls.username, password = cls.password)
        cls.file_manager = file_util.FileUtil()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.src.client.close()

    def test_get_ko_by_name(self):
        result_0 = self.src.get_ko_by_name('gyar')
        self.assertEqual('K00015', result_0)
        result_1 = self.src.get_ko_by_name('gyaR')
        self.assertEqual(result_1, result_0)
        result_2 = self.src.get_ko_by_name('yuyyyyyy')
        self.assertEqual(None, result_2)