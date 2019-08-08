import unittest
from datanator_flask_REST.query import query_pax
import tempfile
import shutil
import configparser
import os

class TestQueryPax(unittest.TestCase):

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
        cls.src = query_pax.QueryPax(
            cache_dirname=cls.cache_dirname, MongoDB=cls.MongoDB, db=cls.db,
                 verbose=True, max_entries=20, username = cls.username, password = cls.password)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    def test_get_all_species(self):
        result = self.src.get_all_species()
        self.assertTrue('Synechocystis.sp. 6803' in result)