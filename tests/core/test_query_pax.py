import unittest
from datanator.core import query_pax
from datanator.util import file_util
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
        cls.file_manager = file_util.FileUtil()
        cls.src = query_pax.QueryPax(
            cache_dirname=cls.cache_dirname, MongoDB=cls.MongoDB, db=cls.db,
                 verbose=True, max_entries=20, username = cls.username, password = cls.password)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    def test_get_all_species(self):
        result = self.src.get_all_species()
        self.assertTrue('Synechocystis.sp. 6803' in result)

    def test_get_abundance_from_uniprot(self):
        uniprot_id = 'F4KDK1'
        result = self.src.get_abundance_from_uniprot(uniprot_id)
        dic = self.file_manager.search_dict_list(result, 'abundance', '15.2')
        exp = [{'organ': 'COTYLEDON', 'abundance': '15.2'}]
        self.assertEqual(dic, exp)
        self.assertEqual({'ncbi_taxonomy_id': 3702, 'species_name': 'A.thaliana'}, result[0] )
        uniprot_id_1 = 'asdfasdf'
        result_1 = self.src.get_abundance_from_uniprot(uniprot_id_1)
        self.assertEqual(result_1, [])