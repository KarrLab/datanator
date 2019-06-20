import unittest
from datanator.data_source import pax_nosql
import datanator.config.core
import tempfile
import shutil

class TestCorumNoSQL(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.db = 'test'
        self.username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        self.password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        self.MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
        port = datanator.config.core.get_config()['datanator']['mongodb']['port']
        replSet = datanator.config.core.get_config()['datanator']['mongodb']['replSet']

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    # only loads partial content because it takes too long to load everything
    def test_load_content(self):
        src = pax_nosql.PaxNoSQL(
            self.cache_dirname, self.MongoDB, self.db, verbose=True, max_entries = 5,
            password = self.password, username = self.username)
        collection = src.load_content()
        self.assertEqual(collection.count(), 5)
        cursor = collection.find({'file_name': '882/882-WHOLE_ORGANISM-integrated.txt'})
        self.assertEqual(cursor.count(), 1)
        self.assertEqual(cursor[0]['species_name'], 'D.vulgaris')
        self.assertEqual(cursor[0]['observation'][0]['string_id'], '882.DVU0949')
        cursor = collection.find({'file_name': '882/882-Desulfo_Lac_Exp_SC_zhang_2006.txt'})
        self.assertEqual(cursor.count(), 1)
        self.assertEqual(cursor[0]['weight'], 20)
        self.assertEqual(cursor[0]['observation'][1]['string_id'], '882.DVU0142')
