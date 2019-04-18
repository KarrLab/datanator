import unittest
from datanator.data_source import pax_nosql
import tempfile
import shutil

class TestCorumNoSQL(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        # '../../datanator/data_source/cache'
        self.db = 'test'
        self.MongoDB = 'mongodb://mongo:27017/'

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_con_db(self):
        src = pax_nosql.PaxNoSQL(
            self.cache_dirname, self.MongoDB, self.db, verbose=True, max_entries=20)
        collection = src.con_db()
        self.assertNotEqual(collection, 'Server not available')
        collection.drop()

    # only loads partial content because it takes too long to load everything
    def test_load_content(self):
        src = pax_nosql.PaxNoSQL(
            self.cache_dirname, self.MongoDB, self.db, verbose=True, max_entries = 5)
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
        collection.drop()
