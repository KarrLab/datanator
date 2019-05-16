from datanator.data_source import intact_nosql
import shutil
import unittest
import tempfile


class TestCorumNoSQL(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.db = 'test'
        self.MongoDB = 'mongodb://mongo:27017/'
        self.src = intact_nosql.IntActNoSQL(
            self.cache_dirname, self.MongoDB, self.db, None, verbose=True, max_entries=20)
        self.src.load_content()

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)
        self.src.client_interaction.close()
        self.src.client_complex.close()

    # @unittest.skip("loading everything")
    def test_load_complex(self):
        int_complex = self.src.collection_complex
        self.assertEqual(int_complex.find().count(), 20)
        cursor = int_complex.find({'identifier': 'CPX-1928'})
        self.assertEqual(cursor.count(), 1)
        self.assertEqual(cursor[0]['ncbi_id'], 83333)
        self.assertEqual(cursor[0]['subunits'], [{'uniprot_id': 'P0AGE0', 'count': '4'}])

    # @unittest.skip('loaded')
    def test_load_interaction(self):
        int_int = self.src.collection_interaction
        self.assertEqual(int_int.count(), 20)
        cursor = int_int.find({'interaction_id': 'intact:EBI-526288'})
        self.assertEqual(cursor.count(), 1)
        self.assertEqual(cursor[0]['method'], 'anti tag coimmunoprecipitation')
        self.assertEqual(cursor[0]['confidence'], 'intact-miscore:0.51')
