from datanator.data_source import intact_nosql
import shutil
import unittest
import tempfile
import datanator.config.core

class TestCorumNoSQL(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.db = 'test'
        self.username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        self.password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        self.MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
        port = datanator.config.core.get_config()['datanator']['mongodb']['port']
        self.replSet = datanator.config.core.get_config()['datanator']['mongodb']['replSet']
        self.src = intact_nosql.IntActNoSQL(
            cache_dirname = self.cache_dirname, MongoDB = self.MongoDB, 
            db = self.db, verbose=True, max_entries=20,
            username  = self.username, password = self.password)
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
