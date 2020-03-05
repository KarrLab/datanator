from datanator.data_source import intact_nosql
import shutil
import os
import unittest
import tempfile
import datanator.config.core

class TestCorumNoSQL(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'test'
        cls.username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        cls.password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        cls.MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
        cls.replSet = datanator.config.core.get_config()['datanator']['mongodb']['replSet']
        cls.src = intact_nosql.IntActNoSQL(
            cache_dirname = cls.cache_dirname, MongoDB = cls.MongoDB, 
            db = cls.db, verbose=True, max_entries=20,
            username  = cls.username, password = cls.password)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.src.collection_interaction.drop()
        cls.src.collection_complex.drop()
        cls.src.client_interaction.close()
        cls.src.client_complex.close()

    def test_download_content(self):
        self.src.download_content()
        file = os.path.join(self.cache_dirname, 'intact', 'complextab', 'bos_taurus.tsv')
        self.assertTrue(os.path.exists(file))

    # @unittest.skip("loading everything")
    def test_load_complex(self):
        self.src.add_complexes()
        int_complex = self.src.collection_complex
        self.assertTrue(int_complex.find().count() > 11)
        # cursor = int_complex.find({'identifier': 'CPX-3140'})
        # self.assertEqual(cursor.count(), 1)
        # self.assertEqual(cursor[0]['ncbi_id'], 7227)
        # self.assertEqual(cursor[0]['subunits'], [{'uniprot_id': 'P48607-PRO_0000022407', 'count': '1'},
        #                                         {'uniprot_id': 'P48607-PRO_0000022407', 'count': '1'}])

    # @unittest.skip('loaded')
    def test_load_interaction(self):
        self.src.add_interactions()
        int_int = self.src.collection_interaction
        self.assertTrue(int_int.count() in [18, 19, 20])
        cursor = int_int.find({'interaction_id': 'intact:EBI-526288'})
        self.assertEqual(cursor.count(), 1)
        self.assertEqual(cursor[0]['method'], 'anti tag coimmunoprecipitation')
        self.assertEqual(cursor[0]['confidence'], 'intact-miscore:0.51')
