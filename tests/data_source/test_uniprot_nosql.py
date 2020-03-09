import unittest
import shutil
import tempfile
from datanator.data_source import uniprot_nosql
import datanator.config.core


class TestUniprotNoSQL(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        db = 'test'
        username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
        cls.src = uniprot_nosql.UniprotNoSQL(MongoDB=MongoDB, db=db, max_entries=20,
                                            username=username, password=password, collection_str='test_uniprot')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.src.db.drop_collection(cls.src.collection_str)

    # @unittest.skip('large single file download')
    def test_proper_loading(self):
        self.src.load_uniprot()
        # count = uni.count()
        # self.assertEqual(count, 10)
        # self.assertNotEqual(uni.find_one()['gene_name'], None)

    def test_fill_species_name(self):
        self.src.fill_species_name()