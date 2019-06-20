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
        port = datanator.config.core.get_config()['datanator']['mongodb']['port']
        replSet = datanator.config.core.get_config()['datanator']['mongodb']['replSet']
        cls.src = uniprot_nosql.UniprotNoSQL(MongoDB = MongoDB, db = db, max_entries=10,
                                            username = username, password = password)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)


    @unittest.skip('large single file download')
    def test_proper_loading(self):
        uni = self.src.load_uniprot()
        count = uni.count()
        self.assertEqual(count, 10)
        self.assertNotEqual(uni.find_one()['gene_name'], None)
