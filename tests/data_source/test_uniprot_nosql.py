import unittest
import shutil
import tempfile
from datanator.data_source import uniprot_nosql
from datanator.util import server_util


class TestUniprotNoSQL(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        db = 'test'
        config_file = '/root/host/karr_lab/datanator/.config/config.ini'
        username, password, MongoDB, port = server_util.ServerUtil(
            config_file=config_file).get_user_config()
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
