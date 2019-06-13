import unittest
from datanator.core import mongo_util
import tempfile
import shutil
from datanator.util import server_util

class TestCalcTanimoto(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'datanator'
        config_file = '/root/host/karr_lab/datanator/.config/config.ini'
        username, password, server, port = server_util.ServerUtil(config_file = config_file).get_user_config()
        cls.src = mongo_util.MongoUtil(
            cache_dirname=cls.cache_dirname, MongoDB=server, replicaSet=None, db=cls.db,
                 verbose=True, max_entries=20, password = password, username = username)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    def test_get_tanimoto(self):
    	pass