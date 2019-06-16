import unittest
from datanator.util import server_util

class TestServerUtil(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.config_file = '/root/host/karr_lab/datanator/.config/config_test.ini'
        cls.username = None
        cls.password = None
        cls.port = None
        cls.src = server_util.ServerUtil(
            config_file=cls.config_file, username = cls.username, password = cls.password,
    			port = cls.port)

    # @classmethod
    # def tearDownClass(cls):
    #     shutil.rmtree(cls.cache_dirname)

    def test_get_admin_config(self):
    	username, password, server, port = self.src.get_user_config()
    	self.assertEqual(username, 'some-user')
    	self.assertEqual(password, 'some-password')
    	self.assertEqual(server, '35.173.159.185')
    	self.assertEqual(port, '27017')