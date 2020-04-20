import unittest
import shutil
import tempfile
from datanator.data_source.brenda import reaction
from datanator_query_python.config import config


class TestBrendaRxn(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        conf = config.TestConfig()
        cls.collection_str = 'brenda_reaction'
        username = conf.MONGO_TEST_USERNAME
        password = conf.MONGO_TEST_PASSWORD
        MongoDB = conf.SERVER
        cls.src = reaction.BrendaRxn(MongoDB=MongoDB, db='test', collection_str=cls.collection_str,
                                     username=username, password=password, authSource='admin',
                                     max_entries=20, verbose=True)

    @classmethod
    def tearDownClass(cls):
        cls.src.db_obj.drop_collection(cls.collection_str)
        cls.src.client.close()

    def test_download_and_read(self):
        result = self.src.download_and_read()
        print(result)