import unittest
from datanator.data_source import ec
import datanator.config.core


class TestEC(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        db = 'test'
        username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
        cls.src = ec.EC(server=MongoDB, db=db, username=username, password=password, authSource='admin',
                        readPreference='nearest', max_entries=20)

    @classmethod
    def tearDownClass(cls):
        cls.src.db.drop_collection(cls.src.collection_str)
        cls.src.client.close()