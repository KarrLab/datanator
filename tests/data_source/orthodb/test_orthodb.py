from datanator.data_source.orthodb import orthodb
from datanator_query_python.config import config
import unittest
import numpy as np


class TestTransform(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        conf = config.SchemaMigration()
        cls.des_col = "orthodb"
        cls.src = orthodb.OrthoDB(MongoDB=conf.SERVER,
                                    db="test",
                                    des_col=cls.des_col,
                                    username=conf.USERNAME,
                                    password=conf.PASSWORD,
                                    max_entries=10,
                                    verbose=True)

    @classmethod
    def tearDownClass(cls):
        cls.src.db_obj.drop_collection(cls.des_col)

    def test_pairwise_name_group(self):
        url = './docs/orthodb/odb10v1_OGs.tab'
        self.src.pairwise_name_group(url)