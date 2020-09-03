from datanator.schema_2 import add_orthodb
from datanator_query_python.config import config
import unittest


class TestTransform(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        conf = config.SchemaMigration()
        cls.des_col = "uniprot"
        cls.src = add_orthodb.AddOrtho(MongoDB=conf.SERVER,
                                                db="test",
                                                des_col=cls.des_col,
                                                username=conf.USERNAME,
                                                password=conf.PASSWORD,
                                                max_entries=10,
                                                verbose=True)

    @classmethod
    def tearDownClass(cls):
        pass

    def test_add_ortho(self):
        self.src.add_ortho()