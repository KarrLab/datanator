import unittest
from datanator.data_source.metabolite_concentration import query_demo
from datanator_query_python.config import config


class TestQueryDemo(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        conf = config.Justin()
        cls.src = query_demo.QueryDemo(MongoDB=conf.SERVER,
                                       db='datanator-test',
                                       collection_str="taxon_tree",
                                       username=conf.USERNAME,
                                       password=conf.PASSWORD)

    @classmethod
    def tearDownClass(cls):
        cls.src.client.close()

    def test_get_canon_ancestors(self):
        tax_id = 1280
        result = self.src.get_canon_ancestors(tax_id)
        self.assertEqual(result[0], {'ncbi_taxonomy_id': 131567, 'name': 'cellular organisms'})
        tax_id = 0
        result = self.src.get_canon_ancestors(tax_id)
        self.assertEqual(result, [])