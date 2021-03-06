import unittest
from datanator.data_source.metabolite_concentration import doi_10_1038_nchembio_2077
from datanator_query_python.config import config


class TestMetaboliteConcentration(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        des_db = 'test'
        collection_str = 'metabolite_concentrations'
        conf = config.TestConfig()
        username = conf.USERNAME
        password = conf.PASSWORD
        MongoDB = conf.SERVER    
        cls.src = doi_10_1038_nchembio_2077.Concentration(MongoDB=MongoDB, db=des_db, collection_str=collection_str, 
        username=username, password=password, authSource='admin', readPreference='nearest',
        verbose=True, max_entries=float('inf'))

    @classmethod
    def tearDownClass(cls):
        cls.src.client.close()

    def test_flatten_conc_obj(self):
        _input = {'a': 0, 'b': 1, 'c': 2, "internal": "false"}
        result = self.src._flatten_conc_obj(_input, 1000, "XXXXX",'name')
        self.assertEqual(result, [{'a': 0, 'b': 1, 'c': 2,'internal': 'false', 'ncbi_taxonomy_id': 1000,'species_name': 'name', 'species_name': 'XXXXX'}])
        _input = {'a': [0, 1, 2, 3], 'b': [0, 10, 20, 30], 'c': [0, 100, 200, 300], "internal": "true"}
        result = self.src._flatten_conc_obj(_input, 1000, "YYYYY", 'name')
        self.assertEqual(len(result), 4)
        self.assertEqual(result[1], {'a': 1, 'b': 10, 'c': 100, 'internal': 'r', 'ncbi_taxonomy_id': 1000,'species_name': 'name','species_name': 'YYYYY'})