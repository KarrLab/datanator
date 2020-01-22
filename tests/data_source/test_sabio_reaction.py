import unittest
from datanator.data_source import sabio_reaction
import tempfile
import shutil
import json
import os
import datanator.config.core


class TestProteinAggregate(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cache_dir = os.path.join(cls.cache_dirname, 'logs.txt')
        src_db = 'datanator'
        des_db = 'test'
        cls.collection_str = 'sabio_reaction'
        username = datanator.config.core.get_config()[
            'datanator']['mongodb']['user']
        password = datanator.config.core.get_config(
        )['datanator']['mongodb']['password']
        server = datanator.config.core.get_config(
        )['datanator']['mongodb']['server']
        port = datanator.config.core.get_config(
        )['datanator']['mongodb']['port']        
        cls.src = sabio_reaction.RxnAggregate(username=username, password=password, server=server, 
                                                authSource='admin', src_database=src_db, max_entries=20, 
                                                verbose=True, collection=cls.collection_str, destination_database=des_db,
                                                cache_dir=cache_dir)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        # cls.src.db.drop_collection(cls.collection_str)
        cls.src.client.close()

    def test_get_id(self):
        input_0 = {'resource': [{'namespace': 'something'}, {'id': '2'}, {'namespace': 'sabiork.reaction', 'id': '6570'}]}
        result_0 = self.src.get_rxn_id(input_0)
        self.assertEqual(result_0, 6570)

    def test_create_reactants(self):
        input_0 = {'reaction_participant': [{}, {}, {}, {'substrate_aggregate': '123'}, {'product_aggregate': '456'}]}
        result_0 = self.src.create_reactants(input_0)
        self.assertEqual(result_0, {'substrate_aggregate': '123', 'product_aggregate': '456'})

    def test_fill_collection(self):
        self.src.fill_collection()

    def test_extract_reactant_names(self):
        substrates_0 = {'substrate_name': 'a', 'substrate_synonym': ['a1', 'a2', 'a3']}
        substrates_1 = {'substrate_name': 'b', 'substrate_synonym': ['b1', 'b2', 'b3']}
        products_0 = {'product_name': 'c', 'product_synonym': ['c1', 'c2', 'c3']}
        products_1 = {'product_name': 'd', 'product_synonym': ['d1', 'd2', 'd3']}
        products_2 = {'product_name': 'e', 'product_synonym': []}
        input_0 = {'reaction_participant': [{'substrate': [substrates_0, substrates_1]},{'product': [products_0, products_1, products_2]}]}
        sub_0, pro_0 = self.src.extract_reactant_names(input_0)
        sub_exp_0 = [['a1', 'a2', 'a3', 'a'], ['b1', 'b2', 'b3', 'b']]
        pro_exp_0 = [['c1', 'c2', 'c3', 'c'], ['d1', 'd2', 'd3', 'd'], ['e']]
        self.assertEqual(sub_0, sub_exp_0)
        self.assertEqual(pro_0, pro_exp_0)

    def test_extract_enzyme_names(self):
        input_0 = {'enzymes': [{'enzyme':[{'enzyme_name': 'a', 'enzyme_synonym': ['a1', 'a2', 'a3']}]}]}
        input_1 = {'enzymes': [{'enzyme':[{'enzyme_name': 'a', 'enzyme_synonym': ['a1', 'a2', 'a3']},
                                          {'enzyme_name': 'b', 'enzyme_synonym': ['b1', 'b2', 'b3']}]}]}
        result_0 = self.src.extract_enzyme_names(input_0)
        result_1 = self.src.extract_enzyme_names(input_1)
        self.assertEqual(result_0, ['a1', 'a2', 'a3', 'a']) 
        self.assertEqual(result_1[0], ['a1', 'a2', 'a3', 'a'])