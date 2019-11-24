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
        cls.collection_str = 'test_protein_aggregate'
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
        cls.src.db.drop_collection(cls.collection_str)
        cls.src.client.close()

    def test_get_id(self):
        input_0 = {'resource': [{'namespace': 'something'}, {'id': '2'}, {'namespace': 'sabiork.reaction', 'id': '6570'}]}
        result_0 = self.src.get_rxn_id(input_0)
        self.assertEqual(result_0, 6570)

    def test_create_reactants(self):
        input_0 = {'reaction_participant': [{}, {}, {}, {'substrate_aggregate': '123'}, {'product_aggregate': '456'}]}
        result_0 = self.src.create_reactants(input_0)
        self.assertEqual(result_0, {'substrate_aggregate': '123', 'product_aggregate': '456'})