import unittest
from datanator.data_source import protein_aggregate
import tempfile
import shutil
import json
import os
import datanator.config.core


class TestProteinAggregate(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        src_db = 'datanator'
        des_db = 'test'
        collection_str = 'protein_test'
        username = datanator.config.core.get_config()[
            'datanator']['mongodb']['user']
        password = datanator.config.core.get_config(
        )['datanator']['mongodb']['password']
        server = datanator.config.core.get_config(
        )['datanator']['mongodb']['server']
        port = datanator.config.core.get_config(
        )['datanator']['mongodb']['port']        
        cls.src = protein_aggregate.ProteinAggregate(username=username, password=password, server=server, 
                                                authSource='admin', src_database=src_db, max_entries=20, 
                                                verbose=True, collection=collection_str, destination_database=des_db)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.src.client.close()

    @unittest.skip('passed')
    def test_load_abundance(self):
        self.src.load_abundance()
        entrez_id = self.src.col.find_one(filter={'uniprot_id': 'Q62433'})['entrez_id']
        self.assertEqual(entrez_id, '17988')

    def test_load_ko(self):
        self.src.col.update_one({'uniprot_id': 'a_mock_value'},
            {'$set': {'gene_name': 'gdh'} }, upsert=True) #insert a mock document
        self.src.load_ko()
        doc = self.src.col.find_one(filter={'uniprot_id': 'a_mock_value'})
        self.assertTrue('ko_number' in doc.keys())

    def test_load_taxon(self):
        self.src.load_taxon()
        doc = self.src.col.find_one(filter={'uniprot_id': 'O88425'})
        self.assertTrue('ancestor_name' in doc.keys())
