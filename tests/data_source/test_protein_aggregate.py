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
        cls.src = protein_aggregate.ProteinAggregate(username=username, password=password, server=server, 
                                                authSource='admin', src_database=src_db, max_entries=20, 
                                                verbose=True, collection=cls.collection_str, destination_database=des_db,
                                                cache_dir=cache_dir)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        # cls.src.db.drop_collection(cls.collection_str)
        cls.src.client.close()

    # # @unittest.skip('passed')
    # def test_load_abundance_from_pax(self):
    #     self.src.load_abundance_from_pax()
    #     doc = self.src.col.find_one(filter={'uniprot_id': 'Q72DI0'})
    #     self.assertTrue('abundances' in doc.keys())
    #     self.assertTrue('ncbi_taxonomy_id' in doc.keys())

    # # @unittest.skip('passed')
    # def test_load_ko(self):
    #     self.src.col.insert_one({'uniprot_id': 'a_mock_value',
    #         'gene_name': 'gdh'}) #insert a mock document
    #     self.src.load_ko()
    #     doc = self.src.col.find_one(filter={'uniprot_id': 'a_mock_value'})
    #     self.assertTrue('ko_number' in doc.keys())

    # # @unittest.skip('passed')
    # def test_load_taxon(self):
    #     self.src.col.insert_one({'ncbi_taxonomy_id': 9606,
    #         'uniprot_id': 'taxon_mock_value'})
    #     self.src.load_taxon()
    #     doc = self.src.col.find_one(filter={'uniprot_id': 'taxon_mock_value'})
    #     self.assertTrue('ancestor_name' in doc.keys())

    # def test_load_unreviewed_abundance(self):
    #     dic_0 = {'observation': [{'protein_id': {'string_id': 'string_mock_0', 'uniprot_id': 'id_mock_0'},
    #     'string_id': 'string_mock_0', 'abundance': 0 }], 'ncbi_id': 0, 'species_name': 'name_mock_0', 'organ': 'organ_0'}
    #     dic_1 = {'observation': [{'protein_id': {'string_id': 'string_mock_1', 'uniprot_id': 'id_mock_1'},
    #     'string_id': 'string_mock_1', 'abundance': 1 }], 'ncbi_id': 1, 'species_name': 'name_mock_1', 'organ': 'organ_1'}
    #     dic_2 = {'uniprot_id': 'id_mock_0', 'abundances': []}
    #     dic_3 = {'uniprot_id': 'Q72DIO'}
    #     self.src.col.insert_many([dic_2, dic_3])
    #     self.src.load_unreviewed_abundance()
    #     doc = self.src.col.find_one(filter={'species_name': 'D.vulgaris'})
    #     self.assertEqual(doc['ncbi_taxonomy_id'], 882)

    @unittest.skip('removed the function')
    def test_loadload_kinlaw_from_sabio(self):
        dic_0 = {'uniprot_id': 'P20932'}
        dic_1 = {'uniprot_id': 'id_mock_1', 'protein_name': 'subtilisin'}
        dic_2 = {'uniprot_id': 'P16064', 'protein_name': 'subtilisin'}
        self.src.col.insert_many([dic_0, dic_1, dic_2])
        self.src.load_kinlaw_from_sabio()
        result_0 = self.src.col.find_one({'uniprot_id': 'P20932'})
        self.assertTrue('kinetics' in list(result_0.keys()))
        self.assertTrue({'kinlaw_id': 17, 'ncbi_taxonomy_id': 303} in result_0['kinetics'])
        result_1 = self.src.col.find_one({'uniprot_id': 'P16064'})
        self.assertTrue({'kinlaw_id': 1, 'ncbi_taxonomy_id': 1467} in result_1['kinetics'])

