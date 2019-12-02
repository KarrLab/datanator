import unittest
from datanator.data_source.rna_halflife import doi_10_1186_s12864_016_3219_8
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
        des_db = 'test'
        cls.collection_str = 'test_rna_halflife'
        username = datanator.config.core.get_config()[
            'datanator']['mongodb']['user']
        password = datanator.config.core.get_config(
        )['datanator']['mongodb']['password']
        server = datanator.config.core.get_config(
        )['datanator']['mongodb']['server']       
        cls.src = doi_10_1186_s12864_016_3219_8.Halflife(username=username, password=password, server=server, 
                                                authDB='admin',max_entries=100, uniprot_col_db=des_db,
                                                verbose=True, collection_str=cls.collection_str, db=des_db,
                                                cache_dir=cache_dir)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        # cls.src.db.drop_collection(cls.collection_str)
        # cls.src.uniprot_col_manager.db.drop_collection('uniprot')
        cls.src.client.close()

    def test_download_xlsx(self):
        result = self.src.download_xlsx('MeOH')
        self.assertEqual(result['gene_fragment'][0], 'MA0001')

    @unittest.skip('passed')
    def test_load_halflife(self):
        df = self.src.download_xlsx('MeOH')
        self.src.load_halflife(df)
        df = self.src.download_xlsx('TMA')
        self.src.add_to_halflife(df)

    @unittest.skip('passed')
    def test_fill_gene_protein_name(self):
        self.src.fill_gene_protein_name()
        result = self.src.collection.find_one({'gene_name': '-'})
        self.assertIsNone(result)

    @unittest.skip('passed')
    def test_fill_protein_name(self):
        self.src.fill_protein_name()
        result = self.src.collection.find_one({'$and':[{'gene_name': {'$ne': '-'}}, 
                                                   {'protein_name': {'$exists': False}},
                                                   {'gene_name': {'$exists': True}}]})
        self.assertIsNone(result)

    def test_fill_uniprot_by_oln(self):
        self.src.fill_uniprot_by_oln('MA0002')