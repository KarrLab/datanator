import unittest
from datanator.data_source.rna_halflife import doi_10_1186_gb_2012_13_4_r30
import tempfile
import shutil
import json
import os
from datanator_query_python.config import config
import pandas as pd


class TestProteinAggregate(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        des_db = 'test'
        src_db = 'datanator'
        cls.protein_col = 'uniprot'
        cls.rna_col = 'rna_halflife'
        conf = config.TestConfig()
        username = conf.MONGO_TEST_USERNAME
        password = conf.MONGO_TEST_PASSWORD
        MongoDB = conf.SERVER    
        cls.src = doi_10_1186_gb_2012_13_4_r30.Halflife(server=MongoDB, src_db=src_db,
        protein_col=cls.protein_col, authDB='admin', readPreference='nearest',
        username=username, password=password, verbose=True, max_entries=20,
        des_db=des_db, rna_col=cls.rna_col)

    @classmethod
    def tearDownClass(cls):
        cls.src.uniprot_collection_manager.db.drop_collection(cls.protein_col)
        cls.src.db_obj.drop_collection(cls.rna_col)
        cls.src.uniprot_collection_manager.client.close()
        cls.src.client.close()
        cls.src.uniprot_query_manager.client.close()

    @unittest.skip('passed')
    def test_load_uniprot(self):
        self.src.load_uniprot()

    def test_fill_rna_half_life(self):
        url = """https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2012-13-4-r30/MediaObjects/13059_2011_2880_MOESM3_ESM.XLSX"""
        names = ['ordered_locus_name', 'half_life_ga_2', 'reads_per_kb_per_mb',
                'transcriptional_start_sites', 'transcriptional_end_sites', 'operon',
                'gene_start', 'gene_end', 'strand', 'gene_name', 'protein_annotation',
                'cog', 'kegg', 'half_life_qpcr', 'half_life_454']
        df_10987 = self.src.make_df(url, 'Bc10987', names=names, usecols='A:O', skiprows=[0,1], nrows=self.src.max_entries)
        self.src.fill_rna_half_life(df_10987, names, ['Bacillus cereus ATCC 10987', 222523])
        self.src.fill_rna_half_life(df_10987, names, ['Bacillus cereus ATCC 10987', 222523], quantification_method='RT-qPCR')
        self.src.fill_rna_half_life(df_10987, names, ['Bacillus cereus ATCC 10987', 222523], quantification_method='Roche 454')