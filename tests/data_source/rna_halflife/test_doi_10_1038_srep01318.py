import unittest
from datanator.data_source.rna_halflife import doi_10_1038_srep01318
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
        username = conf.USERNAME
        password = conf.PASSWORD
        MongoDB = conf.SERVER    
        cls.src = doi_10_1038_srep01318.Halflife(server=MongoDB, src_db=src_db,
        protein_col=cls.protein_col, authDB='admin', readPreference='nearest',
        username=username, password=password, verbose=True, max_entries=20,
        des_db=des_db, rna_col=cls.rna_col)

    @classmethod
    def tearDownClass(cls):
        # cls.src.uniprot_collection_manager.db.drop_collection(cls.protein_col)
        cls.src.db_obj.drop_collection(cls.rna_col)
        cls.src.uniprot_collection_manager.client.close()
        cls.src.client.close()
        cls.src.uniprot_query_manager.client.close()

    def test_fill_rna_halflife(self):
        d = {'probeset_id': 34555, 'gene_symbol': 'test_symbol', 'gm07029_a1': 1.1236, 'accession_id': 'NM_001762'}
        df_0 = pd.DataFrame(d, index=[0])
        self.src.fill_rna_halflife(df_0, start=0)
