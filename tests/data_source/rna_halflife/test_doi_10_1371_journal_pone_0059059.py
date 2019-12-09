import unittest
from datanator.data_source.rna_halflife import doi_10_1371_journal_pone_0059059
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
        cls.cache_dir = tempfile.mkdtemp()
        conf = config.TestConfig()
        username = conf.MONGO_TEST_USERNAME
        password = conf.MONGO_TEST_PASSWORD
        MongoDB = conf.SERVER    
        cls.src = doi_10_1371_journal_pone_0059059.Halflife(server=MongoDB, src_db=src_db,
        protein_col=cls.protein_col, authDB='admin', readPreference='nearest',
        username=username, password=password, verbose=True, max_entries=20,
        des_db=des_db, rna_col=cls.rna_col, cache_dir=cls.cache_dir)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dir)
        cls.src.uniprot_collection_manager.db.drop_collection(cls.protein_col)
        cls.src.db_obj.drop_collection(cls.rna_col)
        cls.src.uniprot_collection_manager.client.close()
        cls.src.client.close()
        cls.src.uniprot_query_manager.client.close()