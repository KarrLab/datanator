import unittest
import shutil
import tempfile
from datanator.data_source import uniprot_nosql
import os
import json

class TestUniprotNoSQL(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.MongoDB = 'mongodb://mongo:27017/'
        cls.db = 'tests'
        cls.output_directory = cls.cache_dirname # directory to store JSON files
        cls.src = uniprot_nosql.UniprotNoSQL(cls.MongoDB, cls.db, max_entries=10)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    def setUp(self):
        self.collection = self.src.con_db()

    def tearDown(self):
        self.collection.drop()

    def test_con_db(self):
        self.assertNotEqual(self.collection, 'Server not available')

    def test_proper_loading(self):
        uni = self.src.load_uniprot()
        count = uni.count()
        self.assertEqual(count, 10)
        self.assertEqual(count[0]['uniprot_id'], 'Q2KA61')
