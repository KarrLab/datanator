import unittest
from datanator.util import standardize_util
import tempfile
import shutil
from bson.objectid import ObjectId


class TestStandardizeUtil(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'test'
        cls.MongoDB = 'mongodb://mongo:27017/'
        cls.src = standardize_util.StandardizeUtil(
            cache_dirname= cls.cache_dirname, MongoDB = cls.MongoDB, 
            db = cls.db, verbose=True, max_entries=20)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    @unittest.skip('done')
    def test_standardize_util(self):
        self.src.standardize_sabio()
    
    @unittest.skip('done')
    def test_standardize_util(self):
        self.src.standardize_metabolite()
