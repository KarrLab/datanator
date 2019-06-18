import unittest
from datanator.util import chem_util
import tempfile
import shutil


class TestChemUtil(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.src = chem_util.ChemUtil()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    def test_simplify_inchi(self):
        inchi = 'InChI=1S/H2O/h1H2'
        self.assertEqual('InChI=1S/H2O', self.src.simplify_inchi(inchi))
        inchi = None
        self.assertEqual('InChI = None', self.src.simplify_inchi(inchi))

    def test_hash_inchi(self):
    	inchi = 'InChI=1S/C6H12N2O4S2/c7-3(5(9)10)1-13-14-2-4(8)6(11)12'
    	hashed = 'e0a402c94a0ecd52ec426756854592f76eece8fd3ffef2e7347fb6c5'
    	self.assertEqual(hashed, self.src.hash_inchi(inchi))
    	self.assertEqual('InChI = None', self.src.hash_inchi(None))