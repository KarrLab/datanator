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

    def test_morphineInChIKey(self):
        key = self.src.inchi_to_inchikey("InChI=1S/C17H19NO3/c1-18-7-6-17-10-3-5-13(20)16(17)21-15-12(19)4-2-9(14(15)17)8-11(10)18/h2-5,10-11,13,16,19-20H,6-8H2,1H3/t10-,11+,13-,16-,17-/m0/s1")
        self.assertEqual(key,'BQJCRHHNABKAKU-KBQPJGBKSA-N')
        key_1 = self.src.inchi_to_inchikey('InChI=1S/H2O/h1H2')
        self.assertEqual(key_1, 'XLYOFNOQVPJJNP-UHFFFAOYSA-N')