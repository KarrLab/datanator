import unittest
from datanator.util import calc_tanimoto
import tempfile
import shutil
import datanator.config.core
from datanator.util import mongo_util


class TestCalcTanimoto(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'datanator'
        cls.username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        cls.password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        cls.server = datanator.config.core.get_config()['datanator']['mongodb']['server']
        port = datanator.config.core.get_config()['datanator']['mongodb']['port']
        replSet = datanator.config.core.get_config()['datanator']['mongodb']['replSet']
        cls.src = calc_tanimoto.CalcTanimoto(
            cache_dirname=cls.cache_dirname, MongoDB=cls.server, replicaSet=replSet, db=cls.db,
            verbose=True, max_entries=5, password=cls.password, username=cls.username)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    # @unittest.skip('passed')
    def test_get_tanimoto(self):
        mol1 = 'InChI=1S/C17H21N4O9P/c1-7-3-9-10(4-8(7)2)21(15-13(18-9)16(25)20-17(26)19-15)5-11(22)14(24)12(23)6-30-31(27,28)29'
        mol2 = 'InChI=1S/C10H7NO3/c12-9(10(13)14)7-5-11-8-4-2-1-3-6(7)8/h1-5,11H,(H,13,14)'
        mol3 = 'InChI=1S/C17H21N4O9P/c1-7-3-9-10(4-8(7)2)21(15-13(18-9)16(25)20-17(26)19-15)5-11(22)14(24)12(23)6-30-31(27,28)29'
        coe = self.src.get_tanimoto(mol1, mol2)
        coe2 = self.src.get_tanimoto(mol1, mol3)
        self.assertEqual(0.121, coe)
        self.assertEqual(1., coe2)

    # @unittest.skip('passed')
    def test_one_to_many(self):
        inchi = 'InChI=1S/C6H8O6/c7-1-2(8)5-3(9)4(10)6(11)12-5'
        coeff, hashes = self.src.one_to_many(inchi)
        client, _, col = mongo_util.MongoUtil(db = self.db, MongoDB = self.server,
                                        username = self.username, password = self.password).con_db('metabolites_meta')
        inchi1 = col.find_one({'inchi_hashed': hashes[5]})['inchi_deprot']
        inchi2 = col.find_one({'inchi_hashed': hashes[10]})['inchi_deprot']
        self.assertEqual(coeff[5], self.src.get_tanimoto(inchi, inchi1))
        self.assertEqual(coeff[10], self.src.get_tanimoto(inchi, inchi2))
        client.close()

    # @unittest.skip('passed')
    def test_many_to_many(self):
        client, _, col = mongo_util.MongoUtil(db = self.db, MongoDB = self.server,
                                        username = self.username, password = self.password).con_db('metabolites_meta')
        self.src.many_to_many()
