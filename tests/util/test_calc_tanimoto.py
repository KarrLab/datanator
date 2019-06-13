import unittest
from datanator.util import calc_tanimoto
import tempfile
import shutil
from datanator.util import server_util


class TestCalcTanimoto(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'datanator'
        config_file = '/root/host/karr_lab/datanator/.config/config.ini'
        username, password, server, port = server_util.ServerUtil(
            config_file=config_file).get_user_config()
        cls.src = calc_tanimoto.CalcTanimoto(
            cache_dirname=cls.cache_dirname, MongoDB=server, replicaSet=None, db=cls.db,
            verbose=True, max_entries=20, password=password, username=username)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    def test_get_tanimoto(self):
        mol1 = 'InChI=1S/C17H21N4O9P/c1-7-3-9-10(4-8(7)2)21(15-13(18-9)16(25)20-17(26)19-15)5-11(22)14(24)12(23)6-30-31(27,28)29'
        mol2 = 'InChI=1S/C10H7NO3/c12-9(10(13)14)7-5-11-8-4-2-1-3-6(7)8/h1-5,11H,(H,13,14)'
        coe = self.src.get_tanimoto(mol1, mol2)
        self.assertEqual(0.121, coe)
