from datanator.data_source import sabio_rk_trial
import unittest
import tempfile
import shutil

class TestSabioRk(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.src = sabio_rk_trial.SabioRk()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    def test_suds_trial(self):
    	self.src.suds_trial()