import tempfile
import shutil
import unittest
from kinetic_datanator.data_source import ensembl



class TestGetGenes(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.gene = ensembl.GetGenes(cache_dirname=cls.cache_dirname, download_backups=False, load_content=True)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    def test(self):
        pass
