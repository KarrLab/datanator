import tempfile
import shutil
import unittest
from kinetic_datanator.data_source import ensembl



class TestGetGenes(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.gene = ensembl.GetGenes(cache_dirname=self.cache_dirname, download_backups=False, load_content=True)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)

    def test(self):
        pass
