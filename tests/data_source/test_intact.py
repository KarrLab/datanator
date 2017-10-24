import unittest
from kinetic_datanator.data_source import intact
import tempfile
import shutil


class TestIntAct(unittest.TestCase):
    """

    """

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.intact = intact.IntAct(cache_dirname = self.cache_dirname,
                                clear_content = False,
                                download_backup= True, load_content = False,
                                verbose = True)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)

    def test_loading(self):

        q = self.intact.session.query(intact.ProteinInteractions).filter_by(interactor_a = 'uniprotkb:P27986').count()
        self.assertEqual(q, 274)

        q = self.intact.session.query(intact.ProteinInteractions).filter_by(interactor_a = 'uniprotkb:Q61824').first()
        self.assertEqual(q.interactor_b, 'uniprotkb:Q60631')
        self.assertEqual(q.publications, 'pubmed:11127814|mint:MINT-5213342')

    #TODO: Add more Tests
