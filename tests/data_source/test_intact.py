import unittest
from kinetic_datanator.data_source import intact
import tempfile
import shutil


class TestFromServerIntAct(unittest.TestCase):
    """

    """

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.intact = intact.IntAct(cache_dirname = self.cache_dirname,
                                clear_content = False,
                                download_backups= True, load_content = False,
                                verbose = True)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)


    def test_add_complex(self):
        q = self.intact.session.query(intact.ProteinComplex).get('EBI-1256672')

        self.assertEqual(q.name, 'INO80 chromatin remodeling complex')
        self.assertEqual(q.ncbi, '559292')
        self.assertEqual(q.evidence, 'intact:EBI-515508')


    def test_add_interactions(self):
        q = self.intact.session.query(intact.ProteinInteractions).filter_by(interactor_a = 'uniprotkb:P27986').count()
        self.assertEqual(q, 274)

        q = self.intact.session.query(intact.ProteinInteractions).filter_by(interactor_a = 'uniprotkb:Q61824').first()
        self.assertEqual(q.interactor_b, 'uniprotkb:Q60631')
        self.assertEqual(q.publications, 'pubmed:11127814|mint:MINT-5213342')


class TetstLoadingIntAct(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.intact = intact.IntAct(cache_dirname=self.cache_dirname, download_backups=False, load_content=True, max_entries=10)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)

    def testloadedcomplex(self):
        q = self.intact.session.query(intact.ProteinComplex).get('EBI-1256672')

        self.assertEqual(q.name, 'INO80 chromatin remodeling complex')
        self.assertEqual(q.ncbi, '559292')
        self.assertEqual(q.evidence, 'intact:EBI-515508')


    def testloadedinteractions(self):
        q = self.intact.session.query(intact.ProteinInteractions).all()

        self.assertEqual(len(q), self.intact.max_entries)
        self.assertEqual(q[0].interactor_a, 'intact:EBI-7121510')
