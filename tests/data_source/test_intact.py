import unittest
from kinetic_datanator.data_source import intact
import tempfile
import shutil


class TestFromServerIntAct(unittest.TestCase):
    """
    Testing IntAct database from cached server
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
        self.assertEqual(q.go_annot, 'GO:0031011(Ino80 complex)|GO:0043044(ATP-dependent chromatin remodeling)|GO:0016887(ATPase activity)|GO:0043138(3&apos;-5&apos; DNA helicase activity)|GO:0006281(DNA repair)|GO:0045449(regulation of transcription, DNA-templated)')
        self.assertEqual(q.evidence, 'intact:EBI-515508')


    def test_add_interactions(self):
        q = self.intact.session.query(intact.ProteinInteraction).filter_by(protein_a = 'P27986').count()
        self.assertEqual(q, 272)

        q = self.intact.session.query(intact.ProteinInteraction).filter_by(protein_a = 'Q61824').first()
        self.assertEqual(q.protein_b, 'Q60631')
        self.assertEqual(q.gene_a, 'Adam12')
        self.assertEqual(q.gene_b, 'Grb2')
        self.assertEqual(q.type_a, q.type_b)
        self.assertEqual(q.role_a, q.role_b)
        self.assertEqual(q.method, 'coimmunoprecipitation')
        self.assertEqual(q.publication, '11127814')
        self.assertEqual(q.publication_author, 'Suzuki et al. (2000)')


class TestLoadingIntAct(unittest.TestCase):
    """
    Testing loading IntAct database
    """

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.intact = intact.IntAct(cache_dirname=self.cache_dirname , download_backups=False, load_content=True, max_entries=500)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)


    def test_add_complex(self):
        q = self.intact.session.query(intact.ProteinComplex).filter_by(name='CAX1-CAX3 complex').first()

        self.assertEqual(q.ncbi, '3702')
        self.assertEqual(q.evidence, 'intact:EBI-2292448')

    def test_add_interactions(self):

        q = self.intact.session.query(intact.ProteinInteraction).get(99)

        self.assertEqual(q.protein_a, 'Q9Y6M1')
        self.assertEqual(q.protein_b, 'Q14103-2')
        self.assertEqual(q.gene_a, 'IGF2BP2')
        self.assertEqual(q.gene_b, 'HNRNPD')
        self.assertEqual(q.type_a, q.type_b)
        self.assertEqual(q.role_a, q.role_b)
        self.assertEqual(q.method, 'two hybrid')
        self.assertEqual(q.publication, '12674497')
        self.assertEqual(q.publication_author, 'Moraes et al. (2003)')
