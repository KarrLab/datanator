"""
Test for IntAct database module

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-08-13
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""

from datanator.data_source import intact
import shutil
import unittest
import tempfile


class TestFromServerIntAct(unittest.TestCase):
    """
    Testing IntAct database from cached server
    """

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.intact = intact.IntAct(cache_dirname=cls.cache_dirname,
                                    clear_content=False,
                                    download_backups=True,
                                    load_content=False,
                                    verbose=True)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    def test_complexes(self):
        q = self.intact.session.query(intact.ProteinComplex).get('CPX-863')
        self.assertEqual(q.name, 'INO80 chromatin remodeling complex')
        self.assertEqual(q.ncbi, '559292')
        self.assertEqual(q.evidence, 'intact:EBI-515508')

    def test_interactions(self):
        q = self.intact.session.query(intact.ProteinInteraction).filter_by(protein_a='P21127', protein_b='Q13153').first()
        self.assertEqual(q.gene_a, 'CDK11B')
        self.assertEqual(q.gene_b, 'PAK1')
        self.assertEqual(q.type_a, 'protein')
        self.assertEqual(q.type_b, 'protein')
        self.assertEqual(q.role_a, 'unspecified role')
        self.assertEqual(q.role_b, 'unspecified role')
        self.assertEqual(q.method, 'anti tag coimmunoprecipitation')
        self.assertEqual(q.publication, '12624090')
        self.assertEqual(q.publication_author, 'Chen et al. (2003)')

        q = self.intact.session.query(intact.ProteinInteraction).filter_by(protein_a='Q8IVG9').count()
        self.assertEqual(q, 6)


class TestLoadingIntAct(unittest.TestCase):
    """
    Testing loading IntAct database
    """

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.intact = intact.IntAct(cache_dirname=cls.cache_dirname, download_backups=False, load_content=True, max_entries=500)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    def test_complexes(self):
        q = self.intact.session.query(intact.ProteinComplex).filter_by(identifier='CPX-1326').first()

        self.assertEqual(q.name, 'CAX1-CAX3 complex')
        self.assertEqual(q.ncbi, '3702')
        self.assertEqual(q.evidence, 'intact:EBI-2292448')

    def test_interactions(self):
        q = self.intact.session.query(intact.ProteinInteraction).filter_by(index=99).first()
        self.assertEqual(q.protein_a, 'Q14103-2')
        self.assertEqual(q.protein_b, 'Q9Y6M1')
        self.assertEqual(q.gene_b, 'IGF2BP2')
        self.assertEqual(q.gene_a, 'HNRNPD')
        self.assertEqual(q.type_a, q.type_b)
        self.assertEqual(q.role_a, q.role_b)
        self.assertEqual(q.method, 'two hybrid')
        self.assertEqual(q.publication, '12674497')
        self.assertEqual(q.publication_author, 'Moraes et al. (2003)')
