""" Tests of bio_portal

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-05
:Copyright: 2017, Karr Lab
:License: MIT
"""

from datanator.data_source import bio_portal
import os
import shutil
import sys
import tempfile
import unittest


class TestBioPortal(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        bio_portal.BioPortal(
            ontologies=['BTO.obo', 'SBO.obo'],
            cache_dirname=cls.cache_dirname, download_backups=False, clear_content=True, load_content=True)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    def test_get_ontologies_filename(self):
        bp = bio_portal.BioPortal(ontologies=['BTO.obo', 'SBO.obo'], cache_dirname=self.cache_dirname)
        self.assertEqual(bp.get_ontologies_filename(), os.path.join(self.cache_dirname, 'bio_portal.ontologies.json'))

    def test_download_ontologies(self):
        bp = bio_portal.BioPortal(ontologies=['BTO.obo', 'SBO.obo'], cache_dirname=self.cache_dirname)
        bp.download_ontologies()
        self.assertTrue(os.path.isfile(os.path.join(self.cache_dirname, 'bio_portal.ontologies.json')))

    def test_get_ontologies(self):
        bp = bio_portal.BioPortal(ontologies=['BTO.obo', 'SBO.obo'], cache_dirname=self.cache_dirname)
        ontologies = bp.get_ontologies()
        self.assertTrue(os.path.isfile(os.path.join(self.cache_dirname, 'bio_portal.ontologies.json')))
        self.assertIn('SBO', ontologies)

    def test_get_ontology_filename(self):
        bp = bio_portal.BioPortal(ontologies=['BTO.obo', 'SBO.obo'], cache_dirname=self.cache_dirname)
        self.assertEqual(bp.get_ontology_filename('SBO'), os.path.join(self.cache_dirname, 'SBO'))

    def test_download_ontology(self):
        bp = bio_portal.BioPortal(ontologies=['BTO.obo', 'SBO.obo'], cache_dirname=self.cache_dirname)
        bp.download_ontology('SBO')
        self.assertTrue(os.path.isfile(os.path.join(self.cache_dirname, 'SBO.obo')))
        self.assertTrue(os.path.isfile(os.path.join(self.cache_dirname, 'SBO.py{}.pkl'.format(sys.version_info[0]))))

    def test_get_ontology(self):
        bp = bio_portal.BioPortal(ontologies=['BTO.obo', 'SBO.obo'], cache_dirname=self.cache_dirname)
        onto = bp.get_ontology('SBO')
        self.assertTrue(os.path.isfile(os.path.join(self.cache_dirname, 'SBO.obo')))
        self.assertTrue(os.path.isfile(os.path.join(self.cache_dirname, 'SBO.py{}.pkl'.format(sys.version_info[0]))))
        self.assertEqual(onto['SBO:0000001'].name, 'rate law')


class TestBioPortalFromBackup(unittest.TestCase):

    @classmethod
    def setUp(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        bio_portal.BioPortal(cache_dirname=cls.cache_dirname, download_backups=True)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    def test_download_backup(self):
        bp = bio_portal.BioPortal(ontologies=['DOID.obo'], cache_dirname=self.cache_dirname)
        onto = bp.get_ontology('DOID')
