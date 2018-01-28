""" Tests of bio_portal

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-05
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.data_source import bio_portal
import os
import shutil
import sys
import tempfile
import unittest


class TestBioPortal(unittest.TestCase):

    def setUp(self):
        self.cache_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.cache_dir)

    def test_get_ontologies_filename(self):
        bp = bio_portal.BioPortal(cache_dir=self.cache_dir)
        self.assertEqual(bp.get_ontologies_filename(), os.path.join(self.cache_dir, 'bio_portal.ontologies.json'))

    def test_download_ontologies(self):
        bp = bio_portal.BioPortal(cache_dir=self.cache_dir)
        bp.download_ontologies()
        self.assertTrue(os.path.isfile(os.path.join(self.cache_dir, 'bio_portal.ontologies.json')))

    def test_get_ontologies(self):
        bp = bio_portal.BioPortal(cache_dir=self.cache_dir)
        ontologies = bp.get_ontologies()
        self.assertTrue(os.path.isfile(os.path.join(self.cache_dir, 'bio_portal.ontologies.json')))
        self.assertIn('SBO', ontologies)

    def test_get_ontology_filename(self):
        bp = bio_portal.BioPortal(cache_dir=self.cache_dir)
        self.assertEqual(bp.get_ontology_filename('SBO'), os.path.join(self.cache_dir, 'SBO'))

    def test_download_ontology(self):
        bp = bio_portal.BioPortal(cache_dir=self.cache_dir)
        bp.download_ontology('SBO')
        self.assertTrue(os.path.isfile(os.path.join(self.cache_dir, 'SBO.obo')))
        self.assertTrue(os.path.isfile(os.path.join(self.cache_dir, 'SBO.py{}.pkl'.format(sys.version_info[0]))))

    def test_get_ontology(self):
        bp = bio_portal.BioPortal(cache_dir=self.cache_dir)
        onto = bp.get_ontology('SBO')
        self.assertTrue(os.path.isfile(os.path.join(self.cache_dir, 'SBO.obo')))
        self.assertTrue(os.path.isfile(os.path.join(self.cache_dir, 'SBO.py{}.pkl'.format(sys.version_info[0]))))
        self.assertEqual(onto['SBO:0000001'].name, 'rate law')
