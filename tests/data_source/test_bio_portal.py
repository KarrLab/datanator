""" Tests of bio_portal

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-05
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.data_source import bio_portal
import os
import shutil
import tempfile
import unittest


class TestBioPortal(unittest.TestCase):

    def setUp(self):
        self.cache_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.cache_dir)

    def test_get_ontology_filename(self):
        bp = bio_portal.BioPortal(cache_dir=self.cache_dir)
        self.assertEqual(bp.get_ontology_filename('SBO'), os.path.join(self.cache_dir, 'SBO.obo'))

    def test_download_ontology(self):
        bp = bio_portal.BioPortal(cache_dir=self.cache_dir)
        onto = bp.download_ontology('SBO')
        self.assertTrue(os.path.isfile(os.path.join(self.cache_dir, 'SBO.obo')))

    def test_load_ontology(self):
        bp = bio_portal.BioPortal(cache_dir=self.cache_dir)
        onto = bp.load_ontology('SBO')
        self.assertTrue(os.path.isfile(os.path.join(self.cache_dir, 'SBO.obo')))
        self.assertEqual(onto['SBO:0000001'].name, 'rate law')
