# -*- coding: utf-8 -*-

""" Test of UniProt module

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-08-15
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""
from kinetic_datanator.data_source import uniprot
import shutil
import tempfile
import unittest


class TestServerDownloadUniprot(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.uni = uniprot.Uniprot(cache_dirname=self.cache_dirname, download_backups=False, load_content=True,
                                   max_entries=10)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)

    def test_proper_loading(self):
        count = self.uni.session.query(uniprot.UniprotData).count()
        self.assertEqual(count, 10)

    def test_load(self):
        prot = self.uni.session.query(uniprot.UniprotData).filter_by(uniprot_id='Q12181').all()
        self.assertEqual(len(prot), 1)
        self.assertEqual(prot[0].entry_name, 'NDOR1_YEAST')
        self.assertEqual(prot[0].protein_name,
                         'NADPH-dependent diflavin oxidoreductase 1 (EC 1.18.1.-) (NADPH-dependent FMN and FAD-containing oxidoreductase)')
        self.assertEqual(prot[0].entrez_id, 856161)

    def test_load1(self):
        prot = self.uni.session.query(uniprot.UniprotData).filter_by(uniprot_id='Q8GHV6').all()
        self.assertEqual(len(prot), 1)
        self.assertEqual(prot[0].mass, 37570)
        self.assertEqual(prot[0].ec_number, None)

    def test_load2(self):
        prot = self.uni.session.query(uniprot.UniprotData).filter_by(uniprot_id='B5BAY2').all()
        self.assertEqual(len(prot), 1)
        self.assertEqual(prot[0].entrez_id, None)
        self.assertEqual(prot[0].ec_number, '2.7.4.6')


class TestUniprotFromBackup(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.uni = uniprot.Uniprot(cache_dirname=self.cache_dirname, download_backups=True)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)

    def test_proper_loading(self):
        count = self.uni.session.query(uniprot.UniprotData).count()
        self.assertGreater(count, 500000)
