# -*- coding: utf-8 -*-

""" Test of common schema

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-07-31
:Copyright: 2017, Karr Lab
:License: MIT
"""
import unittest
from kinetic_datanator.data_source import uniprot
import tempfile
import shutil


class TestPortionUniprot(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.uni = uniprot.Uniprot(cache_dirname=self.cache_dirname)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)

    def test_load(self):
        test = self.uni.session.query(uniprot.UniprotData).filter_by(uniprot_id = 'B5F7M5').all()
        self.assertEqual(len(test),1)
        self.assertEqual(test[0].entry_name, 'AAEA_SALA4')
        self.assertEqual(test[0].protein_name,
        'p-hydroxybenzoic acid efflux pump subunit AaeA (pHBA efflux pump protein A)')

    def test_load1(self):
        test = self.uni.session.query(uniprot.UniprotData).filter_by(uniprot_id = 'Q39RX3').all()
        self.assertEqual(len(test),1)
        self.assertEqual(test[0].mass, '41,582')
        self.assertEqual(test[0].ec_number, '2.3.1.31')

    def test_load2(self):
        test = self.uni.session.query(uniprot.UniprotData).filter_by(uniprot_id = 'A6USA0').all()
        self.assertEqual(len(test),1)
        self.assertEqual(test[0].entrez_id, 5324858)
        self.assertEqual(test[0].ec_number, '2.1.1.166')
