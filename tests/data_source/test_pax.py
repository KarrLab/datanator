# -*- coding: utf-8 -*-

""" Test of pax database

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-08-1
:Copyright: 2017, Karr Lab
:License: MIT
"""
import unittest
from kinetic_datanator.data_source import pax
from sqlalchemy.orm import sessionmaker
import tempfile
import shutil


class TestPaxDBDownload(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.src = pax.Pax(cache_dirname=self.cache_dirname,
                           load_content=False, clear_content=False,
                           download_backups=True, verbose=True,
                           max_entries=5)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)

    def test_downloaded_content(self):
        session = self.src.session

        obs = session.query(pax.Observation).get(2)
        self.assertIsInstance(obs.dataset_id, int)

        data = obs.dataset

        self.assertEqual(data.taxon_ncbi_id, 882)
        self.assertIsInstance(data.score, float)
        self.assertIsInstance(data.weight, int)
        self.assertIsInstance(data.coverage, int)

        prot = obs.protein

        self.assertEqual(prot.protein_id, obs.protein_id)

        refined_data = session.query(pax.Dataset).filter(
            pax.Dataset.file_name == '882/882-Desulfo_Form_Exp_SC_zhang_2006.txt').first()
        self.assertEqual(refined_data.score, 2.47)
        self.assertEqual(refined_data.weight, 100)
        self.assertEqual(refined_data.taxon_ncbi_id, 882)


class TestPaxDBCreation(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.src = pax.Pax(cache_dirname=self.cache_dirname,
                           load_content=True, clear_content=True,
                           download_backups=False, verbose=True,
                           max_entries=2)

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_load_some_content(self):
        session = self.src.session

        obs = session.query(pax.Observation).get(2)
        self.assertIsInstance(obs.dataset_id, int)

        data = obs.dataset

        self.assertEqual(data.taxon_ncbi_id, 882)
        self.assertIsInstance(data.score, float)
        self.assertIsInstance(data.weight, int)
        self.assertIsInstance(data.coverage, int)

        prot = obs.protein

        self.assertEqual(prot.protein_id, obs.protein_id)

        refined_data = session.query(pax.Dataset).filter(
            pax.Dataset.file_name == '882/882-Desulfo_Form_Stat_SC_zhang_2006.txt').first()
        self.assertEqual(refined_data.score, 0.61)
        self.assertEqual(refined_data.weight, 20)
        self.assertEqual(refined_data.taxon_ncbi_id, 882)
