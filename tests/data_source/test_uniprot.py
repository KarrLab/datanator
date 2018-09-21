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
