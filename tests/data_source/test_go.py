""" Tests of GO
:Author:  Saahith Pochiraju  <saahith116@gmail.com>
:Date: 2018-04-09
:Copyright: 2018, Karr Lab
:License: MIT
"""

from kinetic_datanator.data_source import go
import datetime
import dateutil
import os
import shutil
import tempfile
import unittest


class TestChebiDownload(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.go_ontology = go.GO(download_backups=False, load_content=True, clear_content=False)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)

    def test_get_names(self):
        self.assertEqual(self.go_ontology.get_name(self.go_ontology.graph, "GO:2001141"), 'regulation of RNA biosynthetic process')
        self.assertEqual(self.go_ontology.get_name(self.go_ontology.graph, 'GO:2000479'), 'regulation of cAMP-dependent protein kinase activity')
