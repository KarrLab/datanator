""" Tests of Chebi
:Author:  Saahith Pochiraju  <saahith116@gmail.com>
:Date: 2018-04-09
:Copyright: 2018, Karr Lab
:License: MIT
"""

from kinetic_datanator.data_source import chebi
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
        self.cheb = chebi.Chebi( download_backups=False, load_content=True, clear_content=False)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)

    def test_chebi_names(self):
        self.assertEqual(self.cheb.get_name(self.cheb.graph, 'CHEBI:78665'), 'royal jelly')
        self.assertEqual(self.cheb.get_name(self.cheb.graph, 'CHEBI:17203'), 'L-proline')
        self.assertEqual(self.cheb.get_name(self.cheb.graph, 'CHEBI:22658'), 'aspartate family amino acid')
        self.assertEqual(self.cheb.get_name(self.cheb.graph, 'CHEBI:78668'), '(E)-10-hydroxydec-2-enoic acid')
