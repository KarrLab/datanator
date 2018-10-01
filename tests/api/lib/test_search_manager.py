""" Test of text search manager

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2018-08-08
:Copyright: 2018, Karr Lab
:License: MIT
"""

from datanator.api.lib.search.manager import search_manager
import unittest

class TestSearchManager(unittest.TestCase):

    def test_search(self):
        dict= search_manager.search('2-Oxopentanoate')
        self.assertGreater(len(dict['Metabolite']), 0)
        self.assertGreater(len(dict['Reaction']), 0)
        self.assertEqual(len(dict['ProteinSubunit']), 0)
        self.assertEqual(len(dict['ProteinComplex']), 0)
