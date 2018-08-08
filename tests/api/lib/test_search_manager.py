""" Test of text search manager

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2018-08-08
:Copyright: 2018, Karr Lab
:License: MIT
"""

from kinetic_datanator.api.lib.search.manager import search_manager
import unittest

class TestSearchManagern(unittest.TestCase):

    def test_search(self):
        dict, list = search_manager.search('2-Oxopentanoate')
        self.assertGreater(len(dict['Compound']), 0)
        self.assertGreater(len(list), 0)
