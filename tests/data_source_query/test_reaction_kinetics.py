# -*- coding: utf-8 -*-

""" Tests of sabio_rk

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-26
:Copyright: 2017, Karr Lab
:License: MIT
"""

import unittest


class TestSabioRkUtil(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.sabio_rk = sabio_rk.SabioRkUtil()

    def test_get_compounds_by_structure(self):
        sabio_rk = self.sabio_rk

        self.assertRaises(ValueError, sabio_rk.get_compounds_by_structure, '')

        structure = 'InChI=1S/C6H13O9P/c7-3-2(1-14-16(11,12)13)15-6(10)5(9)4(3)8/h2-10H,1H2,(H2,11,12,13)/t2-,3-,4+,5-,6+/m1/s1'
        compounds = sabio_rk.get_compounds_by_structure(structure,
                                                        check_protonation=True,
                                                        check_double_bonds=True,
                                                        check_stereochemistry=True)
        self.assertEqual(len(compounds), 1)
        self.assertEqual(compounds[0].id, 24)

        structure = 'InChI=1S/C6H13O9P/c7-3-2(1-14-16(11,12)13)15-6(10)5(9)4(3)8/h2-10H,1H2,(H2,11,12,13)/t2-,3-,4+,5-,6+/m1/s1'
        compounds = sabio_rk.get_compounds_by_structure(structure)
        self.assertEqual(len(compounds), 10)
        self.assertEqual(set(c.id for c in compounds), set([24, 1323, 1357, 1366, 1404, 1405, 1493, 1522, 1524, 24709]))
