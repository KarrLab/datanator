""" Tests of the compartment utilities

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.util import compartment_util
import unittest


class TestCompartmentUtil(unittest.TestCase):

    def test_Compartment(self):
        c = compartment_util.Compartment('c')
        self.assertEqual(c.name, 'c')
