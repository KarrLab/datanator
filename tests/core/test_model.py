""" Test of observation filtering and ordering

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-10
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import model
import unittest


class TestModel(unittest.TestCase):

    def test_Parameter(self):
        p = model.Parameter(component='AtpSynthase', attribute='v_max', value=1.0, units='U/mg', consensus_method='mean')

        p.validate()

        self.assertEqual(p.component, 'AtpSynthase')
        self.assertEqual(p.attribute, 'v_max')
        self.assertEqual(p.value, 1.0)
        self.assertEqual(p.units, 'U/mg')
        self.assertEqual(p.evidence, [])
        self.assertEqual(p.consensus_method, 'mean')
