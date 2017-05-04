""" Tests of data structures

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-14
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_structs
import unittest


class TestCrossReferenceAssignmentMethod(unittest.TestCase):

    def test(self):
        self.assertEqual(data_structs.CrossReferenceAssignmentMethod['manual'].value, 0)
        self.assertEqual(data_structs.CrossReferenceAssignmentMethod['predicted'].value, 1)


class TestCrossReference(unittest.TestCase):

    def test_init(self):
        xr = data_structs.CrossReference(source='src', id='identifier', relevance=2.,
                                         assignment_method=data_structs.CrossReferenceAssignmentMethod.manual)
        self.assertEqual(xr.source, 'src')
        self.assertEqual(xr.id, 'identifier')
        self.assertEqual(xr.relevance, 2.)
        self.assertEqual(xr.assignment_method, data_structs.CrossReferenceAssignmentMethod.manual)
