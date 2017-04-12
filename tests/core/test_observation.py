# -*- coding: utf-8 -*-

""" Test of observation filtering and ordering

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-10
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import model
from kinetic_datanator.core import observation
import unittest


class TestObservation(unittest.TestCase):

    def test_Observation(self):
        o = observation.Observation(component='AtpSynthase', attribute='v_max', value=1.0, units='U/mg')
        o.taxon = observation.Taxon(name='Mycoplasma pneumoniae', perturbations='ΔMPN001')
        o.environment = observation.Environment(temperature=37, ph=7., media='Hayflick')
        o.method = observation.ExperimentalMethod(name='assay', description='description of assay')
        o.reference = observation.Reference(title='title', author='author', year=2017, volume=1, number=1, pages='1-10')
        o.parameters.create(component='AtpSynthase', attribute='k_cat', value=2.0, units='1/s', consensus_method='mean')

        o.validate()

        self.assertEqual(o.component, 'AtpSynthase')
        self.assertEqual(o.attribute, 'v_max')
        self.assertEqual(o.value, 1.0)
        self.assertEqual(o.units, 'U/mg')

        self.assertEqual(o.taxon.name, 'Mycoplasma pneumoniae')
        self.assertEqual(o.taxon.perturbations, 'ΔMPN001')

        self.assertEqual(o.environment.temperature, 37.)
        self.assertEqual(o.environment.ph, 7.)
        self.assertEqual(o.environment.media, 'Hayflick')

        self.assertEqual(o.method.name, 'assay')
        self.assertEqual(o.method.description, 'description of assay')

        self.assertEqual(o.reference.title, 'title')
        self.assertEqual(o.reference.author, 'author')
        self.assertEqual(o.reference.year, 2017)
        self.assertEqual(o.reference.volume, 1)
        self.assertEqual(o.reference.number, 1)
        self.assertEqual(o.reference.pages, '1-10')

        self.assertEqual(len(o.parameters), 1)
        p = o.parameters[0]
        self.assertEqual(p.component, 'AtpSynthase')
        self.assertEqual(p.attribute, 'k_cat')
        self.assertEqual(p.value, 2.0)
        self.assertEqual(p.units, '1/s')
        self.assertEqual(p.evidence, [o])
        self.assertEqual(p.consensus_method, 'mean')

    def test_Taxon(self):
        taxon = observation.Taxon(name='Mycoplasma pneumoniae', perturbations='ΔMPN001')
        self.assertFalse(taxon.is_wildtype())
        self.assertTrue(taxon.is_mutant())

        taxon = observation.Taxon(name='Mycoplasma pneumoniae', perturbations='')
        self.assertTrue(taxon.is_wildtype())
        self.assertFalse(taxon.is_mutant())
