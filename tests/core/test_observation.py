# -*- coding: utf-8 -*-

""" Test of observation filtering and ordering

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-10
:Copyright: 2017, Karr Lab
:License: MIT
"""

from datetime import datetime
from kinetic_datanator.core import observation
import unittest


class TestConsensusMethod(unittest.TestCase):

    def test(self):
        self.assertEqual(observation.ConsensusMethod['manual'].value, 0)
        self.assertEqual(observation.ConsensusMethod['mean'].value, 1)


class TestConsensus(unittest.TestCase):

    def test(self):
        c = observation.Consensus(
            observable=observation.Property(id='volume'),
            value=1.1,
            error=0.3,
            units='uM',
            method=observation.ConsensusMethod.mean,
            evidence=[
                observation.Evidence(value=observation.ObservedValue(), relevance=10.),
                observation.Evidence(value=observation.ObservedValue(), relevance=5.),
            ],
            user='jonrkarr',
            date=datetime.utcnow(),
        )


class TestEvidence(unittest.TestCase):

    def test(self):
        value = observation.ObservedValue()
        evidence = observation.Evidence(
            value=value,
            relevance=10.,
        )
        self.assertEqual(evidence.value, value)
        self.assertEqual(evidence.relevance, 10.)


class TestObservation(unittest.TestCase):

    def test_Observation(self):
        o = observation.Observation()
        o.genetics = observation.Genetics(taxon='Mycoplasma pneumoniae', variation='ΔMPN001')
        o.environment = observation.Environment(temperature=37, ph=7., media='Hayflick')
        o.reference = observation.Reference(title='title', author='author', year=2017, volume=1, number=1, pages='1-10')

        observable = observation.Property(
            id='K_m',
            parent=observation.Specie(
                id='ATP',
                parent=observation.Compartment(
                    id='c',
                    parent=observation.Reaction(
                        id='AtpSynthase'
                    )
                )
            )
        )

        ov = observation.ObservedValue(
            observable=observable,
            value=1.0,
            error=0.5,
            units='U/mg',
            method=observation.ExperimentalMethod(name='assay', description='description of assay'),
        )
        o.values.append(ov)

        o.validate()

        self.assertEqual(o.values, [ov])
        self.assertEqual(ov.observation, o)
        self.assertEqual(ov.observable.id, 'K_m')
        self.assertEqual(ov.observable.parent.id, 'ATP')
        self.assertEqual(ov.observable.parent.parent.id, 'c')
        self.assertEqual(ov.observable.parent.parent.parent.id, 'AtpSynthase')
        self.assertEqual(ov.value, 1.0)
        self.assertEqual(ov.error, 0.5)
        self.assertEqual(ov.units, 'U/mg')
        self.assertEqual(ov.method.name, 'assay')
        self.assertEqual(ov.method.description, 'description of assay')

        self.assertEqual(o.genetics.taxon, 'Mycoplasma pneumoniae')
        self.assertEqual(o.genetics.variation, 'ΔMPN001')

        self.assertEqual(o.environment.temperature, 37.)
        self.assertEqual(o.environment.ph, 7.)
        self.assertEqual(o.environment.media, 'Hayflick')

        self.assertEqual(o.reference.title, 'title')
        self.assertEqual(o.reference.author, 'author')
        self.assertEqual(o.reference.year, 2017)
        self.assertEqual(o.reference.volume, 1)
        self.assertEqual(o.reference.number, 1)
        self.assertEqual(o.reference.pages, '1-10')


class TestReaction(unittest.TestCase):
    adp = 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'
    atp = 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'
    h = '[H+]'
    h2o = 'O'
    pi = 'OP([O-])([O-])=O'

    def make_reaction(self):
        return observation.Reaction(participants=[
            observation.ReactionParticipant(
                specie=observation.Specie(structure=self.atp, id='atp'),
                compartment=observation.Compartment(id='c'),
                coefficient=-1),
            observation.ReactionParticipant(
                specie=observation.Specie(structure=self.h2o, id='h2o'),
                compartment=observation.Compartment(id='c'),
                coefficient=-1),
            observation.ReactionParticipant(
                specie=observation.Specie(structure=self.adp, id='adp'),
                compartment=observation.Compartment(id='c'),
                coefficient=1),
            observation.ReactionParticipant(
                specie=observation.Specie(structure=self.pi, id='pi'),
                compartment=observation.Compartment(id='c'),
                coefficient=1),
            observation.ReactionParticipant(
                specie=observation.Specie(structure=self.h, id='h'),
                compartment=observation.Compartment(id='c'),
                coefficient=1),
        ])

    def test_init(self):
        rxn = self.make_reaction()
        self.assertEqual(len(rxn.participants), 5)

    def test_id(self):
        rxn = observation.Reaction(id='rxn')
        self.assertEqual(rxn.id, 'rxn')

    def test_name(self):
        rxn = observation.Reaction(name='rxn')
        self.assertEqual(rxn.name, 'rxn')

    def test_normalize(self):
        rxn = observation.Reaction(participants=[
            observation.ReactionParticipant(
                specie=observation.Specie(structure=self.atp, id='atp'),
                compartment=observation.Compartment(id='c'),
                coefficient=-1),
            observation.ReactionParticipant(
                specie=observation.Specie(structure=self.atp, id='atp'),
                compartment=observation.Compartment(id='c'),
                coefficient=-1),
            observation.ReactionParticipant(
                specie=observation.Specie(structure=self.atp, id='atp'),
                compartment=observation.Compartment(id='c'),
                coefficient=1),
            observation.ReactionParticipant(
                specie=observation.Specie(structure=self.atp, id='atp'),
                compartment=observation.Compartment(id='c'),
                coefficient=0),
            observation.ReactionParticipant(
                specie=observation.Specie(structure=self.atp, id='atp'),
                compartment=observation.Compartment(id='e'),
                coefficient=-1),
            observation.ReactionParticipant(
                specie=observation.Specie(structure=self.adp, id='adp'),
                compartment=observation.Compartment(id='c'),
                coefficient=1),
            observation.ReactionParticipant(
                specie=observation.Specie(structure=self.adp, id='adp'),
                compartment=observation.Compartment(id='c'),
                coefficient=3),
        ])

        ordered_participants = rxn.get_ordered_participants()

        self.assertEqual(len(ordered_participants), 4)

        part = list(filter(lambda p: p.specie.id == 'atp' and p.compartment.id == 'c' and p.coefficient < 0, ordered_participants))
        self.assertEqual(len(part), 1)
        self.assertEqual(part[0].coefficient, -2)

        part = list(filter(lambda p: p.specie.id == 'atp' and p.compartment.id == 'c' and p.coefficient > 0, ordered_participants))
        self.assertEqual(len(part), 1)
        self.assertEqual(part[0].coefficient, 1)

        part = list(filter(lambda p: p.specie.id == 'atp' and p.compartment.id == 'e', ordered_participants))
        self.assertEqual(len(part), 1)
        self.assertEqual(part[0].coefficient, -1)

        part = list(filter(lambda p: p.specie.id == 'adp', ordered_participants))
        self.assertEqual(len(part), 1)
        self.assertEqual(part[0].coefficient, 4)

    def test_get_reactants(self):
        rxn = self.make_reaction()
        self.assertEqual(rxn.get_reactants(), rxn.participants[0:2])

    def test_get_products(self):
        rxn = self.make_reaction()
        self.assertEqual(rxn.get_products(), rxn.participants[2:])

    def test_get_reactant_product_pairs(self):
        rxn = self.make_reaction()
        pairs = rxn.get_reactant_product_pairs()

        self.assertEqual(pairs[0][0].specie.id, 'atp')
        self.assertEqual(pairs[0][1].specie.id, 'adp')

        self.assertEqual(pairs[1][0].specie.id, 'h2o')
        self.assertIn(pairs[1][1].specie.id, ['h', 'pi'])

        self.assertEqual(pairs[2][0], None)
        self.assertIn(pairs[2][1].specie.id, ['h', 'pi'])

    def test_get_ec_number(self):
        rxn = observation.Reaction(cross_references=[
            observation.Resource(namespace='xx', id='yy', relevance=20.,
                                 assignment_method=observation.ResourceAssignmentMethod.predicted),
            observation.Resource(namespace='ec-code', id='1.1.1.1', relevance=20.,
                                 assignment_method=observation.ResourceAssignmentMethod.predicted),
            observation.Resource(namespace='ec-code', id='1.1.1.2', relevance=30.,
                                 assignment_method=observation.ResourceAssignmentMethod.predicted),
            observation.Resource(namespace='ec-code', id='1.1.1.3', relevance=10.,
                                 assignment_method=observation.ResourceAssignmentMethod.predicted),
        ])
        self.assertEqual(rxn.get_ec_number(), '1.1.1.2')
        self.assertEqual(rxn.get_ec_numbers(), rxn.cross_references[1:])
        self.assertEqual(rxn.get_manual_ec_numbers(), [])
        self.assertEqual(rxn.get_predicted_ec_numbers(), rxn.cross_references[1:])

        rxn = observation.Reaction(cross_references=[
            observation.Resource(namespace='ec-code', id='1.1.1.1', relevance=20.,
                                 assignment_method=observation.ResourceAssignmentMethod.manual),
            observation.Resource(namespace='ec-code', id='1.1.1.2', relevance=30.,
                                 assignment_method=observation.ResourceAssignmentMethod.predicted),
            observation.Resource(namespace='ec-code', id='1.1.1.3', relevance=10.,
                                 assignment_method=observation.ResourceAssignmentMethod.predicted),
            observation.Resource(namespace='xx', id='yy', relevance=20.,
                                 assignment_method=observation.ResourceAssignmentMethod.predicted),
        ])
        self.assertEqual(rxn.get_ec_number(), '1.1.1.1')
        self.assertEqual(rxn.get_ec_numbers(), rxn.cross_references[0:-1])
        self.assertEqual(rxn.get_manual_ec_numbers(), rxn.cross_references[0:1])
        self.assertEqual(rxn.get_predicted_ec_numbers(), rxn.cross_references[1:3])

        rxn = observation.Reaction(cross_references=[
            observation.Resource(namespace='ec-code', id='1.1.1.1', relevance=20.,
                                 assignment_method=observation.ResourceAssignmentMethod.manual),
            observation.Resource(namespace='ec-code', id='1.1.1.2', relevance=30.,
                                 assignment_method=observation.ResourceAssignmentMethod.predicted),
            observation.Resource(namespace='ec-code', id='1.1.1.3', relevance=10.,
                                 assignment_method=observation.ResourceAssignmentMethod.manual),
        ])
        self.assertRaises(ValueError, rxn.get_ec_number)

        rxn = observation.Reaction(cross_references=[
            observation.Resource(namespace='ec2', id='1.1.1.1', relevance=20.,
                                 assignment_method=observation.ResourceAssignmentMethod.manual),
            observation.Resource(namespace='ec-code', id='1.1.1.2', relevance=30.,
                                 assignment_method=observation.ResourceAssignmentMethod.predicted),
            observation.Resource(namespace='ec2', id='1.1.1.3', relevance=10.,
                                 assignment_method=observation.ResourceAssignmentMethod.manual),
        ])
        self.assertEqual(rxn.get_ec_number(), '1.1.1.2')

        rxn = observation.Reaction(cross_references=[])
        self.assertEqual(rxn.get_ec_number(), '')


class TestGenetics(unittest.TestCase):

    def test(self):
        genetics = observation.Genetics(taxon='Mycoplasma pneumoniae', variation='ΔMPN001')
        self.assertFalse(genetics.is_wildtype())
        self.assertTrue(genetics.is_variant())

        genetics = observation.Genetics(taxon='Mycoplasma pneumoniae', variation='')
        self.assertTrue(genetics.is_wildtype())
        self.assertFalse(genetics.is_variant())


class TestResourceAssignmentMethod(unittest.TestCase):

    def test(self):
        self.assertEqual(observation.ResourceAssignmentMethod['manual'].value, 0)
        self.assertEqual(observation.ResourceAssignmentMethod['predicted'].value, 1)


class TestResource(unittest.TestCase):

    def test_init(self):
        xr = observation.Resource(namespace='src', id='identifier', relevance=2.,
                                  assignment_method=observation.ResourceAssignmentMethod.manual)
        self.assertEqual(xr.namespace, 'src')
        self.assertEqual(xr.id, 'identifier')
        self.assertEqual(xr.relevance, 2.)
        self.assertEqual(xr.assignment_method, observation.ResourceAssignmentMethod.manual)
