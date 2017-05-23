# -*- coding: utf-8 -*-

""" Test of observation filtering and ordering

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-10
:Copyright: 2017, Karr Lab
:License: MIT
"""

from datetime import datetime
from kinetic_datanator.core import data_model
import unittest


class TestConsensusMethod(unittest.TestCase):

    def test(self):
        self.assertEqual(data_model.ConsensusMethod['manual'].value, 0)
        self.assertEqual(data_model.ConsensusMethod['mean'].value, 1)


class TestConsensus(unittest.TestCase):

    def test(self):
        c = data_model.Consensus(
            observable=data_model.Observable(),
            value=1.1,
            error=0.3,
            units='uM',
            method=data_model.ConsensusMethod.mean,
            evidence=[
                data_model.Evidence(value=data_model.ObservedValue(), relevance=10.),
                data_model.Evidence(value=data_model.ObservedValue(), relevance=5.),
            ],
            user='jonrkarr',
            date=datetime.utcnow(),
        )


class TestEvidence(unittest.TestCase):

    def test(self):
        value = data_model.ObservedValue()
        evidence = data_model.Evidence(
            value=value,
            relevance=10.,
        )
        self.assertEqual(evidence.value, value)
        self.assertEqual(evidence.relevance, 10.)


class TestObservation(unittest.TestCase):

    def test_Observation(self):
        o = data_model.Observation()
        o.genetics = data_model.Genetics(taxon='Mycoplasma pneumoniae', variation='ΔMPN001')
        o.environment = data_model.Environment(temperature=37, ph=7., media='Hayflick')
        o.reference = data_model.Reference(title='title', author='author', year=2017, volume=1, number=1, pages='1-10')

        observable = data_model.Observable(
            interaction=data_model.Reaction(id='AtpSynthase'),
            specie=data_model.Specie(id='ATP'),
            compartment=data_model.Compartment(id='c'),
            property='K_m',
        )

        ov = data_model.ObservedValue(
            observable=observable,
            value=1.0,
            error=0.5,
            units='U/mg',
            method=data_model.ExperimentalMethod(name='assay', description='description of assay'),
        )
        o.values.append(ov)

        o.validate()

        self.assertEqual(o.values, [ov])
        self.assertEqual(ov.observation, o)
        self.assertEqual(ov.observable.interaction.id, 'AtpSynthase')
        self.assertEqual(ov.observable.specie.id, 'ATP')
        self.assertEqual(ov.observable.compartment.id, 'c')
        self.assertEqual(ov.observable.property, 'K_m')

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


class TestSpecie(unittest.TestCase):
    adp = 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'
    atp = 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'
    pi = 'InChI=1S/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)/p-2'

    def test_to_inchi(self):
        s = data_model.Specie(structure=self.pi)
        self.assertEqual(s.to_inchi(), self.pi)
        self.assertEqual(s.to_inchi(only_formula_and_connectivity=True), 'H3O4P/c1-5(2,3)4')

    def test_to_mol(self):
        s = data_model.Specie(structure=self.pi)
        self.assertIsInstance(s.to_mol(), str)
        self.assertNotEqual(s.to_mol(), '')

    def test_to_openbabel(self):
        s = data_model.Specie(structure=self.pi)
        self.assertEqual(s.to_openbabel().GetFormula(), 'HO4P--')

    def test_to_pybel(self):
        s = data_model.Specie(structure=self.pi)
        self.assertEqual(s.to_pybel().formula, 'HO4P--')

    def test_to_smiles(self):
        s = data_model.Specie(structure=self.pi)
        self.assertEqual(s.to_smiles(), '[O-]P(=O)(O)[O-]')

    def get_similarity(self):
        adp = data_model.Specie(structure=seld.adp)
        atp = data_model.Specie(structure=seld.atp)
        self.assertAlmostEqual(adp.get_similarity(atp), 0.955, places=3)


class TestPolymerSpecie(unittest.TestCase):

    def test(self):
        s = data_model.PolymerSpecie(sequence='ACGT')
        self.assertEqual(s.sequence, 'ACGT')

        s = data_model.DnaSpecie(sequence='ACGT')
        self.assertEqual(s.sequence, 'ACGT')

        s = data_model.ProteinSpecie(sequence='ACGT')
        self.assertEqual(s.sequence, 'ACGT')

        s = data_model.RnaSpecie(sequence='ACGT')
        self.assertEqual(s.sequence, 'ACGT')


class TestReaction(unittest.TestCase):
    adp = 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'
    atp = 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'
    h = '[H+]'
    h2o = 'O'
    pi = 'OP([O-])([O-])=O'

    def make_reaction(self):
        return data_model.Reaction(participants=[
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.atp, id='atp'),
                compartment=data_model.Compartment(id='c'),
                coefficient=-1),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.h2o, id='h2o'),
                compartment=data_model.Compartment(id='c'),
                coefficient=-1),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.adp, id='adp'),
                compartment=data_model.Compartment(id='c'),
                coefficient=1),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.pi, id='pi'),
                compartment=data_model.Compartment(id='c'),
                coefficient=1),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.h, id='h'),
                compartment=data_model.Compartment(id='c'),
                coefficient=1),
        ])

    def test_init(self):
        rxn = self.make_reaction()
        self.assertEqual(len(rxn.participants), 5)

    def test_id(self):
        rxn = data_model.Reaction(id='rxn')
        self.assertEqual(rxn.id, 'rxn')

    def test_name(self):
        rxn = data_model.Reaction(name='rxn')
        self.assertEqual(rxn.name, 'rxn')

    def test_normalize(self):
        rxn = data_model.Reaction(participants=[
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.atp, id='atp'),
                compartment=data_model.Compartment(id='c'),
                coefficient=-1),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.atp, id='atp'),
                compartment=data_model.Compartment(id='c'),
                coefficient=-1),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.atp, id='atp'),
                compartment=data_model.Compartment(id='c'),
                coefficient=1),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.atp, id='atp'),
                compartment=data_model.Compartment(id='c'),
                coefficient=0),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.atp, id='atp'),
                compartment=data_model.Compartment(id='e'),
                coefficient=-1),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.adp, id='adp'),
                compartment=data_model.Compartment(id='c'),
                coefficient=1),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.adp, id='adp'),
                compartment=data_model.Compartment(id='c'),
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
        rxn = data_model.Reaction(cross_references=[
            data_model.Resource(namespace='xx', id='yy', relevance=20.,
                                assignment_method=data_model.ResourceAssignmentMethod.predicted),
            data_model.Resource(namespace='ec-code', id='1.1.1.1', relevance=20.,
                                assignment_method=data_model.ResourceAssignmentMethod.predicted),
            data_model.Resource(namespace='ec-code', id='1.1.1.2', relevance=30.,
                                assignment_method=data_model.ResourceAssignmentMethod.predicted),
            data_model.Resource(namespace='ec-code', id='1.1.1.3', relevance=10.,
                                assignment_method=data_model.ResourceAssignmentMethod.predicted),
        ])
        self.assertEqual(rxn.get_ec_number(), '1.1.1.2')
        self.assertEqual(rxn.get_ec_numbers(), rxn.cross_references[1:])
        self.assertEqual(rxn.get_manual_ec_numbers(), [])
        self.assertEqual(rxn.get_predicted_ec_numbers(), rxn.cross_references[1:])

        rxn = data_model.Reaction(cross_references=[
            data_model.Resource(namespace='ec-code', id='1.1.1.1', relevance=20.,
                                assignment_method=data_model.ResourceAssignmentMethod.manual),
            data_model.Resource(namespace='ec-code', id='1.1.1.2', relevance=30.,
                                assignment_method=data_model.ResourceAssignmentMethod.predicted),
            data_model.Resource(namespace='ec-code', id='1.1.1.3', relevance=10.,
                                assignment_method=data_model.ResourceAssignmentMethod.predicted),
            data_model.Resource(namespace='xx', id='yy', relevance=20.,
                                assignment_method=data_model.ResourceAssignmentMethod.predicted),
        ])
        self.assertEqual(rxn.get_ec_number(), '1.1.1.1')
        self.assertEqual(rxn.get_ec_numbers(), rxn.cross_references[0:-1])
        self.assertEqual(rxn.get_manual_ec_numbers(), rxn.cross_references[0:1])
        self.assertEqual(rxn.get_predicted_ec_numbers(), rxn.cross_references[1:3])

        rxn = data_model.Reaction(cross_references=[
            data_model.Resource(namespace='ec-code', id='1.1.1.1', relevance=20.,
                                assignment_method=data_model.ResourceAssignmentMethod.manual),
            data_model.Resource(namespace='ec-code', id='1.1.1.2', relevance=30.,
                                assignment_method=data_model.ResourceAssignmentMethod.predicted),
            data_model.Resource(namespace='ec-code', id='1.1.1.3', relevance=10.,
                                assignment_method=data_model.ResourceAssignmentMethod.manual),
        ])
        self.assertRaises(ValueError, rxn.get_ec_number)

        rxn = data_model.Reaction(cross_references=[
            data_model.Resource(namespace='ec2', id='1.1.1.1', relevance=20.,
                                assignment_method=data_model.ResourceAssignmentMethod.manual),
            data_model.Resource(namespace='ec-code', id='1.1.1.2', relevance=30.,
                                assignment_method=data_model.ResourceAssignmentMethod.predicted),
            data_model.Resource(namespace='ec2', id='1.1.1.3', relevance=10.,
                                assignment_method=data_model.ResourceAssignmentMethod.manual),
        ])
        self.assertEqual(rxn.get_ec_number(), '1.1.1.2')

        rxn = data_model.Reaction(cross_references=[])
        self.assertEqual(rxn.get_ec_number(), '')


class TestGenetics(unittest.TestCase):

    def test(self):
        genetics = data_model.Genetics(taxon='Mycoplasma pneumoniae', variation='ΔMPN001')
        self.assertFalse(genetics.is_wildtype())
        self.assertTrue(genetics.is_variant())

        genetics = data_model.Genetics(taxon='Mycoplasma pneumoniae', variation='')
        self.assertTrue(genetics.is_wildtype())
        self.assertFalse(genetics.is_variant())


class TestResourceAssignmentMethod(unittest.TestCase):

    def test(self):
        self.assertEqual(data_model.ResourceAssignmentMethod['manual'].value, 0)
        self.assertEqual(data_model.ResourceAssignmentMethod['predicted'].value, 1)


class TestResource(unittest.TestCase):

    def test_init(self):
        xr = data_model.Resource(namespace='src', id='identifier', relevance=2.,
                                 assignment_method=data_model.ResourceAssignmentMethod.manual)
        self.assertEqual(xr.namespace, 'src')
        self.assertEqual(xr.id, 'identifier')
        self.assertEqual(xr.relevance, 2.)
        self.assertEqual(xr.assignment_method, data_model.ResourceAssignmentMethod.manual)
