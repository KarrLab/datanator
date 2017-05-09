""" Tests of the reaction utilities

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_structs
from kinetic_datanator.util import compartment_util
from kinetic_datanator.util import molecule_util
from kinetic_datanator.util import reaction_util
from kinetic_datanator.util import warning_util
import unittest

warning_util.disable_warnings()


class TestReaction(unittest.TestCase):
    adp = 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'
    atp = 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'
    h = '[H+]'
    h2o = 'O'
    pi = 'OP([O-])([O-])=O'

    def make_reaction(self):
        return reaction_util.Reaction(participants=[
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(structure=self.atp, id='atp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(structure=self.h2o, id='h2o'),
                compartment=compartment_util.Compartment('c'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(structure=self.adp, id='adp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(structure=self.pi, id='pi'),
                compartment=compartment_util.Compartment('c'),
                coefficient=1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(structure=self.h, id='h'),
                compartment=compartment_util.Compartment('c'),
                coefficient=1),
        ])

    def test_init(self):
        rxn = self.make_reaction()
        self.assertEqual(len(rxn.participants), 5)

    def test_id(self):
        rxn = reaction_util.Reaction(id='rxn')
        self.assertEqual(rxn.id, 'rxn')

    def test_name(self):
        rxn = reaction_util.Reaction(name='rxn')
        self.assertEqual(rxn.name, 'rxn')

    def test_normalize(self):
        rxn = reaction_util.Reaction(participants=[
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(structure=self.atp, id='atp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(structure=self.atp, id='atp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(structure=self.atp, id='atp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(structure=self.atp, id='atp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=0),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(structure=self.atp, id='atp'),
                compartment=compartment_util.Compartment('e'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(structure=self.adp, id='adp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(structure=self.adp, id='adp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=3),
        ])

        rxn.normalize()

        self.assertEqual(len(rxn.participants), 4)

        part = list(filter(lambda p: p.molecule.id == 'atp' and p.compartment.id == 'c' and p.coefficient < 0, rxn.participants))
        self.assertEqual(len(part), 1)
        self.assertEqual(part[0].coefficient, -2)

        part = list(filter(lambda p: p.molecule.id == 'atp' and p.compartment.id == 'c' and p.coefficient > 0, rxn.participants))
        self.assertEqual(len(part), 1)
        self.assertEqual(part[0].coefficient, 1)

        part = list(filter(lambda p: p.molecule.id == 'atp' and p.compartment.id == 'e', rxn.participants))
        self.assertEqual(len(part), 1)
        self.assertEqual(part[0].coefficient, -1)

        part = list(filter(lambda p: p.molecule.id == 'adp', rxn.participants))
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

        self.assertEqual(pairs[0][0].molecule.id, 'atp')
        self.assertEqual(pairs[0][1].molecule.id, 'adp')

        self.assertEqual(pairs[1][0].molecule.id, 'h2o')
        self.assertIn(pairs[1][1].molecule.id, ['h', 'pi'])

        self.assertEqual(pairs[2][0], None)
        self.assertIn(pairs[2][1].molecule.id, ['h', 'pi'])

    def test_calc_reactant_product_pairs(self):
        rxn = self.make_reaction()
        pairs = rxn.calc_reactant_product_pairs()

        self.assertEqual(pairs[0][0].molecule.id, 'atp')
        self.assertEqual(pairs[0][1].molecule.id, 'adp')

        self.assertEqual(pairs[1][0].molecule.id, 'h2o')
        self.assertIn(pairs[1][1].molecule.id, ['h', 'pi'])

        self.assertEqual(pairs[2][0], None)
        self.assertIn(pairs[2][1].molecule.id, ['h', 'pi'])

    def test_get_ec_number(self):
        rxn = reaction_util.Reaction(cross_references=[
            data_structs.CrossReference(namespace='xx', id='yy', relevance=20.,
                                        assignment_method=data_structs.CrossReferenceAssignmentMethod.predicted),
            data_structs.CrossReference(namespace='ec-code', id='1.1.1.1', relevance=20.,
                                        assignment_method=data_structs.CrossReferenceAssignmentMethod.predicted),
            data_structs.CrossReference(namespace='ec-code', id='1.1.1.2', relevance=30.,
                                        assignment_method=data_structs.CrossReferenceAssignmentMethod.predicted),
            data_structs.CrossReference(namespace='ec-code', id='1.1.1.3', relevance=10.,
                                        assignment_method=data_structs.CrossReferenceAssignmentMethod.predicted),
        ])
        self.assertEqual(rxn.get_ec_number(), '1.1.1.2')
        self.assertEqual(rxn.get_ec_numbers(), rxn.cross_references[1:])

        rxn = reaction_util.Reaction(cross_references=[
            data_structs.CrossReference(namespace='ec-code', id='1.1.1.1', relevance=20.,
                                        assignment_method=data_structs.CrossReferenceAssignmentMethod.manual),
            data_structs.CrossReference(namespace='ec-code', id='1.1.1.2', relevance=30.,
                                        assignment_method=data_structs.CrossReferenceAssignmentMethod.predicted),
            data_structs.CrossReference(namespace='ec-code', id='1.1.1.3', relevance=10.,
                                        assignment_method=data_structs.CrossReferenceAssignmentMethod.predicted),
            data_structs.CrossReference(namespace='xx', id='yy', relevance=20.,
                                        assignment_method=data_structs.CrossReferenceAssignmentMethod.predicted),
        ])
        self.assertEqual(rxn.get_ec_number(), '1.1.1.1')
        self.assertEqual(rxn.get_ec_numbers(), rxn.cross_references[0:-1])

        rxn = reaction_util.Reaction(cross_references=[
            data_structs.CrossReference(namespace='ec-code', id='1.1.1.1', relevance=20.,
                                        assignment_method=data_structs.CrossReferenceAssignmentMethod.manual),
            data_structs.CrossReference(namespace='ec-code', id='1.1.1.2', relevance=30.,
                                        assignment_method=data_structs.CrossReferenceAssignmentMethod.predicted),
            data_structs.CrossReference(namespace='ec-code', id='1.1.1.3', relevance=10.,
                                        assignment_method=data_structs.CrossReferenceAssignmentMethod.manual),
        ])
        self.assertRaises(ValueError, rxn.get_ec_number)

        rxn = reaction_util.Reaction(cross_references=[
            data_structs.CrossReference(namespace='ec2', id='1.1.1.1', relevance=20.,
                                        assignment_method=data_structs.CrossReferenceAssignmentMethod.manual),
            data_structs.CrossReference(namespace='ec-code', id='1.1.1.2', relevance=30.,
                                        assignment_method=data_structs.CrossReferenceAssignmentMethod.predicted),
            data_structs.CrossReference(namespace='ec2', id='1.1.1.3', relevance=10.,
                                        assignment_method=data_structs.CrossReferenceAssignmentMethod.manual),
        ])
        self.assertEqual(rxn.get_ec_number(), '1.1.1.2')

        rxn = reaction_util.Reaction(cross_references=[])
        self.assertEqual(rxn.get_ec_number(), '')
