""" Tests of the reaction utilities

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.util import compartment_util
from kinetic_datanator.util import molecule_util
from kinetic_datanator.util import reaction_util
import unittest


class TestReactionUtil(unittest.TestCase):
    adp = 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'
    atp = 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'
    h = '[H+]'
    h2o = 'O'
    pi = 'OP([O-])([O-])=O'

    def make_reaction(self):
        return reaction_util.Reaction([
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.atp, name='atp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.h2o, name='h2o'),
                compartment=compartment_util.Compartment('c'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.adp, name='adp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.pi, name='pi'),
                compartment=compartment_util.Compartment('c'),
                coefficient=1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.h, name='h'),
                compartment=compartment_util.Compartment('c'),
                coefficient=1),
        ])

    def test_Reaction(self):
        rxn = self.make_reaction()
        self.assertEqual(len(rxn.participants), 5)

    def test_Reaction_name(self):
        rxn = reaction_util.Reaction([], name='rxn')
        self.assertEqual(rxn.name, 'rxn')

    def test_Reaction_normalize(self):
        rxn = reaction_util.Reaction([
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.atp, name='atp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.atp, name='atp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.atp, name='atp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.atp, name='atp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=0),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.atp, name='atp'),
                compartment=compartment_util.Compartment('e'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.adp, name='adp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.adp, name='adp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=3),
        ])

        rxn.normalize()

        self.assertEqual(len(rxn.participants), 4)

        part = filter(lambda p: p.molecule.name == 'atp' and p.compartment.name == 'c' and p.coefficient < 0, rxn.participants)
        self.assertEqual(len(part), 1)
        self.assertEqual(part[0].coefficient, -2)

        part = filter(lambda p: p.molecule.name == 'atp' and p.compartment.name == 'c' and p.coefficient > 0, rxn.participants)
        self.assertEqual(len(part), 1)
        self.assertEqual(part[0].coefficient, 1)

        part = filter(lambda p: p.molecule.name == 'atp' and p.compartment.name == 'e', rxn.participants)
        self.assertEqual(len(part), 1)
        self.assertEqual(part[0].coefficient, -1)

        part = filter(lambda p: p.molecule.name == 'adp', rxn.participants)
        self.assertEqual(len(part), 1)
        self.assertEqual(part[0].coefficient, 4)

    def test_Reaction_get_reactants(self):
        rxn = self.make_reaction()
        self.assertEqual(rxn.get_reactants(), rxn.participants[0:2])

    def test_Reaction_get_products(self):
        rxn = self.make_reaction()
        self.assertEqual(rxn.get_products(), rxn.participants[2:])
