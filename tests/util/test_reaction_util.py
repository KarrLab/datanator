""" Tests of the reaction utilities

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model
from kinetic_datanator.util import reaction_util
import unittest


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

    def test_calc_reactant_product_pairs(self):
        rxn = self.make_reaction()
        pairs = reaction_util.calc_reactant_product_pairs(rxn)

        self.assertEqual(pairs[0][0].specie.id, 'atp')
        self.assertEqual(pairs[0][1].specie.id, 'adp')

        self.assertEqual(pairs[1][0].specie.id, 'h2o')
        self.assertIn(pairs[1][1].specie.id, ['h', 'pi'])

        self.assertEqual(pairs[2][0], None)
        self.assertIn(pairs[2][1].specie.id, ['h', 'pi'])
