""" Tests of io

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-14
:Copyright: 2017, Karr Lab
:License: MIT
"""

from datanator import io
from datanator.core import data_model
from os import path
from wc_utils.util.types import assert_value_equal
import unittest


class TestInputReader(unittest.TestCase):

    def are_reactions_equal(self, rxn1, rxn2):
        if rxn1.reversible != rxn2.reversible:
            return False

        for part1, part2 in zip(rxn1.participants, rxn2.participants):
            if part1.specie != part2.specie:
                return False
            if part1.compartment != part2.compartment:
                return False
            if part1.coefficient != part2.coefficient:
                return False
            if part1.order != part2.order:
                return False

        return True

    def test_parse_reaction_equation(self):
        c = data_model.Compartment(id='c')
        e = data_model.Compartment(id='e')
        compartments = [c, e]

        Complex_Acp = data_model.Specie(id='Complex_Acp')
        ATP = data_model.Specie(id='ATP')
        HDCA = data_model.Specie(id='HDCA')
        PPI = data_model.Specie(id='PPI')
        AMP = data_model.Specie(id='AMP')
        Complex_Acp_hdc = data_model.Specie(id='Complex_Acp_hdc')
        A = data_model.Specie(id='A')
        B2 = data_model.Specie(id='B2')
        B_2 = data_model.Specie(id='B_2')
        C = data_model.Specie(id='C')
        D = data_model.Specie(id='D')
        E = data_model.Specie(id='E')
        F = data_model.Specie(id='F')
        species = [
            Complex_Acp,
            ATP,
            HDCA,
            PPI,
            AMP,
            Complex_Acp_hdc,
            A,
            B2,
            B_2,
            C,
            D,
            E,
            F,
        ]

        eq = "[c]: Complex_Acp + ATP + HDCA ==> PPI + AMP + Complex_Acp_hdc"
        response = io.InputReader().parse_reaction_equation(eq, compartments, species)
        expected_answer = data_model.Reaction(
            participants=[
                data_model.ReactionParticipant(specie=Complex_Acp,     compartment=c, coefficient=-1, order=0),
                data_model.ReactionParticipant(specie=ATP,             compartment=c, coefficient=-1, order=1),
                data_model.ReactionParticipant(specie=HDCA,            compartment=c, coefficient=-1, order=2),
                data_model.ReactionParticipant(specie=PPI,             compartment=c, coefficient=1, order=3),
                data_model.ReactionParticipant(specie=AMP,             compartment=c, coefficient=1, order=4),
                data_model.ReactionParticipant(specie=Complex_Acp_hdc, compartment=c, coefficient=1, order=5),
            ],
            reversible=False,
        )
        self.assertTrue(self.are_reactions_equal(response, expected_answer))

        eq = "A[c] + B2[e] + 2 C[c] ==> 2.5 D[c] + E[c] + F[c]"
        response = io.InputReader().parse_reaction_equation(eq, compartments, species)
        expected_answer = data_model.Reaction(
            participants=[
                data_model.ReactionParticipant(specie=A,  compartment=c, coefficient=-1., order=0),
                data_model.ReactionParticipant(specie=B2, compartment=e, coefficient=-1., order=1),
                data_model.ReactionParticipant(specie=C,  compartment=c, coefficient=-2., order=2),
                data_model.ReactionParticipant(specie=D,  compartment=c, coefficient=2.5, order=3),
                data_model.ReactionParticipant(specie=E,  compartment=c, coefficient=1., order=4),
                data_model.ReactionParticipant(specie=F,  compartment=c, coefficient=1., order=5),
            ],
            reversible=False,
        )
        self.assertTrue(self.are_reactions_equal(response, expected_answer))

        eq = "A[c]+B_2[e]+2 C[c]<==>2.5 D[c] + E[c] + F[c]"
        response = io.InputReader().parse_reaction_equation(eq, compartments, species)
        expected_answer = data_model.Reaction(
            participants=[
                data_model.ReactionParticipant(specie=A,   compartment=c, coefficient=-1., order=0),
                data_model.ReactionParticipant(specie=B_2, compartment=e, coefficient=-1., order=1),
                data_model.ReactionParticipant(specie=C,   compartment=c, coefficient=-2., order=2),
                data_model.ReactionParticipant(specie=D,   compartment=c, coefficient=2.5, order=3),
                data_model.ReactionParticipant(specie=E,   compartment=c, coefficient=1., order=4),
                data_model.ReactionParticipant(specie=F,   compartment=c, coefficient=1., order=5),
            ],
            reversible=True,
        )
        self.assertTrue(self.are_reactions_equal(response, expected_answer))

        # alternative separators
        expected_answer = data_model.Reaction(
            participants=[
                data_model.ReactionParticipant(specie=A,  compartment=c, coefficient=-1., order=0),
                data_model.ReactionParticipant(specie=B2, compartment=e, coefficient=-1., order=1),
                data_model.ReactionParticipant(specie=C,  compartment=c, coefficient=-2., order=2),
                data_model.ReactionParticipant(specie=D,  compartment=c, coefficient=2.5, order=3),
                data_model.ReactionParticipant(specie=E,  compartment=c, coefficient=1., order=4),
                data_model.ReactionParticipant(specie=F,  compartment=c, coefficient=1., order=5),
            ],
            reversible=False,
        )

        eq = "A[c] + B2[e] + 2 C[c] ==> 2.5 D[c] + E[c] + F[c]"
        self.are_reactions_equal(io.InputReader().parse_reaction_equation(eq, compartments, species), expected_answer)

        eq = "A[c] + B2[e] + 2 C[c] => 2.5 D[c] + E[c] + F[c]"
        self.are_reactions_equal(io.InputReader().parse_reaction_equation(eq, compartments, species), expected_answer)

        eq = "A[c] + B2[e] + 2 C[c] --> 2.5 D[c] + E[c] + F[c]"
        self.are_reactions_equal(io.InputReader().parse_reaction_equation(eq, compartments, species), expected_answer)

        eq = "A[c] + B2[e] + 2 C[c] -> 2.5 D[c] + E[c] + F[c]"
        self.are_reactions_equal(io.InputReader().parse_reaction_equation(eq, compartments, species), expected_answer)

        expected_answer = data_model.Reaction(
            participants=[
                data_model.ReactionParticipant(specie=A,  compartment=c, coefficient=-1., order=0),
                data_model.ReactionParticipant(specie=B2, compartment=e, coefficient=-1., order=1),
                data_model.ReactionParticipant(specie=C,  compartment=c, coefficient=-2., order=2),
                data_model.ReactionParticipant(specie=D,  compartment=c, coefficient=2.5, order=3),
                data_model.ReactionParticipant(specie=E,  compartment=c, coefficient=1., order=4),
                data_model.ReactionParticipant(specie=F,  compartment=c, coefficient=1., order=5),
            ],
            reversible=True,
        )

        eq = "A[c] + B2[e] + 2 C[c] <==> 2.5 D[c] + E[c] + F[c]"
        self.are_reactions_equal(io.InputReader().parse_reaction_equation(eq, compartments, species), expected_answer)

        eq = "A[c] + B2[e] + 2 C[c] <=> 2.5 D[c] + E[c] + F[c]"
        self.are_reactions_equal(io.InputReader().parse_reaction_equation(eq, compartments, species), expected_answer)

        eq = "A[c] + B2[e] + 2 C[c] <--> 2.5 D[c] + E[c] + F[c]"
        self.are_reactions_equal(io.InputReader().parse_reaction_equation(eq, compartments, species), expected_answer)

        eq = "A[c] + B2[e] + 2 C[c] <-> 2.5 D[c] + E[c] + F[c]"
        self.are_reactions_equal(io.InputReader().parse_reaction_equation(eq, compartments, species), expected_answer)

        # negative examples
        eq = "[c]: A[c] + B2[e] + 2 C[c] ==> 2.5 D[c] + E[c] + F[c]"
        self.assertRaises(ValueError, io.InputReader().parse_reaction_equation, eq, compartments, species)

        eq = "[c]: A + B2 + 2 C[c] ==> 2.5 D + E + F"
        self.assertRaises(ValueError, io.InputReader().parse_reaction_equation, eq, compartments, species)

        eq = "A + B2 + 2 C[c] ==> 2.5 D + E + F"
        self.assertRaises(ValueError, io.InputReader().parse_reaction_equation, eq, compartments, species)

        eq = "A[c] + B-2[e] + 2 C[c] ==> 2.5 D[c] + E[c] + F[c]"
        self.assertRaises(ValueError, io.InputReader().parse_reaction_equation, eq, compartments, species)

        eq = "A[c] + B2[e] + 2 C[c] ===> 2.5 D[c] + E[c] + F[c]"
        self.assertRaises(ValueError, io.InputReader().parse_reaction_equation, eq, compartments, species)

        eq = "A[c] + B2[e] + 2 C[c] ---> 2.5 D[c] + E[c] + F[c]"
        self.assertRaises(ValueError, io.InputReader().parse_reaction_equation, eq, compartments, species)

    def test_run(self):
        filename = path.join(path.dirname(__file__), "fixtures", "five_reactions.xlsx")
        genetics, compartments, species, reactions = io.InputReader().run(filename)

        #####################################
        # taxon
        #####################################
        self.assertIsInstance(genetics, data_model.Genetics)
        self.assertEqual(genetics.taxon, 'Mycoplasma pneumoniae M129')

        #####################################
        # compartments
        #####################################
        self.assertEqual(len(compartments), 1)
        assert_value_equal(compartments[0].id, 'c')

        #####################################
        # species
        #####################################
        self.assertEqual(len(species), 17)

        self.assertEqual(species[0].id, 'ADP')
        self.assertEqual(species[0].structure, 'C1=NC2=C(C(=N1)N)N=CN2C3C(C(C(O3)COP(=O)(O)OP(=O)(O)O)O)O')
        self.assertEqual(species[0].cross_references, [])

        self.assertEqual(species[3].id, 'dCMP64dTMP')
        self.assertEqual(species[3].structure,
                         'CC1=CN(C2CC(O)C(COP([O-])([O-])=O)O2)C(=O)NC1C1=[N+](C2CC(O)C(COP([O-])([O-])=O)O2)C(=O)N=C(N)[C-]1O')
        self.assertEqual(species[3].cross_references, [])

        self.assertEqual(species[12].id, 'NAD')
        self.assertEqual(species[12].structure,
                         'NC(=O)C1=CC=C[N+](=C1)C1OC(COP([O-])(=O)OP([O-])(=O)OCC2OC(C(O)C2O)N2C=NC3=C(N)N=CN=C23)C(O)C1O')
        self.assertEqual(species[12].cross_references, [])

        #####################################
        # reactions
        #####################################
        c = next(c for c in compartments if c.id == 'c')
        atp = next(s for s in species if s.id == 'ATP')
        ump = next(s for s in species if s.id == 'UMP')
        udp = next(s for s in species if s.id == 'UDP')
        adp = next(s for s in species if s.id == 'ADP')

        self.assertEqual(len(reactions), 5)
        self.assertEqual(reactions[0].id, 'ump_kinase')
        self.assertEqual(reactions[0].reversible, True)

        self.assertEqual(reactions[0].participants[0].specie, ump)
        self.assertEqual(reactions[0].participants[1].specie, atp)
        self.assertEqual(reactions[0].participants[2].specie, udp)
        self.assertEqual(reactions[0].participants[3].specie, adp)

        self.assertEqual(reactions[0].participants[0].compartment, c)
        self.assertEqual(reactions[0].participants[1].compartment, c)
        self.assertEqual(reactions[0].participants[2].compartment, c)
        self.assertEqual(reactions[0].participants[3].compartment, c)

        self.assertEqual(reactions[0].participants[0].coefficient, -1)
        self.assertEqual(reactions[0].participants[1].coefficient, -1)
        self.assertEqual(reactions[0].participants[2].coefficient, 1)
        self.assertEqual(reactions[0].participants[3].coefficient, 1)

        self.assertEqual(reactions[0].participants[0].order, 0)
        self.assertEqual(reactions[0].participants[1].order, 1)
        self.assertEqual(reactions[0].participants[2].order, 2)
        self.assertEqual(reactions[0].participants[3].order, 3)


class TestResultsWriter(unittest.TestCase):
    # todo: implement
    pass
