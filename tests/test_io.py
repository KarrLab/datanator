""" Tests of io

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-14
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator import io
from kinetic_datanator.util import compartment_util
from kinetic_datanator.util import molecule_util
from kinetic_datanator.util import reaction_util
from kinetic_datanator.util import taxonomy_util
from os import path
from wc_utils.util.types import assert_value_equal
import unittest


class TestInputReader(unittest.TestCase):

    def test_parse_reaction_equation(self):
        eq = "[c]: Complex_Acp + ATP + HDCA ==> PPI + AMP + Complex_Acp_hdc"
        response = io.InputReader.parse_reaction_equation(eq)
        expected_answer = reaction_util.Reaction(
            participants=[
                reaction_util.ReactionParticipant(molecule='Complex_Acp',     compartment='c', coefficient=-1, order=0),
                reaction_util.ReactionParticipant(molecule='ATP',             compartment='c', coefficient=-1, order=1),
                reaction_util.ReactionParticipant(molecule='HDCA',            compartment='c', coefficient=-1, order=2),
                reaction_util.ReactionParticipant(molecule='PPI',             compartment='c', coefficient=1, order=3),
                reaction_util.ReactionParticipant(molecule='AMP',             compartment='c', coefficient=1, order=4),
                reaction_util.ReactionParticipant(molecule='Complex_Acp_hdc', compartment='c', coefficient=1, order=5),
            ],
            reversible=False,
        )
        assert_value_equal(response, expected_answer)

        eq = "A[c] + B2[e] + 2 C[c] ==> 2.5 D[c] + E[c] + F[c]"
        response = io.InputReader.parse_reaction_equation(eq)
        expected_answer = reaction_util.Reaction(
            participants=[
                reaction_util.ReactionParticipant(molecule='A',  compartment='c', coefficient=-1., order=0),
                reaction_util.ReactionParticipant(molecule='B2', compartment='e', coefficient=-1., order=1),
                reaction_util.ReactionParticipant(molecule='C',  compartment='c', coefficient=-2., order=2),
                reaction_util.ReactionParticipant(molecule='D',  compartment='c', coefficient=2.5, order=3),
                reaction_util.ReactionParticipant(molecule='E',  compartment='c', coefficient=1., order=4),
                reaction_util.ReactionParticipant(molecule='F',  compartment='c', coefficient=1., order=5),
            ],
            reversible=False,
        )
        assert_value_equal(response, expected_answer)

        eq = "A[c]+B_2[e]+2 C[c]<==>2.5 D[c] + E[c] + F[c]"
        response = io.InputReader.parse_reaction_equation(eq)
        expected_answer = reaction_util.Reaction(
            participants=[
                reaction_util.ReactionParticipant(molecule='A',   compartment='c', coefficient=-1., order=0),
                reaction_util.ReactionParticipant(molecule='B_2', compartment='e', coefficient=-1., order=1),
                reaction_util.ReactionParticipant(molecule='C',   compartment='c', coefficient=-2., order=2),
                reaction_util.ReactionParticipant(molecule='D',   compartment='c', coefficient=2.5, order=3),
                reaction_util.ReactionParticipant(molecule='E',   compartment='c', coefficient=1., order=4),
                reaction_util.ReactionParticipant(molecule='F',   compartment='c', coefficient=1., order=5),
            ],
            reversible=True,
        )
        assert_value_equal(response, expected_answer)

        # alternative separators
        expected_answer = reaction_util.Reaction(
            participants=[
                reaction_util.ReactionParticipant(molecule='A',  compartment='c', coefficient=-1., order=0),
                reaction_util.ReactionParticipant(molecule='B2', compartment='e', coefficient=-1., order=1),
                reaction_util.ReactionParticipant(molecule='C',  compartment='c', coefficient=-2., order=2),
                reaction_util.ReactionParticipant(molecule='D',  compartment='c', coefficient=2.5, order=3),
                reaction_util.ReactionParticipant(molecule='E',  compartment='c', coefficient=1., order=4),
                reaction_util.ReactionParticipant(molecule='F',  compartment='c', coefficient=1., order=5),
            ],
            reversible=False,
        )

        eq = "A[c] + B2[e] + 2 C[c] ==> 2.5 D[c] + E[c] + F[c]"
        assert_value_equal(io.InputReader.parse_reaction_equation(eq), expected_answer)

        eq = "A[c] + B2[e] + 2 C[c] => 2.5 D[c] + E[c] + F[c]"
        assert_value_equal(io.InputReader.parse_reaction_equation(eq), expected_answer)

        eq = "A[c] + B2[e] + 2 C[c] --> 2.5 D[c] + E[c] + F[c]"
        assert_value_equal(io.InputReader.parse_reaction_equation(eq), expected_answer)

        eq = "A[c] + B2[e] + 2 C[c] -> 2.5 D[c] + E[c] + F[c]"
        assert_value_equal(io.InputReader.parse_reaction_equation(eq), expected_answer)

        expected_answer = reaction_util.Reaction(
            participants=[
                reaction_util.ReactionParticipant(molecule='A',  compartment='c', coefficient=-1., order=0),
                reaction_util.ReactionParticipant(molecule='B2', compartment='e', coefficient=-1., order=1),
                reaction_util.ReactionParticipant(molecule='C',  compartment='c', coefficient=-2., order=2),
                reaction_util.ReactionParticipant(molecule='D',  compartment='c', coefficient=2.5, order=3),
                reaction_util.ReactionParticipant(molecule='E',  compartment='c', coefficient=1., order=4),
                reaction_util.ReactionParticipant(molecule='F',  compartment='c', coefficient=1., order=5),
            ],
            reversible=True,
        )

        eq = "A[c] + B2[e] + 2 C[c] <==> 2.5 D[c] + E[c] + F[c]"
        assert_value_equal(io.InputReader.parse_reaction_equation(eq), expected_answer)

        eq = "A[c] + B2[e] + 2 C[c] <=> 2.5 D[c] + E[c] + F[c]"
        assert_value_equal(io.InputReader.parse_reaction_equation(eq), expected_answer)

        eq = "A[c] + B2[e] + 2 C[c] <--> 2.5 D[c] + E[c] + F[c]"
        assert_value_equal(io.InputReader.parse_reaction_equation(eq), expected_answer)

        eq = "A[c] + B2[e] + 2 C[c] <-> 2.5 D[c] + E[c] + F[c]"
        assert_value_equal(io.InputReader.parse_reaction_equation(eq), expected_answer)

        # negative examples
        eq = "[c]: A[c] + B2[e] + 2 C[c] ==> 2.5 D[c] + E[c] + F[c]"
        self.assertRaises(ValueError, io.InputReader.parse_reaction_equation, eq)

        eq = "[c]: A + B2 + 2 C[c] ==> 2.5 D + E + F"
        self.assertRaises(ValueError, io.InputReader.parse_reaction_equation, eq)

        eq = "A + B2 + 2 C[c] ==> 2.5 D + E + F"
        self.assertRaises(ValueError, io.InputReader.parse_reaction_equation, eq)

        eq = "A[c] + B-2[e] + 2 C[c] ==> 2.5 D[c] + E[c] + F[c]"
        self.assertRaises(ValueError, io.InputReader.parse_reaction_equation, eq)

        eq = "A[c] + B2[e] + 2 C[c] ===> 2.5 D[c] + E[c] + F[c]"
        self.assertRaises(ValueError, io.InputReader.parse_reaction_equation, eq)

        eq = "A[c] + B2[e] + 2 C[c] ---> 2.5 D[c] + E[c] + F[c]"
        self.assertRaises(ValueError, io.InputReader.parse_reaction_equation, eq)

    def test_run(self):
        filename = path.join(path.dirname(__file__), "fixtures", "five_reactions.xlsx")
        taxon, compartments, molecules, reactions = io.InputReader.run(filename)

        #####################################
        # taxon
        #####################################
        self.assertIsInstance(taxon, taxonomy_util.Taxon)
        self.assertEqual(taxon.name, 'Mycoplasma pneumoniae M129')

        #####################################
        # compartments
        #####################################
        self.assertEqual(len(compartments), 1)
        assert_value_equal(compartments[0], compartment_util.Compartment(id='c'))

        #####################################
        # molecules
        #####################################
        self.assertEqual(len(molecules), 17)

        assert_value_equal(molecules[0], molecule_util.Molecule(
            id='ADP',
            structure='C1=NC2=C(C(=N1)N)N=CN2C3C(C(C(O3)COP(=O)(O)OP(=O)(O)O)O)O',
            cross_references=[],
        ))

        assert_value_equal(molecules[3], molecule_util.Molecule(
            id='dCMP64dTMP',
            structure='CC1=CN(C2CC(O)C(COP([O-])([O-])=O)O2)C(=O)NC1C1=[N+](C2CC(O)C(COP([O-])([O-])=O)O2)C(=O)N=C(N)[C-]1O',
            cross_references=[],
        ))

        assert_value_equal(molecules[12], molecule_util.Molecule(
            id='NAD',
            structure='NC(=O)C1=CC=C[N+](=C1)C1OC(COP([O-])(=O)OP([O-])(=O)OCC2OC(C(O)C2O)N2C=NC3=C(N)N=CN=C23)C(O)C1O',
            cross_references=[],
        ))

        #####################################
        # reactions
        #####################################
        c = next(c for c in compartments if c.id == 'c')
        atp = next(m for m in molecules if m.id == 'ATP')
        ump = next(m for m in molecules if m.id == 'UMP')
        udp = next(m for m in molecules if m.id == 'UDP')
        adp = next(m for m in molecules if m.id == 'ADP')

        self.assertEqual(len(reactions), 5)
        assert_value_equal(reactions[0], reaction_util.Reaction(id='ump_kinase', reversible=True, participants=[
            reaction_util.ReactionParticipant(molecule=ump, compartment=c, coefficient=-1, order=0),
            reaction_util.ReactionParticipant(molecule=atp, compartment=c, coefficient=-1, order=1),
            reaction_util.ReactionParticipant(molecule=udp, compartment=c, coefficient=1, order=2),
            reaction_util.ReactionParticipant(molecule=adp, compartment=c, coefficient=1, order=3),
        ]))


class TestResultsWriter(unittest.TestCase):
    # todo
    pass
