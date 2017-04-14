""" Tests of io

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-14
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator import io
import unittest

class TestIo(unittest.TestCase):

    def test_parse_reaction_stoichiometry(self):
        reaction_string = "[c]: Complex_Acp + ATP + HDCA ==> PPI + AMP + Complex_Acp_hdc"
        response = io.InputReader.parse_reaction_stoichiometry(reaction_string)
        expectedAnswer = [['Complex_Acp', 'ATP', 'HDCA'], ['PPI', 'AMP', 'Complex_Acp_hdc']]
        self.assertEqual(response, expectedAnswer)

        reaction_string = "A + B + C ==> D + E + F"
        response = io.InputReader.parse_reaction_stoichiometry(reaction_string)
        expectedAnswer = [['A', 'B', 'C'], ['D', 'E', 'F']]
        self.assertEqual(response, expectedAnswer)

        # make sure it removes hyrdrogens
        reaction_string = "[c]: Complex_Acp + H2O ==> ProtMon_MPN406 + H + Pantetheine4Phosphate"
        response = io.InputReader.parse_reaction_stoichiometry(reaction_string)
        expectedAnswer = [['Complex_Acp', 'H2O'], ['ProtMon_MPN406', 'Pantetheine4Phosphate']]
        self.assertEqual(response, expectedAnswer)

        reaction_string = "cyclobutane_dTdC[c] ==> cyclobutane_dTdC[e]"
        response = io.InputReader.parse_reaction_stoichiometry(reaction_string)
        expectedAnswer = [['cyclobutane_dTdC'], ['cyclobutane_dTdC']]
        self.assertEqual(response, expectedAnswer)

        reaction_string = "UDPG[c] + DAG161[m] ==> UDP[c] + H[c] + GlcDAG161[m]"
        response = io.InputReader.parse_reaction_stoichiometry(reaction_string)
        expectedAnswer = [['UDPG', 'DAG161'], ['UDP', 'GlcDAG161']]
        self.assertEqual(response, expectedAnswer)

        reaction_string = 'h_m + 1a25dhvitd2_m + o2_m + nadph_m --> h2o_m + nadp_m + 1a2425thvitd2_m'
        response = io.InputReader.parse_reaction_stoichiometry(reaction_string)
        expectedAnswer = [['1a25dhvitd2_m', 'o2_m', 'nadph_m'], ['h2o_m', 'nadp_m', '1a2425thvitd2_m']]
        self.assertEqual(response, expectedAnswer)

        reaction_string = 'h_m + 1a25dhvitd2_m + o2_m + nadph_m <-> h2o_m + nadp_m + 1a2425thvitd2_m'
        response = io.InputReader.parse_reaction_stoichiometry(reaction_string)
        expectedAnswer = [['1a25dhvitd2_m', 'o2_m', 'nadph_m'], ['h2o_m', 'nadp_m', '1a2425thvitd2_m']]
        self.assertEqual(response, expectedAnswer)

        reaction_string = 'h_m + 1a25dhvitd2_m + o2_m + nadph_m <=> h2o_m + nadp_m + 1a2425thvitd2_m'
        response = io.InputReader.parse_reaction_stoichiometry(reaction_string)
        expectedAnswer = [['1a25dhvitd2_m', 'o2_m', 'nadph_m'], ['h2o_m', 'nadp_m', '1a2425thvitd2_m']]
        self.assertEqual(response, expectedAnswer)
