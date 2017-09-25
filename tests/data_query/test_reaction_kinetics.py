# -*- coding: utf-8 -*-

""" Tests of sabio_rk

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-26
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model
from kinetic_datanator.data_source import sabio_rk
from kinetic_datanator.core import common_schema
from kinetic_datanator.data_query import reaction_kinetics
from kinetic_datanator.util import molecule_util
from kinetic_datanator.util import taxonomy_util
import unittest


class TestwPartialCommonSchemaReactionKineticsQueryGenerator(unittest.TestCase):
    """
    Tests for 10000 entry limited Common Schema on Karr Lab Server

    Used for development purposes
    """

    def setUp(self):
        self.reaction = data_model.Reaction(
            participants = [
                data_model.ReactionParticipant(
                    specie = data_model.Specie(
                        id = 'Dihydrofolate',
                        structure = 'InChI=1S/C19H21N7O6/c20-19-25-15-14(17(30)26-19)23-11(8-22-15)'
                        '7-21-10-3-1-9(2-4-10)16(29)24-12(18(31)32)5-6-13(27)28/h1-4,12,21H,5-8H2,'
                        '(H,24,29)(H,27,28)(H,31,32)(H4,20,22,25,26,30)'),
                    coefficient = -1),
                data_model.ReactionParticipant(
                    specie = data_model.Specie(
                        id = 'NADPH',
                        structure = 'InChI=1S/C21H30N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16'
                        '(44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)40-5-10-13(29)15(31)'
                        '20(42-10)27-3-1-2-9(4-27)18(23)32/h1,3-4,7-8,10-11,13-16,20-21,29-31H,2,5-6H2,'
                        '(H2,23,32)(H,36,37)(H,38,39)(H2,22,24,25)(H2,33,34,35)'),
                    coefficient = -1),
                data_model.ReactionParticipant(
                    specie = data_model.Specie(
                        id = 'H+',
                        structure = 'InChI=1S/H'),
                    coefficient = -1),
                data_model.ReactionParticipant(
                    specie = data_model.Specie(
                        id = 'NADP+',
                        structure = 'InChI=1S/C21H28N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16('
                        '44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)40-5-10-13(29)15(31)20'
                        '(42-10)27-3-1-2-9(4-27)18(23)32/h1-4,7-8,10-11,13-16,20-21,29-31H,5-6H2,(H7-,22'
                        ',23,24,25,32,33,34,35,36,37,38,39)/p+1'),
                    coefficient = 1),
                data_model.ReactionParticipant(
                    specie = data_model.Specie(
                        id = '5,6,7,8-Tetrahydrofolate',
                        structure = 'InChI=1S/C19H23N7O6/c20-19-25-15-14(17(30)26-19)23-11(8-22-15)7-21-'
                        '10-3-1-9(2-4-10)16(29)24-12(18(31)32)5-6-13(27)28/h1-4,11-12,21,23H,5-8H2,(H,24,29'
                        ')(H,27,28)(H,31,32)(H4,20,22,25,26,30)'),
                    coefficient = 1)
        ])

        self.reaction_w_resource = data_model.Reaction(
            participants = [
                data_model.ReactionParticipant(
                    specie = data_model.Specie(
                        id = 'Dihydrofolate',
                        structure = 'InChI=1S/C19H21N7O6/c20-19-25-15-14(17(30)26-19)23-11(8-22-15)'
                        '7-21-10-3-1-9(2-4-10)16(29)24-12(18(31)32)5-6-13(27)28/h1-4,12,21H,5-8H2,'
                        '(H,24,29)(H,27,28)(H,31,32)(H4,20,22,25,26,30)'),
                    coefficient = -1),
                data_model.ReactionParticipant(
                    specie = data_model.Specie(
                        id = 'NADPH',
                        structure = 'InChI=1S/C21H30N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16'
                        '(44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)40-5-10-13(29)15(31)'
                        '20(42-10)27-3-1-2-9(4-27)18(23)32/h1,3-4,7-8,10-11,13-16,20-21,29-31H,2,5-6H2,'
                        '(H2,23,32)(H,36,37)(H,38,39)(H2,22,24,25)(H2,33,34,35)'),
                    coefficient = -1),
                data_model.ReactionParticipant(
                    specie = data_model.Specie(
                        id = 'H+',
                        structure = 'InChI=1S/H'),
                    coefficient = -1),
                data_model.ReactionParticipant(
                    specie = data_model.Specie(
                        id = 'NADP+',
                        structure = 'InChI=1S/C21H28N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16('
                        '44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)40-5-10-13(29)15(31)20'
                        '(42-10)27-3-1-2-9(4-27)18(23)32/h1-4,7-8,10-11,13-16,20-21,29-31H,5-6H2,(H7-,22'
                        ',23,24,25,32,33,34,35,36,37,38,39)/p+1'),
                    coefficient = 1),
                data_model.ReactionParticipant(
                    specie = data_model.Specie(
                        id = '5,6,7,8-Tetrahydrofolate',
                        structure = 'InChI=1S/C19H23N7O6/c20-19-25-15-14(17(30)26-19)23-11(8-22-15)7-21-'
                        '10-3-1-9(2-4-10)16(29)24-12(18(31)32)5-6-13(27)28/h1-4,11-12,21,23H,5-8H2,(H,24,29'
                        ')(H,27,28)(H,31,32)(H4,20,22,25,26,30)'),
                    coefficient = 1)
        ],
            cross_references=[
                data_model.Resource(namespace='ec-code', id='1.5.1.3', assignment_method=data_model.ResourceAssignmentMethod.manual),
        ])


    def test_get_kinetic_laws_by_ec_numbers(self):
        q = reaction_kinetics.ReactionKineticsQueryGenerator()

        # single EC, match_levels=4
        laws_62 = q.get_kinetic_laws_by_ec_numbers(['3.4.21.62'], match_levels=4).all()
        ids_62 = set([1987, 1988, 1989, 1990, 1991])
        self.assertEqual(set(l.id for l in laws_62), ids_62)

        laws_73 = q.get_kinetic_laws_by_ec_numbers(['3.4.21.73'], match_levels=4).all()
        ids_73 = set([1996, 1997])
        self.assertEqual(set(l.id for l in laws_73), ids_73)

        laws_89 = q.get_kinetic_laws_by_ec_numbers(['3.4.21.89'], match_levels=4).all()
        ids_89 = set([2112, 2113, 2114, 2115, 2116, 2117, 2118, 2119, 2120, 2121, 2122, 2123, 2106, 2107, 2108, 2109, 2110, 2111])
        self.assertEqual(set(l.id for l in laws_89), ids_89)

        #multiple EC, match_levels=4

        laws = q.get_kinetic_laws_by_ec_numbers(['3.4.21.62', '3.4.21.73'], match_levels=4).all()
        self.assertEqual(set([l.id for l in laws]), ids_62 | ids_73)

        laws = q.get_kinetic_laws_by_ec_numbers(['3.4.21.62', '3.4.21.89'], match_levels=4).all()
        self.assertEqual(set([l.id for l in laws]), ids_62 | ids_89)

        laws = q.get_kinetic_laws_by_ec_numbers(['3.4.21.73', '3.4.21.89'], match_levels=4).all()
        self.assertEqual(set([l.id for l in laws]), ids_73 | ids_89)

        laws = q.get_kinetic_laws_by_ec_numbers(['3.4.21.62', '3.4.21.73', '3.4.21.89'], match_levels=4).all()
        self.assertEqual(set([l.id for l in laws]), ids_62 | ids_73 | ids_89)


        # match_levels=3
        laws = q.get_kinetic_laws_by_ec_numbers(['3.4.21.62'], match_levels=3).all()
        ids = set([l.id for l in laws])
        ids_contains = ids_62 | ids_73 | ids_89
        self.assertEqual(ids_contains.difference(ids), set())



    def test_get_kinetic_laws_by_compound(self):
        q = reaction_kinetics.ReactionKineticsQueryGenerator()

        inchi = "InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2"


        law = q.get_kinetic_laws_by_compound(inchi, role = 'product')
        # print set(c.kineticlaw_id for c in law.all())
        self.assertEqual(set(c.kineticlaw_id for c in law.all()), set([2332,2333]))


    def test_get_compounds_by_structure(self):
        q = reaction_kinetics.ReactionKineticsQueryGenerator()

        inchi = 'InChI=1S/C5H8O3/c1-3(2)4(6)5(7)8/h3H,1-2H3,(H,7,8)'

        compounds = q.get_compounds_by_structure(inchi).all()
        self.assertEqual(set(c.id for c in compounds), set([1962]))


    def test_get_kinetic_laws_by_participants(self):
        q = reaction_kinetics.ReactionKineticsQueryGenerator()

        participants = self.reaction.participants

        law = q.get_kinetic_laws_by_participants(participants, only_formula_and_connectivity = False).order_by(common_schema.KineticLaw.kineticlaw_id)

        self.assertEqual([l.kineticlaw_id for l in law], [2192,2193,2194,2195,2196,2197,2204,2205,2206,
            2207,2208,2209,2210,2211,2212,2213,2214,2215,2216,2217,2218,2219,2220,2221,2222,2223,2224,2225,2226,
            2227,2228,2229,2230,2231,2232,2233,2234,2235,2236,2237, 2238,2239,2240,2241,2242,2243,2244,2245,2246,
            2247,2248,2249,2250,2251,2252,2253,2254,2255])


    def test_get_kinetic_laws_by_reaction(self):
        q = reaction_kinetics.ReactionKineticsQueryGenerator()

        laws = q.get_kinetic_laws_by_reaction(self.reaction) \
            .order_by(common_schema.KineticLaw.kineticlaw_id).all()

        self.assertEqual([l.kineticlaw_id for l in laws], [2192,2193,2194,2195,2196,2197,2204,2205,2206,
            2207,2208,2209,2210,2211,2212,2213,2214,2215,2216,2217,2218,2219,2220,2221,2222,2223,2224,2225,2226,
            2227,2228,2229,2230,2231,2232,2233,2234,2235,2236,2237, 2238,2239,2240,2241,2242,2243,2244,2245,2246,
            2247,2248,2249,2250,2251,2252,2253,2254,2255])

        laws = q.get_kinetic_laws_by_reaction(self.reaction_w_resource) \
            .order_by(common_schema.KineticLaw.kineticlaw_id).all()

        self.assertEqual([l.kineticlaw_id for l in laws], [2192,2193,2194,2195,2196,2197,2204,2205,2206,
            2207,2208,2209,2210,2211,2212,2213,2214,2215,2216,2217,2218,2219,2220,2221,2222,2223,2224,2225,2226,
            2227,2228,2229,2230,2231,2232,2233,2234,2235,2236,2237, 2238,2239,2240,2241,2242,2243,2244,2245,2246,
            2247,2248,2249,2250,2251,2252,2253,2254,2255])

    def test_get_observed_values(self):
        q = reaction_kinetics.ReactionKineticsQueryGenerator()
        vals = q.get_observed_values(self.reaction)

        for val in vals:
            sabiork_id = next(xr.id for xr in val.observable.interaction.cross_references if xr.namespace == 'sabiork.reaction')
            common_schema_kinetic_id = next(xr.id for xr in val.observable.interaction.cross_references if xr.namespace == 'common_schema.kineticlaw_id')
            if sabiork_id == '100' and common_schema_kinetic_id == '2196':
                if val.observable.property == 'Km_A':
                    self.assertEqual(val.observable.specie.name,'NADPH')
                    self.assertEqual(val.value, 7.3e-06)
                    self.assertEqual(val.error, 2.0e-07)
                    self.assertEqual(val.units, 'M')
                    break


        # self.print_observed_values(vals)

    @unittest.skip('implement me')
    def test_filter_observed_values(self):
        q = reaction_kinetics.ReactionKineticsQueryGenerator()
        reaction = self.reaction
        vals = q.get_observed_values(reaction)
        filter_result = q.filter_observed_values(reaction, vals)

        for ov in filter_result.observed_values:
            self.assertIn(ov.observable.property, ['k_cat', 'k_i', 'k_m', 'v_max'])
            if ov.observable.specie:
                self.assertIn(ov.observable.specie.to_inchi(only_formula_and_connectivity=True), [
                    'C3H6O2/c1-3(5)2-4',
                    ('C21H30N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47'
                     '(36,37)40-5-10-13(29)15(31)20(42-10)27-3-1-2-9(4-27)18(23)32'),
                ])
            self.assertEqual(ov.observation.genetics.variation, '')

        # self.print_observed_values(filter_result.observed_values)

    @unittest.skip('implement me')
    def test_get_consensus(self):
        mol = data_model.Specie(structure=inchi)

        ec_numbers = ['1.1.1.52', '1.1.1.55']
        rxn = data_model.Reaction(cross_references=[
            data_model.Resource(namespace='ec-code', id=ec_number),
        ])

### OLD TESTS ####
# class TestwFullCommonSchemaReactionKineticsQueryGenerator(unittest.TestCase):
#
#     def setUp(self):
#         self.reaction_1_1_1_55 = data_model.Reaction(
#             participants=[
#                 data_model.ReactionParticipant(
#                     specie=data_model.Specie(
#                         id='Lactaldehyde',
#                         structure='InChI=1S/C3H6O2/c1-3(5)2-4/h2-3,5H,1H3'),
#                     coefficient=-1),
#                 data_model.ReactionParticipant(
#                     specie=data_model.Specie(
#                         id='NADPH',
#                         structure=(
#                             'InChI=1S/C21H30N7O17P3'
#                             '/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)'
#                             '40-5-10-13(29)15(31)20(42-10)27-3-1-2-9(4-27)18(23)32'
#                             '/h1,3-4,7-8,10-11,13-16,20-21,29-31H,2,5-6H2,(H2,23,32)(H,36,37)(H,38,39)(H2,22,24,25)(H2,33,34,35)'
#                             '/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1')),
#                     coefficient=-1),
#                 data_model.ReactionParticipant(
#                     specie=data_model.Specie(
#                         id='H+',
#                         structure='InChI=1S/p+1'),
#                     coefficient=-1),
#                 data_model.ReactionParticipant(
#                     specie=data_model.Specie(
#                         id='1,2-Propanediol',
#                         structure='InChI=1S/C3H8O2/c1-3(5)2-4/h3-5H,2H2,1H3'),
#                     coefficient=1),
#                 data_model.ReactionParticipant(
#                     specie=data_model.Specie(
#                         id='NADP+',
#                         structure=(
#                             'InChI=1S/C21H28N7O17P3'
#                             '/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)'
#                             '45-47(36,37)40-5-10-13(29)15(31)20(42-10)27-3-1-2-9(4-27)18(23)32'
#                             '/h1-4,7-8,10-11,13-16,20-21,29-31H,5-6H2,(H7-,22,23,24,25,32,33,34,35,36,37,38,39)'
#                             '/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1')),
#                     coefficient=1),
#             ])
#
#         self.reaction_1_1_1_55_rev = data_model.Reaction(
#             participants=[
#                 data_model.ReactionParticipant(
#                     specie=data_model.Specie(
#                         id='Lactaldehyde',
#                         structure='InChI=1S/C3H6O2xxxxxx/c1-3(5)2-4/h2-3,5H,1H3'),
#                     coefficient=-1),
#                 data_model.ReactionParticipant(
#                     specie=data_model.Specie(
#                         id='NADPH',
#                         structure=(
#                             'InChI=1S/C21H30N7O17P3'
#                             '/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)'
#                             '40-5-10-13(29)15(31)20(42-10)27-3-1-2-9(4-27)18(23)32'
#                             '/h1,3-4,7-8,10-11,13-16,20-21,29-31H,2,5-6H2,(H2,23,32)(H,36,37)(H,38,39)(H2,22,24,25)(H2,33,34,35)'
#                             '/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1')),
#                     coefficient=-1),
#                 data_model.ReactionParticipant(
#                     specie=data_model.Specie(
#                         id='H+',
#                         structure='InChI=1S/p+1'),
#                     coefficient=-1),
#                 data_model.ReactionParticipant(
#                     specie=data_model.Specie(
#                         id='1,2-Propanediol',
#                         structure='InChI=1S/C3H8O2/c1-3(5)2-4/h3-5H,2H2,1H3'),
#                     coefficient=1),
#                 data_model.ReactionParticipant(
#                     specie=data_model.Specie(
#                         id='NADP+',
#                         structure=(
#                             'InChI=1S/C21H28N7O17P3'
#                             '/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)'
#                             '45-47(36,37)40-5-10-13(29)15(31)20(42-10)27-3-1-2-9(4-27)18(23)32'
#                             '/h1-4,7-8,10-11,13-16,20-21,29-31H,5-6H2,(H7-,22,23,24,25,32,33,34,35,36,37,38,39)'
#                             '/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1')),
#                     coefficient=1),
#             ],
#             cross_references=[
#                 data_model.Resource(namespace='ec-code', id='1.1.1.55', assignment_method=data_model.ResourceAssignmentMethod.manual),
#             ])
#
#     def test_get_compounds_by_structure(self):
#         q = reaction_kinetics.ReactionKineticsQueryGenerator()
#
#         inchi = 'InChI=1S/C6H13O9P/c7-3-2(1-14-16(11,12)13)15-6(10)5(9)4(3)8/h2-10H,1H2,(H2,11,12,13)/t2-,3-,4+,5-,6+/m1/s1'
#
#         compounds = q.get_compounds_by_structure(inchi).all()
#         self.assertEqual(set(c.id for c in compounds), set([24, 1323, 1357, 1404, 1405, 1493, 1517, 1522, 24709]))
#         # 1517 is a match based on the InChI conversion of its SMILES representation
#         # 1524 is not in the sqlite database beacuse it it reachable from any of reaction
#
#         compounds = q.get_compounds_by_structure(inchi, only_formula_and_connectivity=False).all()
#         self.assertEqual([c.id for c in compounds], [24])
#
#         # select=sabio_rk.Compound.id
#         compounds = q.get_compounds_by_structure(inchi, only_formula_and_connectivity=False, select=sabio_rk.Compound.id).all()
#         self.assertEqual([c[0] for c in compounds], [24])
#
    # def test_get_kinetic_laws_by_compound(self):
    #     q = reaction_kinetics.ReactionKineticsQueryGenerator()
    #
    #     d_Lactaldehyde = 'InChI=1S/C3H6O2/c1-3(5)2-4/h2-3,5H,1H3/t3-/m1/s1'
    #
    #     # whole structure, reactant, select id
    #     q_law = q.get_kinetic_laws_by_compound(
    #         d_Lactaldehyde, only_formula_and_connectivity=False,
    #         role='reactant', select=sabio_rk.KineticLaw.id) \
    #         .order_by(sabio_rk.KineticLaw.id)
    #     self.assertEqual([l[0] for l in q_law.all()], [22870, 22877, 22882, 38452, 44597, 50464, 50469, 50470])
    #
    #     # whole structure, reactant, select object
    #     q_law = q.get_kinetic_laws_by_compound(
    #         d_Lactaldehyde, only_formula_and_connectivity=False,
    #         role='reactant', select=sabio_rk.KineticLaw) \
    #         .order_by(sabio_rk.KineticLaw.id)
    #     self.assertEqual([l.id for l in q_law.all()], [22870, 22877, 22882, 38452, 44597, 50464, 50469, 50470])
    #
    #     rxn_ids = set()
    #     for law in q_law.all():
    #         rxn_ids.add(next(xr.id for xr in law.cross_references if xr.namespace == 'sabiork.reaction'))
    #     self.assertEqual(rxn_ids, set(['548', '10424']))
    #
    #     # formula and connectivity, reactant, select object
    #     q_law = q.get_kinetic_laws_by_compound(
    #         molecule_util.InchiMolecule(d_Lactaldehyde).get_formula_and_connectivity(), only_formula_and_connectivity=True,
    #         role='reactant', select=sabio_rk.KineticLaw)
    #     law_ids = set([l.id for l in q_law.all()])
    #     self.assertGreater(len(law_ids.difference(set([22870, 22877, 22882, 38452, 44597, 50464, 50469, 50470]))), 0)
    #
    #     rxn_ids = set()
    #     for law in q_law.all():
    #         rxn_ids.add(next(xr.id for xr in law.cross_references if xr.namespace == 'sabiork.reaction'))
    #     self.assertGreater(len(rxn_ids.difference(set(['548', '10424']))), 0)
    #
    #     # whole structure, product, select object
    #     q_law = q.get_kinetic_laws_by_compound(
    #         d_Lactaldehyde, only_formula_and_connectivity=False,
    #         role='product', select=sabio_rk.KineticLaw)
    #     self.assertEqual([l.id for l in q_law.all()], [44603])
    #
    #     rxn_ids = set()
    #     for law in q_law.all():
    #         rxn_ids.add(next(xr.id for xr in law.cross_references if xr.namespace == 'sabiork.reaction'))
    #     self.assertEqual(rxn_ids, set(['10424']))
#
#     def test_get_kinetic_laws_by_participants(self):
#         q = reaction_kinetics.ReactionKineticsQueryGenerator()
#
#         participants = self.reaction_1_1_1_55.participants
#
#         # only_formula_and_connectivity=True, include_water_hydrogen=False
#         laws = q.get_kinetic_laws_by_participants(participants) \
#             .order_by(sabio_rk.KineticLaw.id)
#         self.assertEqual([l.id for l in laws], [22870, 22871, 22874, 22876, 22877,
#                                                 22880, 22882, 28503, 38452, 38455, 44597, 46425, 46426, 46428])
#
#         rxn_ids = set()
#         for law in laws:
#             rxn_ids.add(next(xr.id for xr in law.cross_references if xr.namespace == 'sabiork.reaction'))
#         self.assertEqual(rxn_ids, set(['554', '2227', '10424', '10434']))
#
#         # only_formula_and_connectivity=False
#         laws = q.get_kinetic_laws_by_participants(participants, only_formula_and_connectivity=False) \
#             .order_by(sabio_rk.KineticLaw.id)
#         self.assertEqual([l.id for l in laws], [22870, 22871, 22877, 22882, 38452, 38455, 44597, 46425, 46426, 46428])
#
#         rxn_ids = set()
#         for law in laws:
#             rxn_ids.add(next(xr.id for xr in law.cross_references if xr.namespace == 'sabiork.reaction'))
#         self.assertEqual(rxn_ids, set(['554', '2227', '10424']))
#
#         # include_water_hydrogen=True
#         laws = q.get_kinetic_laws_by_participants(participants, include_water_hydrogen=True) \
#             .order_by(sabio_rk.KineticLaw.id)
#         self.assertEqual([l.id for l in laws], [22870, 22871, 22874, 22876, 22877,
#                                                 22880, 22882, 28503, 38452, 38455, 44597, 46425, 46426, 46428])
#
#         rxn_ids = set()
#         for law in laws:
#             rxn_ids.add(next(xr.id for xr in law.cross_references if xr.namespace == 'sabiork.reaction'))
#         self.assertEqual(rxn_ids, set(['554', '2227', '10424', '10434']))
#
#         # select=sabio_rk.KineticLaw.id
#         laws = q.get_kinetic_laws_by_participants(participants, select=sabio_rk.KineticLaw.id) \
#             .order_by(sabio_rk.KineticLaw.id)
#         self.assertEqual([l[0] for l in laws], [22870, 22871, 22874, 22876, 22877,
#                                                 22880, 22882, 28503, 38452, 38455, 44597, 46425, 46426, 46428])
#
#
#     def test_get_kinetic_laws_by_reaction(self):
#         q = reaction_kinetics.ReactionKineticsQueryGenerator()
#
#         laws = q.get_kinetic_laws_by_reaction(self.reaction_1_1_1_55) \
#             .order_by(sabio_rk.KineticLaw.id) \
#             .all()
#         self.assertEqual([l.id for l in laws], [22870, 22871, 22874, 22876, 22877,
#                                                 22880, 22882, 28503, 38452, 38455, 44597, 46425, 46426, 46428])
#
#         # select=sabio_rk.KineticLaw.id
#         laws = q.get_kinetic_laws_by_reaction(self.reaction_1_1_1_55, select=sabio_rk.KineticLaw.id) \
#             .order_by(sabio_rk.KineticLaw.id) \
#             .all()
#         self.assertEqual([l[0] for l in laws], [22870, 22871, 22874, 22876, 22877,
#                                                 22880, 22882, 28503, 38452, 38455, 44597, 46425, 46426, 46428])
#
#         # get by EC
#         laws = q.get_kinetic_laws_by_reaction(self.reaction_1_1_1_55_rev) \
#             .order_by(sabio_rk.KineticLaw.id) \
#             .all()
#         self.assertEqual([l.id for l in laws], [46425, 46426, 46427, 46428])
#
#
#     def print_observed_values(self, vals):
#         print('\n')
#         print('{:<9}  {:<16}  {:<20}  {:<20}  {:<5}  {:<22}  {:<9}  {:<11}  {:<3}  {:<20}  {:<17}'.format(
#             'Parameter', 'Species', 'Compartment', 'Value', 'Units', 'Taxon', 'Variation',
#             'Temperature', 'pH', 'SABIO-RK kinetic law', 'SABIO-RK reaction'))
#         print('{:<9}  {:<16}  {:<20}  {:<20}  {:<5}  {:<22}  {:<9}  {:<11}  {:<3}  {:<20}  {:<17}'.format(
#             '=' * 9, '=' * 16, '=' * 20, '=' * 20, '=' * 5, '=' * 22, '=' * 9, '=' * 11, '=' * 3, '=' * 20, '=' * 17))
#         for v in vals:
#             print('{:<9}  {:<16}  {:<20}  {:>20}  {:<5}  {:<22}  {:<9}  {:>11}  {:>0.1f}  {:>20}  {:>17}'.format(
#                 v.observable.property,
#                 v.observable.specie.name if v.observable.specie else '',
#                 v.observable.compartment.id if v.observable.compartment else '',
#                 v.value, v.units or '',
#                 taxonomy_util.Taxon(ncbi_id=v.observation.genetics.taxon).name, v.observation.genetics.variation or '',
#                 v.observation.environment.temperature or '', v.observation.environment.ph or '',
#                 next(xr.id for xr in v.observable.interaction.cross_references if xr.namespace == 'sabiork.kineticrecord'),
#                 next(xr.id for xr in v.observable.interaction.cross_references if xr.namespace == 'sabiork.reaction'),
#             ))
