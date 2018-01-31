# -*- coding: utf-8 -*-

""" Tests of sabio_rk

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-04-26
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model
from kinetic_datanator.data_source import sabio_rk
from kinetic_datanator.data_query import reaction_kinetics
from kinetic_datanator.util import taxonomy_util, molecule_util
from kinetic_datanator.app import models, flask_common_schema
import unittest
import random
import tempfile
import shutil

# #TODO: Make Print and Filtering functions
# @unittest.skip('skip')
# class TestwPartialCommonSchemaReactionKineticsQueryGenerator(unittest.TestCase):
#     """
#     Tests for 10000 entry limited Common Schema on Karr Lab Server
#
#     Used for development purposes
#     """
#
#     def setUp(self):
#         self.reaction = data_model.Reaction(
#             participants = [
#                 data_model.ReactionParticipant(
#                     specie = data_model.Specie(
#                         id = 'Dihydrofolate',
#                         structure = 'InChI=1S/C19H21N7O6/c20-19-25-15-14(17(30)26-19)23-11(8-22-15)'+\
#                         '7-21-10-3-1-9(2-4-10)16(29)24-12(18(31)32)5-6-13(27)28/h1-4,12,21H,5-8H2,'+\
#                         '(H,24,29)(H,27,28)(H,31,32)(H4,20,22,25,26,30)'),
#                     coefficient = -1),
#                 data_model.ReactionParticipant(
#                     specie = data_model.Specie(
#                         id = 'NADPH',
#                         structure = 'InChI=1S/C21H30N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16'+\
#                         '(44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)40-5-10-13(29)15(31)'+\
#                         '20(42-10)27-3-1-2-9(4-27)18(23)32/h1,3-4,7-8,10-11,13-16,20-21,29-31H,2,5-6H2,'+\
#                         '(H2,23,32)(H,36,37)(H,38,39)(H2,22,24,25)(H2,33,34,35)'),
#                     coefficient = -1),
#                 data_model.ReactionParticipant(
#                     specie = data_model.Specie(
#                         id = 'H+',
#                         structure = 'InChI=1S/H'),
#                     coefficient = -1),
#                 data_model.ReactionParticipant(
#                     specie = data_model.Specie(
#                         id = 'NADP+',
#                         structure = 'InChI=1S/C21H28N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16('+\
#                         '44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)40-5-10-13(29)15(31)20'+\
#                         '(42-10)27-3-1-2-9(4-27)18(23)32/h1-4,7-8,10-11,13-16,20-21,29-31H,5-6H2,(H7-,22'+\
#                         ',23,24,25,32,33,34,35,36,37,38,39)/p+1'),
#                     coefficient = 1),
#                 data_model.ReactionParticipant(
#                     specie = data_model.Specie(
#                         id = '5,6,7,8-Tetrahydrofolate',
#                         structure = 'InChI=1S/C19H23N7O6/c20-19-25-15-14(17(30)26-19)23-11(8-22-15)7-21-'+\
#                         '10-3-1-9(2-4-10)16(29)24-12(18(31)32)5-6-13(27)28/h1-4,11-12,21,23H,5-8H2,(H,24,29'+\
#                         ')(H,27,28)(H,31,32)(H4,20,22,25,26,30)'),
#                     coefficient = 1)
#         ])
#
#         self.reaction_w_resource = data_model.Reaction(
#             participants = [
#                 data_model.ReactionParticipant(
#                     specie = data_model.Specie(
#                         id = 'Dihydrofolate',
#                         structure = 'InChI=1S/C19H21N7O6/c20-19-25-15-14(17(30)26-19)23-11(8-22-15)'+\
#                         '7-21-10-3-1-9(2-4-10)16(29)24-12(18(31)32)5-6-13(27)28/h1-4,12,21H,5-8H2,'+\
#                         '(H,24,29)(H,27,28)(H,31,32)(H4,20,22,25,26,30)'),
#                     coefficient = -1),
#                 data_model.ReactionParticipant(
#                     specie = data_model.Specie(
#                         id = 'NADPH',
#                         structure = 'InChI=1S/C21H30N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16'+\
#                         '(44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)40-5-10-13(29)15(31)'+\
#                         '20(42-10)27-3-1-2-9(4-27)18(23)32/h1,3-4,7-8,10-11,13-16,20-21,29-31H,2,5-6H2,'+\
#                         '(H2,23,32)(H,36,37)(H,38,39)(H2,22,24,25)(H2,33,34,35)'),
#                     coefficient = -1),
#                 data_model.ReactionParticipant(
#                     specie = data_model.Specie(
#                         id = 'H+',
#                         structure = 'InChI=1S/H'),
#                     coefficient = -1),
#                 data_model.ReactionParticipant(
#                     specie = data_model.Specie(
#                         id = 'NADP+',
#                         structure = 'InChI=1S/C21H28N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16('+\
#                         '44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)40-5-10-13(29)15(31)20'+\
#                         '(42-10)27-3-1-2-9(4-27)18(23)32/h1-4,7-8,10-11,13-16,20-21,29-31H,5-6H2,(H7-,22'+\
#                         ',23,24,25,32,33,34,35,36,37,38,39)/p+1'),
#                     coefficient = 1),
#                 data_model.ReactionParticipant(
#                     specie = data_model.Specie(
#                         id = '5,6,7,8-Tetrahydrofolate',
#                         structure = 'InChI=1S/C19H23N7O6/c20-19-25-15-14(17(30)26-19)23-11(8-22-15)7-21-'+\
#                         '10-3-1-9(2-4-10)16(29)24-12(18(31)32)5-6-13(27)28/h1-4,11-12,21,23H,5-8H2,(H,24,29'+\
#                         ')(H,27,28)(H,31,32)(H4,20,22,25,26,30)'),
#                     coefficient = 1)
#         ],
#             cross_references=[
#                 data_model.Resource(namespace='ec-code', id='1.5.1.3', assignment_method=data_model.ResourceAssignmentMethod.manual),
#         ])
#
#     def test_get_kinetic_laws_by_ec_numbers(self):
#         q = reaction_kinetics.ReactionKineticsQueryGenerator()
#         com = common_schema.CommonSchema()
#         ses = com.session
#
#         # single EC, match_levels=4
#         laws_62 = q.get_kinetic_laws_by_ec_numbers(['3.4.21.62'], match_levels=4).all()
#         compare = ses.query(common_schema.KineticLaw).filter_by(enzyme_id = laws_62[0].enzyme_id).all()
#         ids_62 = set([c.kinetic_law_id for c in compare])
#         self.assertEqual(set(l.id for l in laws_62), ids_62)
#
#         laws_73 = q.get_kinetic_laws_by_ec_numbers(['3.4.21.73'], match_levels=4).all()
#         compare = ses.query(common_schema.KineticLaw).filter_by(enzyme_id = laws_73[0].enzyme_id).all()
#         ids_73 = set([c.kinetic_law_id for c in compare])
#         self.assertEqual(set(l.id for l in laws_73), ids_73)
#
#         laws_89 = q.get_kinetic_laws_by_ec_numbers(['3.4.21.89'], match_levels=4).all()
#         compare = ses.query(common_schema.KineticLaw).filter_by(enzyme_id = laws_89[0].enzyme_id).all()
#         ids_89 = set([c.kinetic_law_id for c in compare])
#         self.assertEqual(set(l.id for l in laws_89), ids_89)
#
#         #multiple EC, match_levels=4
#
#         laws = q.get_kinetic_laws_by_ec_numbers(['3.4.21.62', '3.4.21.73'], match_levels=4).all()
#         self.assertEqual(set([l.id for l in laws]), ids_62 | ids_73)
#
#         laws = q.get_kinetic_laws_by_ec_numbers(['3.4.21.62', '3.4.21.89'], match_levels=4).all()
#         self.assertEqual(set([l.id for l in laws]), ids_62 | ids_89)
#
#         laws = q.get_kinetic_laws_by_ec_numbers(['3.4.21.73', '3.4.21.89'], match_levels=4).all()
#         self.assertEqual(set([l.id for l in laws]), ids_73 | ids_89)
#
#         laws = q.get_kinetic_laws_by_ec_numbers(['3.4.21.62', '3.4.21.73', '3.4.21.89'], match_levels=4).all()
#         self.assertEqual(set([l.id for l in laws]), ids_62 | ids_73 | ids_89)
#
#
#         # match_levels=3
#         laws = q.get_kinetic_laws_by_ec_numbers(['3.4.21.62'], match_levels=3).all()
#         ids = set([l.id for l in laws])
#         ids_contains = ids_62 | ids_73 | ids_89
#         self.assertEqual(ids_contains.difference(ids), set())
#
#     def test_get_kinetic_laws_by_compound(self):
#         q = reaction_kinetics.ReactionKineticsQueryGenerator()
#
#         inchi = 'InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)'
#
#         law = q.get_kinetic_laws_by_compound(inchi, role = 'product')
#         self.assertEqual(set(c.kinetic_law_id for c in law.all() if c.kinetic_law_id <400000 ), set([60256,60252]))
#
#
#         inchi = 'InChI=1S/H5O10P3/c1-11(2,3)9-13(7,8)10-12(4,5)6/h(H,7,8)(H2,1,2,3)(H2,4,5,6)'
#         law = q.get_kinetic_laws_by_compound(inchi)
#         self.assertEqual(set(c.kinetic_law_id for c in law.all() if c.kinetic_law_id <400000 ), set([59398, 59399, 59400,
#             59401, 59406, 59407, 59408, 59409]))
#
#
#     def test_get_compounds_by_structure(self):
#         q = reaction_kinetics.ReactionKineticsQueryGenerator()
#
#         inchi = 'InChI=1S/C5H8O3/c1-3(2)4(6)5(7)8/h3H,1-2H3,(H,7,8)'
#
#         compounds = q.get_compounds_by_structure(inchi).all()
#         self.assertEqual(set(c.id for c in compounds if c.id<400000), set([50387, 51867,]))
#
#     def test_get_kinetic_laws_by_participants(self):
#         q = reaction_kinetics.ReactionKineticsQueryGenerator()
#
#         participants = self.reaction.participants
#
#         law = q.get_kinetic_laws_by_participants(participants, only_formula_and_connectivity = False)
#
#         self.assertEqual([l.kinetic_law_id for l in law if l.kinetic_law_id <400000 ], [59296, 59297, 59298, 59299, 59300, 59301,
#         59308, 59309, 59310, 59311, 59312, 59313, 59314, 59315, 59316, 59317, 59318, 59319, 59320,
#         59321, 59322, 59323, 59324, 59325, 59326, 59327, 59328, 59329, 59330, 59331, 59332, 59333,
#         59334, 59335, 59336, 59337, 59338, 59339, 59340, 59341, 59342, 59343, 59344, 59345, 59346,
#         59347, 59348, 59349, 59350, 59351, 59352, 59353, 59354, 59355, 59356, 59357, 59358, 59359])
#
#         for i in range(5):
#             index1 = random.randint(0,len(law)-1)
#             index2 = random.randint(0,len(law)-1)
#             self.assertEqual(law[index1].enzyme, law[index2].enzyme)
#
#     def test_get_kinetic_laws_by_reaction(self):
#         q = reaction_kinetics.ReactionKineticsQueryGenerator()
#
#         laws = q.get_kinetic_laws_by_reaction(self.reaction)
#
#         self.assertEqual([l.kinetic_law_id for l in laws if l.kinetic_law_id <400000], sorted([59296, 59297, 59298, 59299, 59300, 59301,
#         59308, 59309, 59310, 59311, 59312, 59313, 59314, 59315, 59316, 59317, 59318, 59319, 59320,
#         59321, 59322, 59323, 59324, 59325, 59326, 59327, 59328, 59329, 59330, 59331, 59332, 59333,
#         59334, 59335, 59336, 59337, 59338, 59339, 59340, 59341, 59342, 59343, 59344, 59345, 59346,
#         59347, 59348, 59349, 59350, 59351, 59352, 59353, 59354, 59355, 59356, 59357, 59358, 59359], key=int))
#
#         laws = q.get_kinetic_laws_by_reaction(self.reaction_w_resource)
#
#         self.assertEqual([l.kinetic_law_id for l in laws if l.kinetic_law_id <400000], sorted([59317, 59354, 59323, 59356, 59322, 59311,
#         59300, 59332, 59318, 59333, 59351, 59299, 59346, 59320, 59337, 59315, 59331, 59328, 59308,
#         59301, 59309, 59324, 59347, 59352, 59326, 59345, 59327, 59342, 59359, 59297, 59296, 59310,
#         59329, 59321, 59298, 59312, 59314, 59339, 59353, 59340, 59344, 59313, 59338, 59355, 59334,
#         59350, 59319, 59358, 59343, 59357, 59330, 59336, 59325, 59349, 59341, 59335, 59316, 59348], key=int))
#
#     def test_get_observed_values(self):
#         q = reaction_kinetics.ReactionKineticsQueryGenerator()
#         vals = q.get_observed_values(self.reaction)
#
#         for val in vals:
#             sabiork_id = next(xr.id for xr in val.observable.interaction.cross_references if xr.namespace == 'sabiork.reaction')
#             common_schema_kinetic_id = next(xr.id for xr in val.observable.interaction.cross_references if xr.namespace == 'common_schema.kinetic_law_id')
#             if sabiork_id == '100' and common_schema_kinetic_id == '2205':
#                 if val.observable.property == 'Km_A':
#                     self.assertEqual(val.observable.specie.name,'NADPH')
#                     self.assertEqual(val.value, 7.3e-06)
#                     self.assertEqual(val.error, 2.0e-07)
#                     self.assertEqual(val.units, 'M')
#                     break
#
#
#         # self.print_observed_values(vals)
#
#     @unittest.skip('implement me')
#     def test_filter_observed_values(self):
#         q = reaction_kinetics.ReactionKineticsQueryGenerator()
#         vals = q.get_observed_values(self.reaction)
#         filter_result = q.filter_observed_values(reaction, vals)
#
#         for ov in filter_result.observed_values:
#             self.assertIn(ov.observable.property, ['k_cat', 'k_i', 'k_m', 'v_max'])
#             if ov.observable.specie:
#                 self.assertIn(ov.observable.specie.to_inchi(only_formula_and_connectivity=True), [
#                     'C3H6O2/c1-3(5)2-4',
#                     ('C21H30N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47'
#                      '(36,37)40-5-10-13(29)15(31)20(42-10)27-3-1-2-9(4-27)18(23)32'),
#                 ])
#             self.assertEqual(ov.observation.genetics.variation, '')
#
#         # self.print_observed_values(filter_result.observed_values)
#
#     @unittest.skip('implement me')
#     def test_get_consensus(self):
#         mol = data_model.Specie(structure=inchi)
#
#         ec_numbers = ['1.1.1.52', '1.1.1.55']
#         rxn = data_model.Reaction(cross_references=[
#             data_model.Resource(namespace='ec-code', id=ec_number),
#         ])

class TestFlaskCommonSchemaReactionKineticsQueryGenerator(unittest.TestCase):
    """
    Tests for 10000 entry limited Common Schema on Karr Lab Server

    Used for development purposes
    """

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.flk = flask_common_schema.FlaskCommonSchema(cache_dirname=self.cache_dirname)

        self.q = reaction_kinetics.FlaskReactionKineticsQueryGenerator()


    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)


    def test_get_kinetic_laws_by_ec_numbers(self):


        # single EC, match_levels=4
        laws_62 = self.q.get_kinetic_laws_by_ec_numbers(['3.4.21.62'], match_levels=4).all()
        compare = self.flk.session.query(models.KineticLaw).filter_by(enzyme_id = laws_62[0].enzyme_id).all()
        ids_62 = set([c.kinetic_law_id for c in compare])
        self.assertEqual(set(l.id for l in laws_62), ids_62)

        laws_73 = self.q.get_kinetic_laws_by_ec_numbers(['3.4.21.73'], match_levels=4).all()
        compare = self.flk.session.query(models.KineticLaw).filter_by(enzyme_id = laws_73[0].enzyme_id).all()
        ids_73 = set([c.kinetic_law_id for c in compare])
        self.assertEqual(set(l.id for l in laws_73), ids_73)

        #multiple EC, match_levels=4

        laws = self.q.get_kinetic_laws_by_ec_numbers(['3.4.21.62', '3.4.21.73'], match_levels=4).all()
        self.assertEqual(set([l.id for l in laws]), ids_62 | ids_73)


    def test_get_kinetic_laws_by_compound(self):

        compound = self.flk.session.query(models.Compound).filter_by(compound_name = '2-Hydroxyisocaproate').first()

        law = self.q.get_kinetic_laws_by_compound(compound, role = 'reactant')
        self.assertEqual(set(c.kinetic_law_id for c in law.all() if c.kinetic_law_id <400000 ), set([20331, 20307]))

        compound = self.flk.session.query(models.Compound).filter_by(compound_name = '4-Hydroxyphenylethyl alcohol').first()

        law = self.q.get_kinetic_laws_by_compound(compound, role = 'product')
        self.assertEqual(set(c.kinetic_law_id for c in law.all() if c.kinetic_law_id <400000 ), set([23314]))


    def test_get_compounds_by_structure(self):
        struct = self.flk.session.query(models.Structure).all()

        ans = self.q.get_compounds_by_structure(struct[3414]._value_inchi).all()
        self.assertEqual(ans, struct[3414].compound)


    def test_get_reaction_by_compound(self):
        compound = self.flk.session.query(models.Compound).filter_by(compound_name = '2-Hydroxyisocaproate').first()

        rxn_list = self.q.get_reaction_by_compound(compound)

        self.assertEqual(len(rxn_list[0].participants), 4)
        self.assertEqual(rxn_list[0].participants[0].specie.id, '2-Hydroxyisocaproate')

    def test_get_kinetic_laws_by_reaction(self):
        pass

    def test_get_observed_values(self):
        pass
