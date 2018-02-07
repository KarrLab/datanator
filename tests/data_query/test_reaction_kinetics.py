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
from kinetic_datanator.core import models, flask_common_schema
import unittest
import random
import tempfile
import shutil

class TestReactionKineticsQueryGenerator(unittest.TestCase):
    """
    Tests for 10000 entry limited Common Schema on Karr Lab Server

    Used for development purposes
    """

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.flk = flask_common_schema.FlaskCommonSchema(cache_dirname=self.cache_dirname)

        self.q = reaction_kinetics.ReactionKineticsQueryGenerator()

        self.reaction = data_model.Reaction(
            participants = [
                data_model.ReactionParticipant(
                    specie = data_model.Specie(
                        id = 'Dihydrofolate',
                        structure = 'InChI=1S/C19H21N7O6/c20-19-25-15-14(17(30)26-19)23-11(8-22-15)'+\
                        '7-21-10-3-1-9(2-4-10)16(29)24-12(18(31)32)5-6-13(27)28/h1-4,12,21H,5-8H2,'+\
                        '(H,24,29)(H,27,28)(H,31,32)(H4,20,22,25,26,30)'),
                    coefficient = -1),
                data_model.ReactionParticipant(
                    specie = data_model.Specie(
                        id = 'NADPH',
                        structure = 'InChI=1S/C21H30N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16'+\
                        '(44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)40-5-10-13(29)15(31)'+\
                        '20(42-10)27-3-1-2-9(4-27)18(23)32/h1,3-4,7-8,10-11,13-16,20-21,29-31H,2,5-6H2,'+\
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
                        structure = 'InChI=1S/C21H28N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16('+\
                        '44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)40-5-10-13(29)15(31)20'+\
                        '(42-10)27-3-1-2-9(4-27)18(23)32/h1-4,7-8,10-11,13-16,20-21,29-31H,5-6H2,(H7-,22'+\
                        ',23,24,25,32,33,34,35,36,37,38,39)/p+1'),
                    coefficient = 1),
                data_model.ReactionParticipant(
                    specie = data_model.Specie(
                        id = '5,6,7,8-Tetrahydrofolate',
                        structure = 'InChI=1S/C19H23N7O6/c20-19-25-15-14(17(30)26-19)23-11(8-22-15)7-21-'+\
                        '10-3-1-9(2-4-10)16(29)24-12(18(31)32)5-6-13(27)28/h1-4,11-12,21,23H,5-8H2,(H,24,29'+\
                        ')(H,27,28)(H,31,32)(H4,20,22,25,26,30)'),
                    coefficient = 1)
        ])


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
        laws = self.q.get_kinetic_laws_by_reaction(self.reaction)

        self.assertEqual([l.kinetic_law_id for l in laws], [20453, 20454, 20455, 20456, 20457,
        20458, 20465, 20466, 20467, 20468, 20469, 20470, 20471, 20472, 20473, 20474, 20475,
        20476, 20477, 20478, 20479, 20480, 20481, 20482, 20483, 20484, 20485, 20486, 20487,
        20488, 20489, 20490, 20491, 20492, 20493, 20494, 20495, 20496, 20497, 20498, 20499,
        20500, 20501, 20502, 20503, 20504, 20505, 20506, 20507, 20508, 20509, 20510, 20511,
        20512, 20513, 20514, 20515, 20516])

    def test_get_observed_values(self):
        vals = self.q.get_observed_values(self.reaction)

        for val in vals:
            sabiork_id = next(xr.id for xr in val.observable.interaction.cross_references if xr.namespace == 'sabiork.reaction')
            common_schema_kinetic_id = next(xr.id for xr in val.observable.interaction.cross_references if xr.namespace == 'common_schema.kinetic_law_id')
            if sabiork_id == '100' and common_schema_kinetic_id == '2205':
                if val.observable.property == 'Km_A':
                    self.assertEqual(val.observable.specie.name,'NADPH')
                    self.assertEqual(val.value, 7.3e-06)
                    self.assertEqual(val.error, 2.0e-07)
                    self.assertEqual(val.units, 'M')
                    break
