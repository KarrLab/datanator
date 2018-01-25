""" Tests of common_schemy metabolite concentration queries

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-09-27
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model, common_schema
from kinetic_datanator.data_query import metabolite_concentrations
from kinetic_datanator.flask_datanator import models, flask_common_schema
import tempfile
import shutil
import unittest

class TestMetaboliteConcentrationsQueryGenerator(unittest.TestCase):

    def setUp(self):
        self.proline = data_model.Specie(
            name = 'L-Proline',
            structure = 'InChI=1S/C5H9NO2/c7-5(8)4-2-1-3-6-4/h4,6H,1-3H2,(H,7,8)/t4-/m0/s1'
        )

        self.uridine_tp = data_model.Specie(
            name = 'Uridine triphosphate',
            structure = 'InChI=1S/C9H15N2O15P3/c12-5-1-2-11(9(15)10-5)8-7(14)6(13)4(24-8)3-23'
            '-28(19,20)26-29(21,22)25-27(16,17)18/h1-2,4,6-8,13-14H,3H2,(H,19,20)(H,21,22)(H,10'
            ',12,15)(H2,16,17,18)/t4-,6-,7-,8-/m1/s1'
        )

    def test_filter_observed_values(self):
        q = metabolite_concentrations.MetaboliteConcentrationsQueryGenerator()

        obs = q.get_observed_values(self.proline)
        self.assertEqual(set(c.value for c in obs), set([385.0, 451.0, 361.0, 143.0, 550.0, 531.67]))
        self.assertEqual(set(c.error for c in obs), set([0.0,0.0,0.0,0.0,45.0,11.37]))
        for c in obs:
            self.assertEqual(c.observation.genetics.taxon, 'Escherichia coli')

    def test_get_concentration_by_structure(self):
        q = metabolite_concentrations.MetaboliteConcentrationsQueryGenerator()

        concentrations = q.get_concentration_by_structure(self.proline.structure, only_formula_and_connectivity=False )
        self.assertEqual(set(c.value for c in concentrations.all()), set([385.0, 451.0, 361.0, 143.0, 550.0, 531.67]))

        concentrations = q.get_concentration_by_structure(self.uridine_tp.structure, only_formula_and_connectivity=False )
        self.assertEqual(set(c.value for c in concentrations.all()), set([8290.0, 3990.0, 2370.0, 663.0]))

class TestFlaskMetaboliteConcentrationsQueryGenerator(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        flk = flask_common_schema.FlaskCommonSchema(cache_dirname=self.cache_dirname)

        self.proline = flk.session.query(models.Compound).filter_by(compound_name = 'L-Proline').first()
        self.uridine_tp = flk.session.query(models.Compound).filter_by(compound_name = 'Uridine triphosphate').first()

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)

    def test_filter_observed_values(self):
        q = metabolite_concentrations.FlaskMetaboliteConcentrationsQueryGenerator()

        obs = q.get_observed_values(self.proline)
        self.assertEqual(set(c.value for c in obs), set([385.0, 451.0, 361.0, 143.0, 550.0, 531.67]))
        self.assertEqual(set(c.error for c in obs), set([0.0,0.0,0.0,0.0,45.0,11.37]))
        for c in obs:
            self.assertEqual(c.observation.genetics.taxon, 'Escherichia coli')

    def test_get_concentration_by_structure(self):
        q = metabolite_concentrations.FlaskMetaboliteConcentrationsQueryGenerator()

        concentrations = q.get_concentration_by_structure(self.proline.structure._value_inchi, only_formula_and_connectivity=False )
        self.assertEqual(set(c.value for c in concentrations.all()), set([385.0, 451.0, 361.0, 143.0, 550.0, 531.67]))

        concentrations = q.get_concentration_by_structure(self.uridine_tp.structure._value_inchi, only_formula_and_connectivity=False )
        self.assertEqual(set(c.value for c in concentrations.all()), set([8290.0, 3990.0, 2370.0, 663.0]))
