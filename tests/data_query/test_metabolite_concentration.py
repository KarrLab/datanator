""" Tests of common_schemy metabolite concentration queries

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-09-27
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model
from kinetic_datanator.data_query import metabolite_concentrations
from kinetic_datanator.core import models, flask_common_schema
import tempfile
import shutil
import unittest
import mock

class TestMetaboliteConcentrationsQueryGenerator(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        flk = flask_common_schema.FlaskCommonSchema(cache_dirname=self.cache_dirname)

        self.proline = flk.session.query(models.Compound).filter_by(compound_name = 'L-Proline').first()
        self.uridine_tp = flk.session.query(models.Compound).filter_by(compound_name = 'Uridine triphosphate').first()

        self.q = metabolite_concentrations.MetaboliteConcentrationsQueryGenerator(cache_dirname=self.cache_dirname)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)

    def test_filter_observed_values(self):

        obs = self.q.get_observed_values(self.proline)
        self.assertEqual(set(c.value for c in obs), set([385.0, 451.0, 361.0, 143.0, 550.0, 531.67]))
        self.assertEqual(set(c.error for c in obs), set([0.0,0.0,0.0,0.0,45.0,11.37]))
        for c in obs:
            self.assertEqual(c.observation.genetics.taxon, 'Escherichia coli')

    def test_get_concentration_by_structure(self):

        concentrations = self.q.get_concentration_by_structure(self.proline.structure._value_inchi, only_formula_and_connectivity=False )
        self.assertEqual(set(c.value for c in concentrations.all()), set([385.0, 451.0, 361.0, 143.0, 550.0, 531.67]))

        concentrations = self.q.get_concentration_by_structure(self.uridine_tp.structure._value_inchi, only_formula_and_connectivity=True)
        self.assertEqual(set(c.value for c in concentrations.all()), set([8290.0, 3990.0, 2370.0, 663.0]))
