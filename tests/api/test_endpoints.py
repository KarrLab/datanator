from flask_testing import TestCase
import unittest
from kinetic_datanator import app
from kinetic_datanator.api.views import *

class TestAPIBlueprint(TestCase):

    def create_app(self):
        return app

    def setUp(self):
        self.general_search = 'adenine'
        self.metabolite_search = 'adenine'
        self.subunit_search = 'egfr'
        self.complex_search = 'MCM complex'

    def tearDown(self):
        pass


    #NOTE: Text Search tests
    def test_search_general(self):
        with self.client:
            response = self.client.get('/api/v0/search/{0}'.format(self.general_search)).json
            print(response.keys())

    def test_search_metabolite(self):
        pass

    def test_search_subunit(self):
        pass

    def test_search_complex(self):
        pass

    #NOTE: Object Specifict Tests
    def test_metabolite(self):
        with self.client:
            response = self.client.get('/api/v0/metabolite/{0}'.format(61696)).json
            self.assertEqual(response[0]['id'], 61696)

    def test_subunit(self):
        pass

    def test_complex(self):
        pass

    def test_reaction(self):
        pass

    #NOTE: Data Specific Tests
    def test_metabolite_concentrations(self):
        pass

    def test_protein_abundances(self):
        pass

    def test_protein_interactions(self):
        pass

    def test_reaction_parameters(self):
        pass
