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
            self.assertEqual(set(response.keys()),set(['complexes', 'metabolites', 'reactions', 'subunits']) )

    def test_search_metabolite(self):
        pass

    def test_search_subunit(self):
        pass

    def test_search_complex(self):
        pass

    #NOTE: Object Specifict Tests
    def test_metabolite(self):
        with self.client:
            response = self.client.get('/api/v0/metabolite/{0}'.format(72098)).json
            self.assertEqual(set(response.keys()),set(['object','concentrations','reactions']))

    def test_subunit(self):
        with self.client:
            response = self.client.get('/api/v0/subunit/{0}'.format(5)).json
            self.assertEqual(set(response.keys()),set(['object','abundances','interactions', 'complexes']))

    def test_complex(self):
        with self.client:
            response = self.client.get('/api/v0/complex/{0}'.format(69907)).json
            self.assertEqual(set(response.keys()),set(['object','subunits']))

    def test_reaction(self):
        with self.client:
            response = self.client.get('/api/v0/reaction/{0}'.format(70208)).json
            self.assertEqual(set(response.keys()),set(['object','parameters']))

    #NOTE: Data Specific Tests
    def test_metabolite_concentrations(self):
        pass

    def test_protein_abundances(self):
        pass

    def test_protein_interactions(self):
        pass

    def test_reaction_parameters(self):
        pass
