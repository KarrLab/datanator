from flask_testing import TestCase
import unittest
from datanator import app
from datanator.api.views import *


class TestViews(TestCase):

    def create_app(self):
        return app

    def setUp(self):
        self.general_search = 'adenine'
        self.metabolite_search = 'guanine'
        self.subunit_search = 'egfr'
        self.complex_search = 'MCM complex'
        self.metabolite_id = 72100
        self.subunit_id = 5
        self.complex_id = 69993
        self.reaction_id = 70208

    def tearDown(self):
        pass

    #NOTE: Text Search tests
    def test_search_general(self):
        response = Search().get(self.general_search)
        self.assertEqual(set(response.keys()),set(['complexes', 'metabolites', 'reactions', 'subunits']) )

    def test_search_metabolite(self):
        response = MetaboliteSearch().get(self.metabolite_search)
        self.assertEqual(set(response.keys()),set(['metabolites']) )

    def test_search_subunit(self):
        response = ProteinSubunitSearch().get(self.subunit_search)
        self.assertEqual(set(response.keys()),set(['subunits']) )

    def test_search_complex(self):
        response = ProteinComplexSearch().get(self.complex_search)
        self.assertEqual(set(response.keys()),set(['complexes']) )

    #NOTE: Object Specifict Tests
    def test_metabolite(self):
        response = Metabolite().get(self.metabolite_id)
        self.assertEqual(set(response.keys()),set(['object','concentrations','reactions']))

    def test_subunit(self):
        response = ProteinSubunit().get(self.subunit_id)
        self.assertEqual(set(response.keys()),set(['object','abundances','interactions', 'complexes']))

    def test_complex(self):
        response = ProteinComplex().get(self.complex_id)
        self.assertEqual(set(response.keys()),set(['object','subunits']))

    def test_reaction(self):
        response = Reaction().get(self.reaction_id)
        self.assertEqual(set(response.keys()),set(['object','parameters']))


    #NOTE: Data Specific Tests
    def test_metabolite_concentrations(self):
        response = MetaboliteConcentration().get(self.metabolite_id)
        self.assertEqual(set(response.keys()),set(['concentrations']))

    def test_protein_abundances(self):
        response = ProteinAbundance().get(self.subunit_id)
        self.assertEqual(set(response.keys()),set(['abundances']))

    def test_reaction_parameters(self):
        response = ReactionParameter().get(self.reaction_id)
        self.assertEqual(set(response.keys()),set(['parameters']))

class TestAPIBlueprint(TestCase):

    def create_app(self):
        return app

    def setUp(self):
        self.general_search = 'adenine'
        self.metabolite_search = 'guanine'
        self.subunit_search = 'egfr'
        self.complex_search = 'MCM complex'
        self.metabolite_id = 72100
        self.subunit_id = 5
        self.complex_id = 69993
        self.reaction_id = 70208

    def tearDown(self):
        pass


    #NOTE: Text Search tests
    def test_search_general(self):
        with self.client:
            response = self.client.get('/api/v0/search/{0}'.format(self.general_search)).json
            self.assertEqual(set(response.keys()),set(['complexes', 'metabolites', 'reactions', 'subunits']) )


    def test_search_metabolite(self):
        with self.client:
            response = self.client.get('/api/v0/search/metabolite/{0}'.format(self.metabolite_search)).json
            self.assertEqual(set(response.keys()),set(['metabolites']) )

    def test_search_subunit(self):
            response = self.client.get('/api/v0/search/subunit/{0}'.format(self.subunit_search)).json
            self.assertEqual(set(response.keys()),set(['subunits']) )

    def test_search_complex(self):
            response = self.client.get('/api/v0/search/complex/{0}'.format(self.complex_search)).json
            self.assertEqual(set(response.keys()),set(['complexes']) )

    #NOTE: Object Specifict Tests
    def test_metabolite(self):
        with self.client:
            response = self.client.get('/api/v0/metabolite/{0}'.format(self.metabolite_id)).json
            self.assertEqual(set(response.keys()),set(['object','concentrations','reactions']))

    def test_subunit(self):
        with self.client:
            response = self.client.get('/api/v0/subunit/{0}'.format(self.subunit_id)).json
            self.assertEqual(set(response.keys()),set(['object','abundances','interactions', 'complexes']))

    def test_complex(self):
        with self.client:
            response = self.client.get('/api/v0/complex/{0}'.format(self.complex_id)).json
            self.assertEqual(set(response.keys()),set(['object','subunits']))

    def test_reaction(self):
        with self.client:
            response = self.client.get('/api/v0/reaction/{0}'.format(self.reaction_id)).json
            self.assertEqual(set(response.keys()),set(['object','parameters']))

    #NOTE: Data Specific Tests
    def test_metabolite_concentrations(self):
        with self.client:
            response = self.client.get('/api/v0/concentrations/{0}'.format(self.metabolite_id)).json
            self.assertEqual(set(response.keys()),set(['concentrations']))

    def test_protein_abundances(self):
        with self.client:
            response = self.client.get('/api/v0/abundances/{0}'.format(self.subunit_id)).json
            self.assertEqual(set(response.keys()),set(['abundances']))

    def test_reaction_parameters(self):
        with self.client:
            response = self.client.get('/api/v0/parameters/{0}'.format(self.reaction_id)).json
            self.assertEqual(set(response.keys()),set(['parameters']))
