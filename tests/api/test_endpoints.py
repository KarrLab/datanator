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


    def test_search_general(self):
        with self.client:
            response = self.client.get('/api/v0/search/{0}'.format(self.general_search)).json


    def test_metabolite(self):
        with self.client:
            response = self.client.get('/api/v0/metabolite/{0}'.format(61696)).json
            self.assertEqual(response[0]['id'], 61696)
