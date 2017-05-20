from kinetic_datanator.new_database_schema.util import observable_util
import unittest


class TestObservable(unittest.TestCase):

	#ATP


	atp = {
		'name':'ATP',
		'id':{
			"KEGG": "C00002",
			"PUBCHEM": "3304",
			"SABIO-RK": "34",
			"CHEBI": "15422"
		},
		'obs_type': 'metabolite',
	}

	def test_init(self):
		ATP = observable_util.Observable(name=self.atp['name'], id=self.atp['id'], obs_type=self.atp['obs_type'])
		self.assertEqual(ATP.name, self.atp['name'])
		self.assertEqual(ATP.obs_type, self.atp['obs_type'])
		self.assertEqual(ATP.id["KEGG"], self.atp['id']["KEGG"])
		self.assertEqual(ATP.id["PUBCHEM"], self.atp['id']["PUBCHEM"])
		self.assertEqual(ATP.id["SABIO-RK"], self.atp['id']["SABIO-RK"])
		self.assertEqual(ATP.id["CHEBI"], self.atp['id']["CHEBI"])

		with self.assertRaises(ValueError):
			apple = observable_util.Observable(obs_type="fruit")

	def test_get_id(self):
		ATP = observable_util.Observable(name=self.atp['name'], id=self.atp['id'], obs_type=self.atp['obs_type'])
		self.assertEqual(ATP.get_id("KEGG"), self.atp['id']["KEGG"])
		self.assertEqual(ATP.get_id("PUBCHEM"), self.atp['id']["PUBCHEM"])
		self.assertEqual(ATP.get_id("fake_database"), "")

		
