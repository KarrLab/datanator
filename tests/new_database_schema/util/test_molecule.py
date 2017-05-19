from kinetic_datanator.new_database_schema.util import molecule_util
import unittest
import numpy



class TestMolecule(unittest.TestCase):
	adp = {'smiles': 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'}
	h = {'smiles': '[H+]'}

	h2o = {
		'name':'water',
		'id':{
			"KEGG": "C00001",
			"PUBCHEM": "3303",
			"SABIO-RK": "34",
			"CHEBI": "15377"
		},
		'inchi': 'InChI=1S/H2O/h1H2',
		'smiles': 'O',
		'mol': (
			'\n'
			' OpenBabel04121720292D\n'
			'\n'
			'  1  0  0  0  0  0  0  0  0  0999 V2000\n'
			'    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n'
			'M  END'
		)

	}

	atp = {
		'name':'ATP',
		'id':{
			"KEGG": "C00002",
			"PUBCHEM": "3304",
			"SABIO-RK": "34",
			"CHEBI": "15422"
		},
		'smiles': 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1',
	}


	def test_init(self):
		c = molecule_util.Molecule(name = self.atp['name'], id=self.atp['id'], structure=self.atp['smiles'])

		self.assertEqual(c.id, self.atp['id'])
		self.assertEqual(c.name, self.atp['name'])
		self.assertEqual(c.structure, self.atp['smiles'])
		self.assertEqual(c.get_format(), 'can')
		

		#c = molecule_util.Molecule(id='h2o', name='water', structure=self.h2o['inchi'])
		c = molecule_util.Molecule(name = self.h2o['name'], id=self.h2o['id'], structure=self.h2o['inchi'])
		self.assertEqual(c.id, self.h2o['id'])
		self.assertEqual(c.name, self.h2o['name'])
		self.assertEqual(c.structure, self.h2o['inchi'])
		self.assertEqual(c.get_format(), 'inchi')

		self.assertEqual(c.to_openbabel().GetFormula(), 'H2O')
		self.assertEqual(c.to_pybel().formula, 'H2O')
		self.assertEqual(c.to_inchi(), self.h2o['inchi'])
		self.assertEqual(c.to_mol().split('\n')[2:], self.h2o['mol'].split('\n')[2:])
		self.assertEqual(c.to_smiles(), self.h2o['smiles'])

		self.assertEqual(molecule_util.Molecule(structure=self.h2o['inchi'][6:]).get_format(), None)

		c = molecule_util.Molecule(structure=self.h2o['inchi'] + '\n')
		self.assertEqual(c.structure, self.h2o['inchi'] + '\n')

		c = molecule_util.Molecule(structure=self.h2o['smiles'])
		self.assertEqual(c.structure, self.h2o['smiles'])

		c = molecule_util.Molecule(structure=self.h2o['mol'])
		self.assertEqual(c.structure, self.h2o['mol'])

	def test_get_id(self):
		ATP = molecule_util.Molecule(name = self.atp['name'], id=self.atp['id'], structure=self.atp['smiles'])
		self.assertEqual(ATP.get_id("KEGG"), self.atp['id']["KEGG"])
		self.assertEqual(ATP.get_id("PUBCHEM"), self.atp['id']["PUBCHEM"])
		self.assertEqual(ATP.get_id("fake_database"), "")

	def test_get_fingerprint_types(self):
		self.assertIsInstance(molecule_util.Molecule.get_fingerprint_types(), list)

	def test_get_fingerprint(self):
		c = molecule_util.Molecule(structure=self.h2o['inchi'])
		self.assertIsInstance(c.get_fingerprint('fp2'), pybel.Fingerprint)

	def test_get_fingerprint(self):
		adp = molecule_util.Molecule(structure=self.adp['smiles'])
		atp = molecule_util.Molecule(structure=self.atp['smiles'])

		self.assertEqual(adp.get_similarity(adp), 1.)

		numpy.testing.assert_almost_equal(adp.get_similarity(atp), 0.955, decimal=3)
		numpy.testing.assert_almost_equal(atp.get_similarity(adp), 0.955, decimal=3)
		
		