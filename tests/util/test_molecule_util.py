""" Tests of the molecule utilities

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.util import molecule_util
from kinetic_datanator.util import warning_util
#from rdkit.Chem import rdMolDescriptors
import numpy
import pybel
import unittest

warning_util.set_warnings()


class TestMoleculeUtil(unittest.TestCase):
    adp = {'smiles': 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'}
    atp = {'smiles': 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'}
    h = {'smiles': '[H+]'}
    h2o = {
        'inchi': 'InChI=1S/H2O/h1H2',
        'smiles': 'O',
        'mol': (
            '\n'
            ' OpenBabel04121720292D\n'
            '\n'
            '  1  0  0  0  0  0  0  0  0  0999 V2000\n'
            '    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n'
            'M  END'
        ),
    }

    def test_Molecule(self):
        c = molecule_util.Molecule(self.h2o['inchi'], name='h2o')
        self.assertEqual(c.name, 'h2o')
        self.assertEqual(c.structure, self.h2o['inchi'])
        self.assertEqual(c.get_format(), 'inchi')

        self.assertEqual(c.to_openbabel().GetFormula(), 'H2O')
        self.assertEqual(c.to_pybel().formula, 'H2O')
        #self.assertEqual(rdMolDescriptors.CalcMolFormula(c.to_rdkit()), 'H2O')
        self.assertEqual(c.to_inchi(), self.h2o['inchi'])
        self.assertEqual(c.to_mol().split('\n')[2:], self.h2o['mol'].split('\n')[2:])
        self.assertEqual(c.to_smiles(), self.h2o['smiles'])

        self.assertEqual(molecule_util.Molecule(self.h2o['inchi'][6:]).get_format(), None)

        c = molecule_util.Molecule(self.h2o['inchi'] + '\n')
        self.assertEqual(c.structure, self.h2o['inchi'] + '\n')

        c = molecule_util.Molecule(self.h2o['smiles'])
        self.assertEqual(c.structure, self.h2o['smiles'])

        c = molecule_util.Molecule(self.h2o['mol'])
        self.assertEqual(c.structure, self.h2o['mol'])

    def test_Molecule_get_fingerprint_types(self):
        self.assertIsInstance(molecule_util.Molecule.get_fingerprint_types(), list)

    def test_Molecule_get_fingerprint(self):
        c = molecule_util.Molecule(self.h2o['inchi'])
        self.assertIsInstance(c.get_fingerprint('fp2'), pybel.Fingerprint)

    def test_Molecule_get_fingerprint(self):
        adp = molecule_util.Molecule(self.adp['smiles'])
        atp = molecule_util.Molecule(self.atp['smiles'])

        self.assertEqual(adp.get_similarity(adp), 1.)

        numpy.testing.assert_almost_equal(adp.get_similarity(atp), 0.955, decimal=3)
        numpy.testing.assert_almost_equal(atp.get_similarity(adp), 0.955, decimal=3)
