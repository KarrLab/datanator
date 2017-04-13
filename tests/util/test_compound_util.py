""" Tests of the compound utilities

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.util import compound_util
import unittest


class TestCompoundUtil(unittest.TestCase):

    def test_Compound(self):
        inchi = 'InChI=1S/H2O/h1H2'
        smiles = 'O'
        mol = (
            '\n'
            ' OpenBabel04121720292D\n'
            '\n'
            '  1  0  0  0  0  0  0  0  0  0999 V2000\n'
            '    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n'
            'M  END'
        )

        c = compound_util.Compound(inchi)
        self.assertEqual(c.input_structure, inchi)
        self.assertEqual(c.input_structure_format, 'inchi')
        self.assertEqual(c.structure, inchi)

        self.assertEqual(c.to_openbabel().GetFormula(), 'H2O')
        self.assertEqual(c.to_inchi(), inchi)
        self.assertEqual(c.to_mol().split('\n')[2:], mol.split('\n')[2:])
        self.assertEqual(c.to_smiles(), smiles)

        self.assertRaises(ValueError, compound_util.Compound, inchi[6:])

        c = compound_util.Compound(inchi + '\n')
        self.assertEqual(c.input_structure, inchi + '\n')
        self.assertEqual(c.input_structure_format, 'inchi')
        self.assertEqual(c.structure, inchi)

        c = compound_util.Compound(smiles)
        self.assertEqual(c.input_structure, smiles)
        self.assertEqual(c.input_structure_format, 'can')
        self.assertEqual(c.structure, inchi)

        c = compound_util.Compound(mol)
        self.assertEqual(c.input_structure, mol)
        self.assertEqual(c.input_structure_format, 'mol')
        self.assertEqual(c.structure, inchi)
