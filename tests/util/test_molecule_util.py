""" Tests of the molecule utilities

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

from datanator.util import molecule_util
from wc_utils.util.types import assert_value_equal
import numpy
import pybel
import unittest


class TestMolecule(unittest.TestCase):
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

    def test_init(self):
        c = molecule_util.Molecule(id='h2o', name='water', structure=self.h2o['inchi'])
        self.assertEqual(c.id, 'h2o')
        self.assertEqual(c.name, 'water')
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


class TestInchiMolecule(unittest.TestCase):

    def test(self):
        inchi = 'InChI=1S/C3H4O3/c1-2(4)3(5)6/h1H3,(H,5,6)'
        layers = molecule_util.InchiMolecule(inchi)
        self.assertEqual(layers.__dict__, {
            'formula': 'C3H4O3',
            'connections': '1-2(4)3(5)6',
            'hydrogens': '1H3,(H,5,6)',
            'protons': '',
            'charge': '',
            'double_bonds': '',
            'stereochemistry': '',
            'stereochemistry_parity': '',
            'stereochemistry_type': '',
            'isotopes': '',
            'fixed_hydrogens': '',
            'reconnected_metals': '',
        })
        self.assertEqual(str(layers), inchi)

        inchi = 'InChI=1S/C3H4O3/c1-2(4)3(5)6'
        layers = molecule_util.InchiMolecule(inchi)
        self.assertEqual(layers.__dict__, {
            'formula': 'C3H4O3',
            'connections': '1-2(4)3(5)6',
            'hydrogens': '',
            'protons': '',
            'charge': '',
            'double_bonds': '',
            'stereochemistry': '',
            'stereochemistry_parity': '',
            'stereochemistry_type': '',
            'isotopes': '',
            'fixed_hydrogens': '',
            'reconnected_metals': '',
        })
        self.assertEqual(str(layers), inchi)

        inchi = 'InChI=1S/Ni/q+2'
        self.assertEqual(str(molecule_util.InchiMolecule(inchi)), inchi)

        inchi = 'InChI=1S/BrH/h1H/p-1'
        self.assertEqual(str(molecule_util.InchiMolecule(inchi)), inchi)

        inchi = 'InChI=1S/p+1'
        self.assertEqual(str(molecule_util.InchiMolecule(inchi)), inchi)

        inchi = 'InChI=1S/4O.V/q;3*-1;'
        self.assertEqual(str(molecule_util.InchiMolecule(inchi)), inchi)

        inchi_1 = 'InChI=1/C6H12O6/c7-1-3-4(9)5(10)6(11,2-8)12-3/h3-5,7-11H,1-2H2/t3-,4-,5+,6u/m0/s1'
        inchi_2 = 'InChI=1S/C6H12O6/c7-1-3-4(9)5(10)6(11,2-8)12-3/h3-5,7-11H,1-2H2/t3-,4-,5+,6u/m0/s1'
        self.assertEqual(str(molecule_util.InchiMolecule(inchi_1)), inchi_2)

    def test_remove_layer(self):
        a = molecule_util.InchiMolecule('InChI=1S/BrH/h1H/p-1')
        a.remove_layer('protons')
        self.assertEqual(str(a), 'InChI=1S/BrH/h1H')

    def test_is_equal(self):
        a = molecule_util.InchiMolecule('InChI=1S/BrH/h1H/p-1')
        c = molecule_util.InchiMolecule('InChI=1S/BrH/h1H')
        self.assertTrue(a.is_equal(a))
        self.assertFalse(a.is_equal(c))

        glc6p = molecule_util.InchiMolecule(
            'InChI=1S/C6H13O9P/c7-3-2(1-14-16(11,12)13)15-6(10)5(9)4(3)8/h2-10H,1H2,(H2,11,12,13)/t2-,3-,4+,5-,6-/m1/s1')
        gal6p = molecule_util.InchiMolecule(
            'InChI=1S/C6H13O9P/c7-3-2(1-14-16(11,12)13)15-6(10)5(9)4(3)8/h2-10H,1H2,(H2,11,12,13)/t2-,3+,4+,5-,6?/m1/s1')
        fru6p = molecule_util.InchiMolecule(
            'InChI=1S/C6H13O9P/c7-1-3(8)5(10)6(11)4(9)2-15-16(12,13)14/h4-7,9-11H,1-2H2,(H2,12,13,14)/t4-,5-,6-/m1/s1')

        self.assertTrue(glc6p.is_equal(glc6p, check_stereochemistry=False))
        self.assertTrue(glc6p.is_equal(gal6p, check_stereochemistry=False))
        self.assertFalse(glc6p.is_equal(fru6p, check_stereochemistry=False))

    def test_is_stereoisomer(self):
        glc6p = molecule_util.InchiMolecule(
            'InChI=1S/C6H13O9P/c7-3-2(1-14-16(11,12)13)15-6(10)5(9)4(3)8/h2-10H,1H2,(H2,11,12,13)/t2-,3-,4+,5-,6-/m1/s1')
        gal6p = molecule_util.InchiMolecule(
            'InChI=1S/C6H13O9P/c7-3-2(1-14-16(11,12)13)15-6(10)5(9)4(3)8/h2-10H,1H2,(H2,11,12,13)/t2-,3+,4+,5-,6?/m1/s1')
        fru6p = molecule_util.InchiMolecule(
            'InChI=1S/C6H13O9P/c7-1-3(8)5(10)6(11)4(9)2-15-16(12,13)14/h4-7,9-11H,1-2H2,(H2,12,13,14)/t4-,5-,6-/m1/s1')

        self.assertTrue(glc6p.is_stereoisomer(glc6p))
        self.assertTrue(glc6p.is_stereoisomer(gal6p))
        self.assertFalse(glc6p.is_stereoisomer(fru6p))

    def test_is_tautomer(self):
        a = molecule_util.InchiMolecule('InChI=1S/C5H5N5O/c6-5-9-3-2(4(11)10-5)7-1-8-3/h1H,(H4,6,7,8,9,10,11)/t1')
        b = molecule_util.InchiMolecule('InChI=1S/C5H5N5O/c6-5-9-3-2(4(11)10-5)7-1-8-3/h1H,(H4,6,7,8,9,10,11)/t2')
        c = molecule_util.InchiMolecule('InChI=1S/C5H5N5O/c6/h1H,(H4,6,7,8,9,10,11)/t2')
        self.assertTrue(a.is_tautomer(a))
        self.assertTrue(a.is_tautomer(a))
        self.assertFalse(a.is_tautomer(c))

    def test_is_protonation_isomer(self):
        a = molecule_util.InchiMolecule('InChI=1S/C6H13O9P/c7-3-2(1-14-16(11,12)13)15-6(10)5(9)4(3)8/h2')
        b = molecule_util.InchiMolecule('InChI=1S/C6H13O9P/c7-3-2(1-14-16(11,12)13)15-6(10)5(9)4(3)8/h3')
        c = molecule_util.InchiMolecule('InChI=1S/C6H13O9P/c7/h4')
        self.assertTrue(a.is_protonation_isomer(a))
        self.assertTrue(b.is_protonation_isomer(b))
        self.assertFalse(a.is_protonation_isomer(c))

    def test_get_formula_and_connectivity(self):
        glc6p = molecule_util.InchiMolecule(
            'InChI=1S/C6H13O9P/c7-3-2(1-14-16(11,12)13)15-6(10)5(9)4(3)8/h2-10H,1H2,(H2,11,12,13)/t2-,3-,4+,5-,6-/m1/s1')
        self.assertEqual(glc6p.get_formula_and_connectivity(), 'C6H13O9P/c7-3-2(1-14-16(11,12)13)15-6(10)5(9)4(3)8')

        water = molecule_util.InchiMolecule('InChI=1S/H2O/h1H2')
        self.assertEqual(water.get_formula_and_connectivity(), 'H2O')
