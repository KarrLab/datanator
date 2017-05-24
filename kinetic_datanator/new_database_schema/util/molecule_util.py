""" Utilities for dealing with molecules

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan <jonrkarr@gmail.com>
:Date: 2017-04-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

import openbabel
import pybel
import re
from kinetic_datanator.new_database_schema.util import observable_util


class Molecule(observable_util.Observable):
    """ Represents a molecule

    Attributes:
        name (:obj:`string`): a canonical name for the observable (e.g. ATP)
        id (:obj:'dict' of :obj:`str`: :obj:`str`): dictionary of ids (keys = name of the database (e.g. KEGG), values = identifier)
        structure (:obj:`str`): structure in InChI, MOL, or canonical SMILES format
    """

    # todo: integrate with :obj:`InchiMolecule`

    def __init__(self, name='', id={}, structure=''):
        """
        Args:
            name (:obj:`string`): a canonical name for the observable (e.g. ATP)
            id (:obj:'dict' of :obj:`str`: :obj:`str`): dictionary of ids (keys = name of the database (e.g. KEGG), values = identifier)
            structure (:obj:`str`): structure in InChI, MOL, or canonical SMILES format

        Raises:
            :obj:`ValueError`: if the structure is not valid
        """


        observable_util.Observable.__init__(self, name, id, "metabolite")
        #self.id = id
        #self.name = name
        self.structure = structure

    def get_format(self):
        """ Get the format of the structure

        Returns:
            :obj:`str`: format
        """
        mol = openbabel.OBMol()
        obConversion = openbabel.OBConversion()
        if len(self.structure) >= 7 and self.structure[0:6] == 'InChI=' and \
                obConversion.SetInFormat('inchi') and \
                obConversion.ReadString(mol, self.structure):
            return 'inchi'
        elif obConversion.SetInFormat('can') and \
                obConversion.ReadString(mol, self.structure):
            return 'can'
        elif obConversion.SetInFormat('smiles') and \
                obConversion.ReadString(mol, self.structure):
            return 'smiles'
        elif obConversion.SetInFormat('mol') and \
                obConversion.ReadString(mol, self.structure):
            return 'mol'
        else:
            return None

    def get_fingerprint(self, type='fp2'):
        """ Calculate a fingerprint

        Args:
            type (:obj:`str`, optional): fingerprint type to calculate

        Returns:
            :obj:`pybel.Fingerprint`: fingerprint
        """
        return self.to_pybel().calcfp(type)

    def get_similarity(self, other, fingerprint_type='fp2'):
        """ Calculate the similarity with another molecule

        Args:
            other (:obj:`Molecule`): a second molecule
            fingerprint_type (:obj:`str`, optional): fingerprint type to use to calculate similarity

        Returns:
            :obj:`float`: the similarity with the other molecule
        """
        self_fp = self.get_fingerprint(fingerprint_type)
        other_fp = other.get_fingerprint(fingerprint_type)
        return self_fp | other_fp

    @staticmethod
    def get_fingerprint_types():
        """ Get list of fingerprint types

        Returns:
            :obj:`list` of :obj:`str`: list of fingerprint types
        """
        return pybel.fps

    def to_openbabel(self):
        """ Create an Open Babel molecule for the molecule

        Returns:
            :obj:`openbabel.OBMol`: Open Babel molecule
        """
        format = self.get_format()
        if format is None:
            raise ValueError('Invalid structure: {}'.format(self.structure))

        mol = openbabel.OBMol()
        obConversion = openbabel.OBConversion()
        obConversion.SetInFormat(format)
        obConversion.ReadString(mol, self.structure)
        return mol

    def to_pybel(self):
        """ Create a pybel molecule for the molecule

        Returns:
            :obj:`pybel.Molecule`: pybel molecule
        """
        return pybel.readstring(self.get_format(), self.structure)

    def to_format(self, format):
        """ Get the structure in a format

        Args:
            :obj:`str`: format such as inchi, mol, smiles

        Returns:
            :obj:`str`: structure in a format
        """
        mol = self.to_openbabel()
        obConversion = openbabel.OBConversion()
        obConversion.SetOutFormat(format)
        return obConversion.WriteString(mol).rstrip()

    def to_inchi(self):
        """ Get the structure in InChI format

        Returns:
            :obj:`str`: structure in InChi format
        """
        return self.to_format('inchi')

    def to_mol(self):
        """ Get the structure in MOL format

        Returns:
            :obj:`str`: structure in MOL format
        """
        return self.to_format('mol')

    def to_smiles(self):
        """ Get the structure in SMILES format

        Returns:
            :obj:`str`: structure in SMILES format
        """
        return self.to_format('can')


class InchiMolecule(object):
    """ Represents the InChI-encoded structure of a molecule

    Attributes:
        formula (:obj:`str`): empirical formula layer
        connections (:obj:`str`): atomic conncetions (c) layer
        hydrogens (:obj:`str`): hydrogen (h) layer
        protons (:obj:`str`): proton (p) layer
        charge (:obj:`str`): charge (q) layer
        double_bonds (:obj:`str`): double bounds (b) layer
        stereochemistry (:obj:`str`): stereochemistry (t) layer
        stereochemistry_parity (:obj:`str`): stereochemistry parity (m) layer
        stereochemistry_type (:obj:`str`): stereochemistry type (s) layer
        isotopes (:obj:`str`): isotype (i) layer
        fixed_hydrogens (:obj:`str`): fixed hydrogens (f) layer
        reconnected_metals (:obj:`str`): reconnected metal (r) layer

        LAYERS (:obj:`tuple`): tuple layers and their prefixes
    """
    LAYERS = (
        {'prefix': '',  'name': 'formula', },
        {'prefix': 'c', 'name': 'connections'},
        {'prefix': 'h', 'name': 'hydrogens'},
        {'prefix': 'p', 'name': 'protons'},
        {'prefix': 'q', 'name': 'charge'},
        {'prefix': 'b', 'name': 'double_bonds'},
        {'prefix': 't', 'name': 'stereochemistry'},
        {'prefix': 'm', 'name': 'stereochemistry_parity'},
        {'prefix': 's', 'name': 'stereochemistry_type'},
        {'prefix': 'i', 'name': 'isotopes'},
        {'prefix': 'f', 'name': 'fixed_hydrogens'},
        {'prefix': 'r', 'name': 'reconnected_metals'},
    )

    def __init__(self, structure):
        """
        Args:
            structure (:obj:`str`): InChI-encoded structure of a molecule
        """
        for layer in self.LAYERS:
            start = structure.find('/' + layer['prefix'])
            if start == -1:
                val = ''
            else:
                end = structure.find('/', start + 1)
                if end == -1:
                    val = structure[start + 1 + len(layer['prefix']):]
                else:
                    val = structure[start + 1 + len(layer['prefix']):end]

            setattr(self, layer['name'], val)

    def __str__(self):
        """ Generate an InChI string representation of the molecule

        Returns:
            :obj:`str`: InChI string representation of the molecule
        """
        vals = []
        for layer in self.LAYERS:
            val = getattr(self, layer['name'])
            if val:
                vals.append('/' + layer['prefix'] + val)

        return 'InChI=1S' + ''.join(vals)

    def remove_layer(self, layer):
        """ Remove a layer from a structure

        Args:
            layer (:obj:`str`): name of the layer
        """
        setattr(self, layer, '')

    def is_equal(self, other,
                 check_protonation=True, check_double_bonds=True, check_stereochemistry=True,
                 check_isotopes=True, check_fixed_hydrogens=True, check_reconnected_metals=True):
        """ Determine if two molecules are semantically equal (all of their layers are equal).

        Args:
            other (:obj:`InchiMolecule`): other molecule
            check_protonation (:obj:`bool`, optional): if obj:`True`, check that the protonation states (h, p, q) are equal
            check_double_bonds (:obj:`bool`, optional): if obj:`True`, check that the doubling bonding layers (b) are equal
            check_stereochemistry (:obj:`bool`, optional): if obj:`True`, check that the stereochemistry layers (t, m, s) are equal
            check_isotopes (:obj:`bool`, optional): if obj:`True`, check that the isotopic layers (i) are equal
            check_fixed_hydrogens (:obj:`bool`, optional): if obj:`True`, check that the fixed hydrogen layers (f) are equal
            check_reconnected_metals (:obj:`bool`, optional): if obj:`True`, check that the reconnected metals layers (r) are equal

        Returns:
            :obj:`bool`: :obj:`True` the molecules are semantically equal
        """
        if not isinstance(other, InchiMolecule):
            return False

        if check_protonation:
            if self.formula != other.formula \
                    or self.connections != other.connections \
                    or self.hydrogens != other.hydrogens \
                    or self.protons != other.protons \
                    or self.charge != other.charge:
                return False
        else:
            if re.sub('H[0-9]*', '', self.formula) != re.sub('H[0-9]*', '', other.formula) or \
                    self.connections != other.connections:
                return False

        if check_double_bonds and self.double_bonds != other.double_bonds:
            return False

        if check_stereochemistry and (
                self.stereochemistry != other.stereochemistry or
                self.stereochemistry_parity != other.stereochemistry_parity or
                self.stereochemistry_type != other.stereochemistry_type):
            return False

        if check_isotopes and self.isotopes != other.isotopes:
            return False

        if check_fixed_hydrogens and self.fixed_hydrogens != other.fixed_hydrogens:
            return False

        if check_reconnected_metals and self.reconnected_metals != other.reconnected_metals:
            return False

        return True

    def is_stereoisomer(self, other):
        """ Determine if two molecules are steroisomers

        Args:
            other (:obj:`InchiMolecule`): other molecule

        Returns:
            :obj:`bool`: :obj:`True` if the molecules are stereoisomers
        """
        return isinstance(other, InchiMolecule) \
            and self.formula == other.formula \
            and self.connections == other.connections \
            and self.hydrogens == other.hydrogens \
            and self.protons == other.protons \
            and self.charge == other.charge \
            and self.double_bonds == other.double_bonds

    def is_tautomer(self, other):
        """ Determine if two molecules are tautomers

        Args:
            other (:obj:`InchiMolecule`): other molecule

        Returns:
            :obj:`bool`: :obj:`True` if the molecules are tautomers
        """
        return isinstance(other, InchiMolecule) \
            and self.formula == other.formula \
            and self.connections == other.connections \
            and self.hydrogens == other.hydrogens \
            and self.protons == other.protons \
            and self.charge == other.charge

    def is_protonation_isomer(self, other):
        """ Determine if two molecules are protonation isomers

        Args:
            other (:obj:`InchiMolecule`): other molecule

        Returns:
            :obj:`bool`: :obj:`True` if the molecules are protonation isomers
        """
        return isinstance(other, InchiMolecule) \
            and re.sub('H[0-9]*', '', self.formula) == re.sub('H[0-9]*', '', other.formula) \
            and self.connections == other.connections

    def get_formula_and_connectivity(self):
        """ Get a string representation of the formula and connectivity

        Returns:
            :obj:`str`: string representation of the formula and connectivity
        """
        val = self.formula
        if self.connections:
            val += '/c' + self.connections
        return val
