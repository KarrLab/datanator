""" Utilities for dealing with molecules

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan <jonrkarr@gmail.com>
:Date: 2017-04-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

import openbabel
import pybel
#import rdkit.Chem


class Molecule(object):
    """ Represents a molecule

    Attributes:
        name (:obj:`str`): name
        structure (:obj:`str`): structure in InChI format
        input_structure (:obj:`str`): structure in input format
        input_structure_format (:obj:`str`): format of the input structure
    """

    def __init__(self, structure, name=''):
        """
        Args:
            structure (:obj:`str`): structure in InChI, MOL, or canonical SMILES format
            name (:obj:`str`, optional): name

        Raises:
            :obj:`ValueError`: if the structure is not valid
        """

        # identifier
        self.name = name

        # read structure
        mol = openbabel.OBMol()
        obConversion = openbabel.OBConversion()
        if len(structure) >= 7 and structure[0:6] == 'InChI=' and \
                obConversion.SetInFormat('inchi') and \
                obConversion.ReadString(mol, structure):
            format = 'inchi'
        elif obConversion.SetInFormat('can') and \
                obConversion.ReadString(mol, structure):
            format = 'can'
        elif obConversion.SetInFormat('smiles') and \
                obConversion.ReadString(mol, structure):
            format = 'smiles'
        elif obConversion.SetInFormat('mol') and \
                obConversion.ReadString(mol, structure):
            format = 'mol'
        else:
            raise ValueError('Invalid structure: {}'.format(structure))

        # store input structure and its format
        self.input_structure = structure
        self.input_structure_format = format

        # convert structure to InChI format
        obConversion.SetOutFormat('inchi')
        self.structure = obConversion.WriteString(mol).rstrip()

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
        mol = openbabel.OBMol()
        obConversion = openbabel.OBConversion()
        obConversion.SetInFormat('inchi')
        obConversion.ReadString(mol, self.structure)
        return mol

    def to_pybel(self):
        """ Create a pybel molecule for the molecule

        Returns:
            :obj:`pybel.Molecule`: pybel molecule
        """
        return pybel.readstring('inchi', self.structure)

    def to_rdkit(self):
        """ Create an RDKit molecule for the molecule

        Returns:
            :obj:`rdkit.Chem.rdchem.Mol`: rdkit molecule
        """
        return rdkit.Chem.MolFromSmiles(self.to_smiles())

    def to_format(self, format):
        """ Get the structure in a format

        Args:
            :obj:`str`: format such as inchi, mol, smiles

        Returns:
            :obj:`str`: structure in a format
        """
        mol = self.to_openbabel()
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats('inchi', format)
        obConversion.ReadString(mol, self.structure)
        return obConversion.WriteString(mol).rstrip()

    def to_inchi(self):
        """ Get the structure in InChI format

        Returns:
            :obj:`str`: structure in InChi format
        """
        return self.structure

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
