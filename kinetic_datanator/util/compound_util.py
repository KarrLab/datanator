""" Utilities for dealing with compounds

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan <jonrkarr@gmail.com>
:Date: 2017-04-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

import openbabel


class Compound(object):
    """ Represents a compound

    Attributes:
        structure (:obj:`str`): structure in InChI format77
        input_structure (:obj:`str`): structure in input format
        input_structure_format (:obj:`str`): format of the input structure
    """

    def __init__(self, structure):
        """
        Args:
            structure (:obj:`str`): structure in InChI, MOL, or canonical SMILES format

        Raises:
            :obj:`ValueError`: if the structure is not valid
        """

        mol = openbabel.OBMol()
        obConversion = openbabel.OBConversion()

        # read structure
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

    def to_openbabel(self):
        """ Create an Open Babel molecule for the compound

        Returns:
            :obj:`openbabel.OBMol`: Open Babel molecule
        """
        mol = openbabel.OBMol()
        obConversion = openbabel.OBConversion()
        obConversion.SetInFormat('inchi')
        obConversion.ReadString(mol, self.structure)
        return mol

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
