""" Utilities for dealing with reactions

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan <jonrkarr@gmail.com>
:Date: 2017-04-13
:Copyright: 2017, Karr Lab
:License: MIT
"""

from . import molecule_util
from . import compartment_util
import numpy


class Reaction(object):
    """ Represents a reaction

    Attributes:
        name (:obj:`str`): name
        participants (:obj:`list` of :obj:`ReactionParticipant`): list of participants in the reaction and
            their compartments and coefficients
        reversible (:obj:`bool`): indicates if the reaction is reversible        
    """

    # todo: calculate :math:`{\Delta}G` and use this instead of the `reversible` boolean-valued attribute

    def __init__(self, participants, reversible=False, name=''):
        """
        Args:
            participants (:obj:`list` of :obj:`ReactionParticipant`): list of participants in the reaction and
                their compartments and coefficients
            reversible (:obj:`bool`, optional): indicates if the reaction is reversible
            name (:obj:`str`, optional): name
        """
        self.participants = participants
        self.reversible = reversible
        self.name = name

    def normalize(self):
        """ Normalize the participant list of a reaction by collapsing repeated participants """
        i_part = 0
        while i_part < len(self.participants):
            part = self.participants[i_part]
            if part.coefficient == 0:
                self.participants.pop(i_part)
                continue

            i_other_part = i_part + 1
            while i_other_part < len(self.participants):
                other_part = self.participants[i_other_part]
                if other_part.molecule.structure == part.molecule.structure and \
                        other_part.molecule.name == part.molecule.name and \
                        other_part.compartment.name == part.compartment.name and \
                        numpy.sign(other_part.coefficient) == numpy.sign(part.coefficient):
                    part.coefficient += other_part.coefficient
                    self.participants.pop(i_other_part)
                else:
                    i_other_part += 1

            i_part += 1

    def get_reactants(self):
        """ Get the reactants of the reaction

        Returns:
            :obj:`list` of :obj:`ReactionParticipant`: list of reactants in the reaction and
                their compartments and coefficients
        """
        return filter(lambda p: p.coefficient < 0, self.participants)

    def get_products(self):
        """ Get the products of the reaction

        Returns:
            :obj:`list` of :obj:`ReactionParticipant`: list of products in the reaction and
                their compartments and coefficients
        """
        return filter(lambda p: p.coefficient > 0, self.participants)


class ReactionParticipant(object):
    """ Represents a participant in a reaction

    Attributes:
        molecule (:obj:`molecule_util.Molecule`): molecule
        compartment (:obj:`compartment_util.Compartment`): compartment
        coefficient (:obj:`float`): coefficient of the molecule/compartment in the reaction
    """

    def __init__(self, molecule, compartment='', coefficient=None):
        """
        Args:
            molecule (:obj:`molecule_util.Molecule`): molecule
            compartment (:obj:`compartment_util.Compartment`): compartment
            coefficient (:obj:`float`): coefficient of the molecule/compartment in the reaction
        """
        self.molecule = molecule
        self.compartment = compartment
        self.coefficient = coefficient
