""" Utilities for dealing with reactions

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan <jonrkarr@gmail.com>
:Date: 2017-04-13
:Copyright: 2017, Karr Lab
:License: MIT
"""

from . import compartment_util
from . import molecule_util
from kinetic_datanator import data_structs
from six.moves import zip_longest
import itertools
import numpy
import re
import requests


class Reaction(object):
    """ Represents a reaction

    Attributes:
        id (:obj:`str`): name
        name (:obj:`str`): name
        participants (:obj:`list` of :obj:`ReactionParticipant`): list of participants in the reaction and
            their compartments and coefficients
        reversible (:obj:`bool`): indicates if the reaction is reversible        
        cross_references (:obj:`list` of :obj:`CrossReference`): list of cross references
    """

    # todo: calculate :math:`{\Delta}G` and use this instead of the `reversible` boolean-valued attribute

    def __init__(self, id='', name='', participants=None, reversible=False, cross_references=None):
        """
        Args:
            id (:obj:`str`, optional): identifier
            name (:obj:`str`, optional): name
            participants (:obj:`list` of :obj:`ReactionParticipant`, optional): list of participants in the reaction and
                their compartments and coefficients
            reversible (:obj:`bool`, optional): indicates if the reaction is reversible            
            cross_references (:obj:`list`, optional): list of cross references
        """
        self.id = id
        self.name = name
        self.participants = participants or None
        self.reversible = reversible
        self.cross_references = cross_references or []

    def normalize(self):
        """ Normalize the participant list by collapsing repeated participants and sorting them by their order """

        # collapse repeated participants
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
                        other_part.molecule.id == part.molecule.id and \
                        ((other_part.compartment is None and part.compartment is None) or
                            other_part.compartment.id == part.compartment.id) and \
                        numpy.sign(other_part.coefficient) == numpy.sign(part.coefficient):
                    part.coefficient += other_part.coefficient
                    self.participants.pop(i_other_part)
                else:
                    i_other_part += 1

            i_part += 1

        # order participants
        self.participants.sort(key=lambda part: part.order if part.order is not None else 1e10)
        for i_part, part in enumerate(self.participants):
            part.order = i_part

    def get_reactants(self):
        """ Get the reactants of the reaction

        Returns:
            :obj:`list` of :obj:`ReactionParticipant`: list of reactants in the reaction and
                their compartments and coefficients
        """
        return list(filter(lambda p: p.coefficient < 0, self.participants))

    def get_products(self):
        """ Get the products of the reaction

        Returns:
            :obj:`list` of :obj:`ReactionParticipant`: list of products in the reaction and
                their compartments and coefficients
        """
        return list(filter(lambda p: p.coefficient > 0, self.participants))

    def get_reactant_product_pairs(self):
        """ Get list of pairs of similar reactants and products

        Note: This requires the modeler to have ordered the reactans and products by their similarity. The modeler is required to 
        specify this pairing because it cannot easily be computed. In particular, we have tried to use Tanitomo similarity to
        predict reactant-product pairings, but this doesn't adequately capture reaction centers.

        Returns:
            :obj:`list` of :obj:`tuple` of obj:`molecule_util.Molecule`, :obj:`molecule_util.Molecule`: list of pairs of similar reactants and products
        """
        self.normalize()

        return list(zip_longest(self.get_reactants(), self.get_products()))

    def calc_reactant_product_pairs(self):
        """ Get list of pairs of similar reactants and products using a greedy algorithm.

        Returns:
            :obj:`list` of :obj:`tuple` of obj:`molecule_util.Molecule`, :obj:`molecule_util.Molecule`: list of pairs of similar reactants and products
        """
        self.normalize()

        # sort by structure to ensure result is reproducible
        key = lambda p: (len(p.molecule.structure), p.molecule.structure)
        reactants = sorted(self.get_reactants(), key=key, reverse=True)
        products = sorted(self.get_products(), key=key, reverse=True)

        # calculate similarities between each reactant and each product
        similarities = numpy.full((len(reactants), len(products)), numpy.nan)
        for i_reactant, reactant in enumerate(reactants):
            for i_product, product in enumerate(products):
                similarities[i_reactant, i_product] = reactant.molecule.get_similarity(product.molecule)

        # initialize pairs of similar reactants and products
        pairs = []

        # iteratively identify the most similar pair of reactants and products
        for i in range(min(len(reactants), len(products))):
            index = numpy.argmax(similarities)
            indices = numpy.unravel_index(index, dims=similarities.shape)
            i_reactant = indices[0]
            i_product = indices[1]
            pairs.append((reactants[i_reactant], products[i_product]))

            reactants.pop(i_reactant)
            products.pop(i_product)
            similarities = numpy.delete(similarities, i_reactant, axis=0)
            similarities = numpy.delete(similarities, i_product, axis=1)

        # unpaired products, reactants
        for reactant in reactants:
            pairs.append((reactant, None))
        for product in products:
            pairs.append((None, product))

        return pairs

    def get_ec_number(self):
        """ Get the most relevant EC number from the list of cross references

        * If the reaction has a single manually-assigned EC number, return that
        * If the reaction has multiple manually-assigned EC numbers, return an error
        * Otherwise, return the most relevant predicted EC number

        Returns:
            :obj:`str`: most relevant EC number
        """

        # return
        manual_ec_numbers = list(filter(lambda xr: xr.source == 'ec-code' and xr.assignment_method ==
                                        data_structs.CrossReferenceAssignmentMethod.manual, self.cross_references))
        if len(manual_ec_numbers) == 1:
            return manual_ec_numbers[0].id
        elif len(manual_ec_numbers) > 1:
            raise ValueError('Reaction {} has multiple EC numbers: {}'.format(self.id, ', '.join(xr.id for xr in manual_ec_numbers)))

        # find most relevant predicted EC number
        max_relevance = float('-inf')
        max_id = ''
        for xr in self.cross_references:
            if xr.assignment_method == data_structs.CrossReferenceAssignmentMethod.predicted and xr.relevance > max_relevance:
                max_relevance = xr.relevance
                max_id = xr.id

        return max_id


class ReactionParticipant(object):
    """ Represents a participant in a reaction

    Attributes:
        molecule (:obj:`molecule_util.Molecule`): molecule
        compartment (:obj:`compartment_util.Compartment`): compartment
        coefficient (:obj:`float`): coefficient of the molecule/compartment in the reaction
        order (:obj:`int`): order of the participant within the reaction (this used to determine pairings of reactants and products)
    """

    def __init__(self, molecule=None, compartment=None, coefficient=float('nan'), order=None):
        """
        Args:
            molecule (:obj:`molecule_util.Molecule`, optional): molecule
            compartment (:obj:`compartment_util.Compartment`, optional): compartment
            coefficient (:obj:`float`, optional): coefficient of the molecule/compartment in the reaction
            order (:obj:`int`, optional): order of the participant within the reaction (this used to determine pairings of reactants and products)
        """
        self.molecule = molecule
        self.compartment = compartment
        self.coefficient = coefficient
        self.order = order
