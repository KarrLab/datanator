""" Utilities for dealing with reactions

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan <jonrkarr@gmail.com>
:Date: 2017-04-13
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import observation
from kinetic_datanator.util import molecule_util
import numpy


def calc_reactant_product_pairs(reaction):
    """ Get list of pairs of similar reactants and products using a greedy algorithm.

    Args:
        reaction (:obj:`observation.Reaction`): reaction

    Returns:
        :obj:`list` of :obj:`tuple` of obj:`observation.Specie`, :obj:`observation.Specie`: list of pairs of similar reactants and products
    """
    participants = reaction.get_ordered_participants()
    reactants = list(filter(lambda p: p.coefficient < 0, participants))
    products = list(filter(lambda p: p.coefficient > 0, participants))

    # sort by structure to ensure result is reproducible
    key = lambda p: (len(p.specie.structure), p.specie.structure)
    reactants = sorted(reactants, key=key, reverse=True)
    products = sorted(products, key=key, reverse=True)

    # create :obj:`molecule_util.Molecule` objects for each reactant and product
    reactant_mols = [molecule_util.Molecule(structure=reactant.specie.structure) for reactant in reactants]
    product_mols = [molecule_util.Molecule(structure=product.specie.structure) for product in products]

    # calculate similarities between each reactant and each product
    similarities = numpy.full((len(reactants), len(products)), numpy.nan)
    for i_reactant, reactant in enumerate(reactant_mols):
        for i_product, product in enumerate(product_mols):
            similarities[i_reactant, i_product] = reactant.get_similarity(product)

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
