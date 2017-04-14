""" Utilities for dealing with reactions

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan <jonrkarr@gmail.com>
:Date: 2017-04-13
:Copyright: 2017, Karr Lab
:License: MIT
"""

from . import compartment_util
from . import molecule_util
import itertools
import numpy
import re
import requests


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
                        ((other_part.compartment is None and part.compartment is None) or other_part.compartment.name == part.compartment.name) and \
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
        return list(filter(lambda p: p.coefficient < 0, self.participants))

    def get_products(self):
        """ Get the products of the reaction

        Returns:
            :obj:`list` of :obj:`ReactionParticipant`: list of products in the reaction and
                their compartments and coefficients
        """
        return list(filter(lambda p: p.coefficient > 0, self.participants))

    def get_reactant_product_pairs(self):
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


class Ezyme(object):
    """ Utilities for using Ezyme to predict EC numbers.

    See Ezyme (http://www.genome.jp/tools-bin/predict_reaction) for more information.

    Attributes:
        REQUEST_URL (:obj:`str`): URL to request Ezyme EC number prediction
        RETRIEVAL_URL (:obj:`str`): URL to retrieve Ezyme results
        EC_PREDICTION_URL (:obj:`str`): URL where predicted EC number is encoded
    """
    REQUEST_URL = 'http://www.genome.jp/tools-bin/predict_view'
    RETRIEVAL_URL = 'http://www.genome.jp/tools-bin/e-zyme2/result.cgi'
    EC_PREDICTION_URL = 'http://www.genome.jp/kegg-bin/get_htext?htext=ko01000.keg&query='

    def __init__(self):
        pass

    @classmethod
    def run(cls, reaction):
        """ Use Ezyme to predict the first three digits of the EC number of a reaction.

        Args:
            :obj:`Reaction`: reaction

        Returns:
            :obj:`list` of :obj:`EzymeResult`: ranked list of predicted EC numbers and their scores
        """
        pairs = reaction.get_reactant_product_pairs()

        reactants = []
        products = []
        for pair in pairs:
            if pair[0]:
                reactants.append(pair[0].molecule.to_mol())
            if pair[1]:
                products.append(pair[1].molecule.to_mol())

        return cls._run(reactants, products)

    @classmethod
    def _run(cls, reactants, products):
        """ Low level interface to use Ezyme to predict the first three digits of the EC number of a reaction.

        Note: The order of the compounds must be the same in the reactants and products list. If you have not manually
        ordered the reactants and products, we recommend that you use the :obj:`run` method which will
        order the reactants and products lists according to their chemical similarity.        

        Args:
            reactants (:obj:`list` of :obj:`str`): list of structures of reactants in MOL format
            products (:obj:`list` of :obj:`str`): list of structures of products in MOL format

        Returns:
            :obj:`list` of :obj:`EzymeResult`: ranked list of predicted EC numbers and their scores
        """

        # Request Ezyme to predict the EC number
        response = requests.post(cls.REQUEST_URL, data={
            'QUERY_MODE': 'MULTI',
            'S_MOLTEXT': reactants,
            'P_MOLTEXT': products,
        })
        response.raise_for_status()
        file = re.findall('<input type=hidden name=file value="(.*?)">', response.text)[0]

        # retrieve the predicted EC number(s)
        response = requests.post(cls.RETRIEVAL_URL, data={
            'file': file,
            'name': ''.join(':' + str(i) for i, c in enumerate(itertools.chain(reactants, products)))
        })
        response.raise_for_status()
        results = re.findall(
            '<td align=center><a href="{}{}">[0-9]+\.[0-9]+\.[0-9]+</a></td>\n<td align=center>([0-9\.]+)</td>'.format(
                re.escape(cls.EC_PREDICTION_URL), '([0-9]+\.[0-9]+\.[0-9]+)\.'),
            response.text, re.MULTILINE)
        return [EzymeResult(r[0], float(r[1])) for r in results]


class EzymeResult(object):
    """ Represents a predicted EC number

    Attributes:
        ec_number (:obj:`str`): EC number
        score (:obj:`float`): score
    """

    def __init__(self, ec_number, score):
        """
        Args:
            ec_number (:obj:`str`): EC number
            score (:obj:`float`): score
        """
        self.ec_number = ec_number
        self.score = score
