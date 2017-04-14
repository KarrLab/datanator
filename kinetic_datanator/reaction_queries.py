""" 
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-06
:Copyright: 2017, Karr Lab
:License: MIT
"""

from . import data_structs
from . import inchi_generator
from . import io
from . import query_string_manipulator
from .util import molecule_util
from .util import reaction_util


class ReactionQuery:

    def __init__(self, id):
        self.id = id
        self.substrates = []  # this is a list of compound objects
        self.products = []  # this is a list of compound objects
        self.num_participants = []  # this is an array of two numbers, first number for substrates, second for products
        self.keggID = ""
        self.ec_number = ""
        self.predicted_ec_number = ""

    def set_predicted_ec_number_from_ezyme_algorithm(self):
        """ Predict EC number and store the value

        Returns:
            :obj:`str`: top-predicted EC number
        """
        parts = []

        for r in self.substrates:
            parts.append(reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(r.inchi_smiles, name=r.id),
                coefficient=-1,
            ))

        for p in self.products:
            parts.append(reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(p.inchi_smiles, name=p.id),
                coefficient=1,
            ))

        rxn = reaction_util.Reaction(parts)
        ec_numbers = reaction_util.Ezyme.run(rxn)
        if ec_numbers:
            ec_number = ec_numbers[0].ec_number
        else:
            ec_number = ''

        self.predicted_ec_number = ec_number


def generate_reaction_queries(filename):
    # Input an excel sheet data object (from openpyxl)
    # It outputs a list of reaction_query Objects

    # this gets a dict that is used to correlate inchiString with Sabio's name for that compount
    sabioNameToInchiDict = inchi_generator.getSabioNameToInchiDict()

    # this instantiates a compound object for each metabolite in the excel sheet
    compounds = io.InputReader.read_compounds(filename)
    reactions = io.InputReader.read_reactions(filename)

    # dictionary with keys = compound ids and values = compounds
    compound_dict = {c.id: c for c in compounds}

    # this is is where the ReactionQuery objects are instantiated
    # each reaction is parsed, each metabolite ID is converted to a Compound object
    # then the fields for each ReactionQuery object are set
    # finally the ReactionQuery object is added to the reaction_queries list
    reaction_queries = []  # this is the only output of this method. Its a list of reaction_queries
    for rxn in reactions:
        query = ReactionQuery(rxn['id'])
        query.num_participants = [len(rxn['stoichiometry'][0]), len(rxn['stoichiometry'][1])]

        compoundBalancedMetab = []
        for side in rxn['stoichiometry']:
            compoundSide = []
            for metab in side:
                if metab in compound_dict:
                    compoundSide.append(compound_dict[metab])
                else:
                    undefined_compound = data_structs.Compound(metab)
                    compoundSide.append(undefined_compound)
            compoundBalancedMetab.append(compoundSide)
        query.substrates = compoundBalancedMetab[0]
        query.products = compoundBalancedMetab[1]
        reaction_queries.append(query)

    return reaction_queries
