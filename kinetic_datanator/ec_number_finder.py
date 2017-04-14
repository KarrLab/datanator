""" 
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-06
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.util import molecule_util
from kinetic_datanator.util import reaction_util


def predict_ec_number(reactants, products):
    """ Unlike this method, Ezyme.run requires a precise matchup of reactants to products
    therefore, this method - predict_ec_number - tries all the combinations in Ezyme.run
    its basically an add-on to Ezyme.run that makes using it much easier

    Args:
        reactants (:obj:`list` of :obj:`str`): list of structures of reactants in InChI or canonical SMILES format
        products (:obj:`list` of :obj:`str`): list of structures of products in InChI or canonical SMILES format

    Returns:
        :obj:`list` of :obj:`str`: ranked list of predicted EC numbers
    """

    parts = []

    for r in reactants:
        parts.append(reaction_util.ReactionParticipant(
            molecule=molecule_util.Molecule(r.inchi_smiles, name=r.id),
            coefficient=-1,
        ))

    for p in products:
        parts.append(reaction_util.ReactionParticipant(
            molecule=molecule_util.Molecule(p.inchi_smiles, name=p.id),
            coefficient=1,
        ))

    rxn = reaction_util.Reaction(parts)
    ec_numbers = reaction_util.Ezyme.run(rxn)
    if ec_numbers:
        return ec_numbers[0].ec_number
    return ''
