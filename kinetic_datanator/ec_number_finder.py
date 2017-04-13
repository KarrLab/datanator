""" 
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-06
:Copyright: 2017, Karr Lab
:License: MIT
"""

import re
import requests
from itertools import chain
from kinetic_datanator.util import compound_util
from requests.packages.urllib3.exceptions import InsecureRequestWarning
requests.packages.urllib3.disable_warnings(InsecureRequestWarning)


def easy_find_ec_number(substrates, products):
    # Unlike this method, predict_ec_number requires a precise matchup of substrates to products
    # therefore, this method - easy_find_ec_number - tries all the combinations in predict_ec_number
    # its basically an add-on to predict_ec_number that makes using it much easier

    ec_num = ""

    # first check to make sure that each compound has a smiles or an inchi
    # if they do not, return an empty string
    every_compound_has_structural_info = True

    for compound in substrates:
        if len(compound.inchi_smiles) == 0:
            every_compound_has_structural_info = False
    for compound in products:
        if len(compound.inchi_smiles) == 0:
            every_compound_has_structural_info = False

    if every_compound_has_structural_info == True:
        sub_inchi_smiles = []
        prod_inchi_smiles = []

        for compound in substrates:
            sub_inchi_smiles.append(compound.inchi_smiles)
        for compound in products:
            prod_inchi_smiles.append(compound.inchi_smiles)

        i = 0
        while i < len(sub_inchi_smiles) and len(ec_num) == 0:
            ec_num = predict_ec_number(sub_inchi_smiles, prod_inchi_smiles)
            sub_inchi_smiles = sub_inchi_smiles[-1:] + sub_inchi_smiles[:-1]
            i += 1

        if len(products) > len(substrates) and len(ec_num) == 0:
            prod_inchi_smiles = prod_inchi_smiles[-1:] + prod_inchi_smiles[:-1]
            ec_num = predict_ec_number(sub_inchi_smiles, prod_inchi_smiles)
            i += 1

    return ec_num

EZYME_REQUEST_URL = 'http://www.genome.jp/tools-bin/predict_view'
# :obj:`str`: URL to request Ezyme EC number prediction

EZYME_RETRIEVAL_URL = 'http://www.genome.jp/tools-bin/e-zyme2/result.cgi'
# :obj:`str`: URL to retrieve Ezyme results


def predict_ec_number(substrates, products):
    """ Predict the first three digits of the EC number of a reaction using Ezyme.

    See Ezyme (http://www.genome.jp/tools-bin/predict_reaction) for more information.

    Args:
        substrates (:obj:`list` of :obj:`str`): list of InChI structures of substrates
        products (:obj:`list` of :obj:`str`): list of InChI structures of products

    Returns:
        :obj:`str`: predicted EC number
    """

    # Request Ezyme to predict the EC number
    response = requests.post(EZYME_REQUEST_URL, data={
        'QUERY_MODE': 'MULTI',
        'S_MOLTEXT': [compound_util.Compound(c).to_mol() for c in substrates],
        'P_MOLTEXT': [compound_util.Compound(c).to_mol() for c in products],
    })
    response.raise_for_status()
    file = re.findall('<input type=hidden name=file value="(.*?)">', response.text)[0]

    # retrieve the predicted EC number(s)
    response = requests.post(EZYME_RETRIEVAL_URL, data={
        'file': file,
        'name': ''.join(':' + str(i) for i, c in enumerate(chain(substrates, products)))
    })
    response.raise_for_status()
    ec_numbers = re.findall(
        '<a href="http://www\.genome\.jp/kegg\-bin/get_htext\?htext=ko01000\.keg&query=([0-9]\.[0-9]\.[0-9])\.">', response.text)
    if ec_numbers:
        return ec_numbers[0]
    return ''
