""" 
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-06
:Copyright: 2017, Karr Lab
:License: MIT
"""

from . import ec_number_finder
from . import query_string_manipulator
import copy

def getSubstrateProductQueryString(reaction_query):
	#input a reaction_query
	#checks to make sure that sabio recognizes at least (all-1) of them. This ensures that the scope of the seaach
	#isn't too broad. If sabio does nor recognzie two or more of the compounds, a blank string is returned
	#outputs a search string based on its substrates and products

	reaction_queryCopy = copy.deepcopy(reaction_query)

	subAndProd = []
	subNames = []
	prodNames = []
	searchString = ""
	numTotalParticipants = len(reaction_queryCopy.substrates) + len(reaction_queryCopy.products)

	for compound in reaction_queryCopy.substrates:
		#print(compound.sabioNames)
		if len(compound.sabioNames)>0:
			subNames.append(compound.sabioNames)
	for compound in reaction_queryCopy.products:
		#print(compound.sabioNames)
		if len(compound.sabioNames)>0:
			prodNames.append(compound.sabioNames)

	numSabioFound = len(subNames)+len(prodNames)
	#print(numSabioFound)
	if numSabioFound >= numTotalParticipants - 1:
		subAndProd.append(subNames)
		subAndProd.append(prodNames)
	#print(subAndProd)
		searchString = query_string_manipulator.getQuerySearchString(subAndProd)
	return searchString

def getGenericECQueryString(reaction_query):
	reaction_queryCopy = copy.deepcopy(reaction_query)
	return format_ec_number_for_sabio(reaction_queryCopy.predicted_ec_number)

def format_ec_number_for_sabio(ec_number, number_four_digit_ec_numbers=100):
    """
    Args:
        ec_number (:obj:`str`): EC number
        number_four_digit_ec_numbers (:obj:`int`): number of four digit EC numbers to generate
    """

    if not ec_number:
        return ''

    if ec_number[-2:] == '.-':
        ec_number = ec_number[0:-2]

    if ec_number.count('.') == 2:
        four_digit_ec_numbers = ['"{}.{}"'.format(ec_number, i) for i in range(1, number_four_digit_ec_numbers + 1)]
        return 'ECNumber: ({})'.format(' OR '.join(four_digit_ec_numbers))

    return 'ECNumber: "{}"'.format(ec_number)