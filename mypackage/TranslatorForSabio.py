import QueryStringManipulator
import ECNumberFinder
import copy


#input a reactionQuery
#outputs a search string based on its substrates and products
def getSubstrateProductQueryString(reactionQuery):

	reactionQueryCopy = copy.deepcopy(reactionQuery)

	subAndProd = []
	subNames = []
	prodNames = []
	searchString = ""


	numTotalParticipants = len(reactionQueryCopy.substrates) + len(reactionQueryCopy.products)

	for compound in reactionQueryCopy.substrates:
		if len(compound.sabioNames)>0:
			subNames.append(compound.sabioNames)
	for compound in reactionQueryCopy.products:
		if len(compound.sabioNames)>0:
			prodNames.append(compound.sabioNames)

	numSabioFound = len(subNames)+len(prodNames)
	if numSabioFound >= numTotalParticipants - 1:
		subAndProd.append(subNames)
		subAndProd.append(prodNames)
	#print subAndProd
		searchString = QueryStringManipulator.getQuerySearchString(subAndProd)
	return searchString

def getGenericECQueryString(reactionQuery):
	reactionQueryCopy = copy.deepcopy(reactionQuery)
	return ECNumberFinder.formatECForSabio(reactionQueryCopy.genericECNumber)