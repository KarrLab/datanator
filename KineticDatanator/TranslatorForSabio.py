import QueryStringManipulator
import ECNumberFinder
import copy


#input a reactionQuery
#checks to make sure that sabio recognizes at least (all-1) of them. This ensures that the scope of the seaach
#isn't too broad. If sabio does nor recognzie two or more of the compounds, a blank string is returned
#outputs a search string based on its substrates and products
def getSubstrateProductQueryString(reactionQuery):

	reactionQueryCopy = copy.deepcopy(reactionQuery)

	subAndProd = []
	subNames = []
	prodNames = []
	searchString = ""
	numTotalParticipants = len(reactionQueryCopy.substrates) + len(reactionQueryCopy.products)

	for compound in reactionQueryCopy.substrates:
		#print compound.sabioNames
		if len(compound.sabioNames)>0:
			subNames.append(compound.sabioNames)
	for compound in reactionQueryCopy.products:
		#print compound.sabioNames
		if len(compound.sabioNames)>0:
			prodNames.append(compound.sabioNames)

	numSabioFound = len(subNames)+len(prodNames)
	#print numSabioFound
	if numSabioFound >= numTotalParticipants - 1:
		subAndProd.append(subNames)
		subAndProd.append(prodNames)
	#print subAndProd
		searchString = QueryStringManipulator.getQuerySearchString(subAndProd)
	return searchString

def getGenericECQueryString(reactionQuery):
	reactionQueryCopy = copy.deepcopy(reactionQuery)
	return ECNumberFinder.formatECForSabio(reactionQueryCopy.genericECNumber)