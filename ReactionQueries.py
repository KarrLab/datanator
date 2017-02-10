import openpyxl
import InchiGenerator
import QueryStringManipulator
import ECNumberFinder


class Compound:
	def __init__(self, id, inchiSmiles = "", sabioNames = []):
		self.id = id
		self.inchiSmiles = inchiSmiles
		self.sabioNames = sabioNames



class ReactionQuery:
	def __init__(self, id):
		self.id = id
		self.reactionString = ""
		self.substrates = [] #this is a list of compound objects
		self.products = [] #this is a list of compound objects
		self.numParticipants = [] #this is an array of two numbers, first number for substrates, second for products
		self.keggID = ""
		self.ECNumber = ""

	def getQueryString(self):
		subAndProd = []
		subNames = []
		prodNames = []
		searchString = ""

		numParticipants = len(self.substrates) + len(self.products)

		for compound in self.substrates:
			if len(compound.sabioNames)>0:
				subNames.append(compound.sabioNames)
		for compound in self.products:
			if len(compound.sabioNames)>0:
				prodNames.append(compound.sabioNames)

		numSabioFound = len(subNames)+len(prodNames)
		print numParticipants
		print numSabioFound
		print numSabioFound >= numParticipants - 1
		if numSabioFound >= numParticipants - 1:
			subAndProd.append(subNames)
			subAndProd.append(prodNames)
		#print subAndProd
			searchString = QueryStringManipulator.getQuerySearchString(subAndProd)
		return searchString

	def generateLiftedReactionQuery(self):#, substrates, products):
		lifted = LiftedReactionQuery(id, self.substrates, self.products)
		return lifted


#to work on: what to do in a reaction like amp+atp==>adp. In that case
#there is really 2 adp molecules, but the fields think there's only one
class LiftedReactionQuery:
	def __init__(self, id, substrates, products):
		self.id = id
		self.substrates = substrates
		self.products = products
		self.numParticipants = ""
		self.genericECNumber = self.findECNumber(self.substrates, self.products)

	def findECNumber(self, substrates, products):
		subInchiSmiles = []
		prodInchiSmiles = []

		for compound in substrates:
			subInchiSmiles.append(compound.inchiSmiles)
		for compound in products:
			prodInchiSmiles.append(compound.inchiSmiles)

		ECNum = ""
		i = 0
		while i<len(subInchiSmiles) and len(ECNum) == 0:
			ECNum = ECNumberFinder.getECNumber(subInchiSmiles, prodInchiSmiles)
			subInchiSmiles = subInchiSmiles[-1:] + subInchiSmiles[:-1]
			i += 1

		if len(products)>len(substrates) and len(ECNum) == 0:
			prodInchiSmiles = prodInchiSmiles[-1:] + prodInchiSmiles[:-1]
			ECNum = ECNumberFinder.getECNumber(subInchiSmiles, prodInchiSmiles)
			i += 1
		return ECNum

	
	def getQueryString(self):
		return ECNumberFinder.formatECForSabio(self.genericECNumber)


#Input an excel sheet
#It outputs a list of reactionQuery Objects
def generateReactionQueries(filename):
	reactionQueries = [] #this is the only output of this method. Its a list of reactionQueries
	idToCompound = {} #this will be a dict used to find compound information about each metab id
	
	wb = openpyxl.load_workbook(filename=filename)
	ws = wb.get_sheet_by_name('Metabolites')

	#this gets a dict that is used to correlate inchiString with Sabio's name for that compount
	sabioNameToInchiDict = InchiGenerator.getSabioNameToInchiDict()

	#this instantiates a compound object for each metabolite in the excel sheet
	i = 0
	while i < len(ws.columns[0]):
		id = ws.columns[0][i].value
		smilesOrInchi = ws.columns[1][i].value
		sabioNames = []
		if smilesOrInchi != None:
			genericInchi = InchiGenerator.generateGenericInchi(smilesOrInchi)
			for name in sabioNameToInchiDict:
				if sabioNameToInchiDict[name] == genericInchi:
					sabioNames.append(name)

		comp = Compound(id, smilesOrInchi, sabioNames)
		idToCompound[id] = comp
		i += 1


	#this creates a dictionary of reaction ID to reaction string
	#to work on: making it possible to only enter one row
	#to work on: make it possible not to have any reaction string data to begin with
	ws = wb.get_sheet_by_name('Reactions')
	idToReactionDict = {}
	i = 0
	while i < len(ws.columns[0]):
		id = "{}".format(ws.columns[0][i].value)
		if len("{}".format(ws.columns[1][i].value))>0:
			idToReactionDict[id] = "{}".format(ws.columns[1][i].value)
		else:
			idToReactionDict[id] = id
		i += 1


	#this is is where the ReactionQuery objects are instantiated
	#each reaction is parsed, each metabolite ID is converted to a Compound object
	#then the fields for each ReactionQuery object are set
	#finally the ReactionQuery object is added to the reactionQueries list
	for id in idToReactionDict:
		query = ReactionQuery(id)
		query.reactionString = idToReactionDict[id]

		balancedMetab = QueryStringManipulator.getParsedReaction(idToReactionDict[id])
		query.numParticipants = [len(balancedMetab[0]), len(balancedMetab[1])]

		compoundBalancedMetab = []
		for side in balancedMetab:
			compoundSide = []
			for metab in side:
				try:
					compoundSide.append(idToCompound[metab])
				except:
					newUnrecognizedCompound = Compound(metab)
					compoundSide.append(newUnrecognizedCompound)
			compoundBalancedMetab.append(compoundSide)
		query.substrates = compoundBalancedMetab[0]
		query.products = compoundBalancedMetab[1]
		reactionQueries.append(query)

	return reactionQueries





		

if __name__ == '__main__':
	filename='SmilesStuff.xlsx'
	reactionQueries = generateReactionQueries(filename)
	for reaction in reactionQueries:
		print reaction.__dict__
		print reaction.getQueryString()