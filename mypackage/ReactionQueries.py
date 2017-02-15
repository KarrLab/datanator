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
		self.genericECNumber = ""

	def setGenericECNumberFromEzymeAlgorithm(self):
		self.genericECNumber = ECNumberFinder.easyFindECNumber(self.substrates, self.products)


#this takes an workbook with "metabolites" as an input. It outputs an array of molecules
def generateCompounds(wb):
	sabioNameToInchiDict = InchiGenerator.getSabioNameToInchiDict()
	ws = wb.get_sheet_by_name('Metabolites')
	#this instantiates a compound object for each metabolite in the excel sheet
	compoundList = []
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
		compoundList.append(comp)
		i += 1
	return compoundList


#Input an excel sheet data object (from openpyxl)
#It outputs a list of reactionQuery Objects
def generateReactionQueries(excelSheetObject):
	reactionQueries = [] #this is the only output of this method. Its a list of reactionQueries
	idToCompound = {} #this will be a dict used to find compound information about each metab id
	wb = excelSheetObject

	#this gets a dict that is used to correlate inchiString with Sabio's name for that compount
	sabioNameToInchiDict = InchiGenerator.getSabioNameToInchiDict()

	#this instantiates a compound object for each metabolite in the excel sheet
	compoundList = generateCompounds(wb)
	#this creates a dictionary of ID to compound to make it quick to find the compound object 
	#with the compound ID. This works because each compound ID must be unique
	for comp in compoundList:
		idToCompound[comp.id] = comp


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