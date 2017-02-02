import openpyxl
from InchiGenerator import generateGenericInchi
from InchiGenerator import getSabioNameToInchiDict


class Compound:
	def __init__(self, ID, inchiSmiles, sabioNames):
		self.ID = ID
		self.inchiSmiles = inchiSmiles
		self.sabioNames = sabioNames


class Reaction:
	def __init__(self, ID):
		self.ID = ID
		self.substrates = []
		self.products = []
		self.numParticipants = ""
		self.keggID = ""
		self.ECNumber = ""

	def generateQueryString():
		return ""

	def generateLiftQueryString():
		return ""



#make a method that takes in the excel sheet, and outputs a list of all the reactions.
#internally, this method should have an ID to compound Dict, but it shouldnt share it. 

#takes in an excel sheet. ouputs reation objects
#def generateReactionQueries():


if __name__ == '__main__':

	idToCompound = {}
	wb = openpyxl.load_workbook(filename='SmilesStuff2.xlsx')
	ws = wb.get_sheet_by_name('Metabolites')


	sabioNameToInchiDict = getSabioNameToInchiDict()

	metabToSmiles = {}

	i = 0
	while i < len (ws.columns[0]):
		id = ws.columns[0][i].value
		smilesOrInchi = ws.columns[1][i].value
		sabioNames = []
		if smilesOrInchi != None:
			genericInchi = generateGenericInchi(smilesOrInchi)
			for name in sabioNameToInchiDict:
				if sabioNameToInchiDict[name] == genericInchi:
					sabioNames.append(name)

		comp = Compound(id, smilesOrInchi, sabioNames)
		idToCompound[id] = comp
		i += 1

	for id in idToCompound:
		print idToCompound[id].__dict__




		#metabToSmiles[thing.value] = ws.columns[1][i].value
		#i += 1

	metabToInchi = {}
	for metab in metabToSmiles:
		metabToInchi[metab] = generateInchi(metabToSmiles[metab])



	ws = wb.get_sheet_by_name('Reactions')
	reactionArray = []
	#for thing in ws.columns[0]:
	#	print thing.value