import openpyxl
from InchiGenerator import generateGenericInchi
from InchiGenerator import getSabioNameToInchiDict


class Compound:
	def __init__(self, ID, inchiSmiles, sabioNames):
		self.ID = ID
		self.inchiSmiles = inchiSmiles
		self.sabioNames = sabioNames


class ReactionQuery:
	def __init__(self, ID):
		self.ID = ID
		self.substrates = []
		self.products = []
		self.numParticipants = ""
		self.keggID = ""
		self.ECNumber = ""

	def generateQueryString():
		return ""

	def generateLiftedReactionQuery():
		return ""

class LiftedReactionQuery:
	def __init__(self, ID):
		self.ID = ID
		self.substrates = []
		self.products = []
		self.numParticipants = ""
		self.genericECNumber = ""
	
	def generateQueryString():
		return ""




#make a method that takes in the excel sheet, and outputs a list of all the reactions.
#internally, this method should have an ID to compound Dict, but it shouldnt share it. 

#takes in an excel sheet. ouputs reation objects
def generateReactionQueries(filename):
	idToCompound = {} #this will be a dict used to find compound information about each metab id
	wb = openpyxl.load_workbook(filename=filename)
	ws = wb.get_sheet_by_name('Metabolites')

	#this gets a dict that is used to correlate inchiString with Sabio's name for that compount
	sabioNameToInchiDict = getSabioNameToInchiDict()

	#this instantiates a compound object for each metabolite in the excel sheet
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

	


	ws = wb.get_sheet_by_name('Reactions')
	reactionArray = []
	for thing in ws.columns[0]:
		reactionArray.append(thing.value)
	print reactionArray



if __name__ == '__main__':
	filename='SmilesStuff2.xlsx'
	generateReactionQueries(filename)