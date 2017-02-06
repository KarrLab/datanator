from SabioInterface import getSabioData
import ReactionQueries
import numpy as np
from CreateExcelSheet import createExcelSheet

#to work on: filter by numparticipants
#add other search criteria
#make sure to use any data if less than (everything-1) was found. BUT!! I can use lift if all 
#smiles/inchi is present

#finish debugging the EC finder

class KineticInfo:
	def __init__(self, sabioResults, name = ""):
		self.name = name
		self.closestEntriesIDs = []
		self.closestEntries = []
		self.medianEntry = []
		self.minEntry = []
		self.maxEntry = []
		self.liftInfo = ""
		self.reactionList = sabioResults.getFieldList(sabioResults.entryList, "reactionStoichiometry")

		#should maybe change these two to only include the values of of the reactions
		#that are in the lowest species
		self.SabioReactionIDs = sabioResults.getFieldList(sabioResults.entryList, "reactionID")
		self.ECNumbers = sabioResults.getFieldList(sabioResults.entryList, "ECNumber")




		proximNums = sabioResults.getFieldList(sabioResults.entryList, "proximity")

		n = 0
		narrowedEntries = []
		while n < len(proximNums): 
			for entry in sabioResults.entryList:
				if entry.proximity == proximNums[n] and len(entry.__dict__[name])>0:
					narrowedEntries.append(entry)
			if len(narrowedEntries)>0:
				break
			else:
				n = n+1
		self.closestEntries = narrowedEntries

		orderedEntries = sorted(narrowedEntries, key=lambda entry: float(entry.__dict__[name]))#, reverse=True)
		#records closest entry Ids - this is currently not working
		self.closestEntriesIDs = sabioResults.getFieldList(orderedEntries, "entryID")

		#to work on: if there are an even number, what should the median be?
		if len(orderedEntries) >= 2:
			self.minEntry = orderedEntries[0]
			self.maxEntry = orderedEntries[len(orderedEntries)-1]

		#to work on: make sure this media is done right
		if len(orderedEntries)>2 or len(orderedEntries)==1:
			number = float(len(orderedEntries))
			number = int(np.around(number/2+.1))
			self.medianEntry = orderedEntries[number-1]

		


	#this checks 
	def hasRelevantData(self, proximLimit = 1000):
		hasRelevantData = False
		if len(self.closestEntries)>0 and self.closestEntries[0].proximity <= proximLimit:
			hasRelevantData = True
		return hasRelevantData


class VmaxInfo(KineticInfo):
	def __init__(self, sabioResults):
		KineticInfo.__init__(self, sabioResults, "vmax")


class KmInfo(KineticInfo):
	def __init__(self, sabioResults):
		KineticInfo.__init__(self, sabioResults, "km")



class FormattedData:
	def __init__(self, id):#, sabioResults):
		self.id = id
		self.reactionIDs = []
		self.KmData = None
		self.VmaxData = None


def createFormattedData(reactionQuery, species, proximLimit = 1000):
	print reactionQuery.numParticipants
	sabioResults = getSabioData(reactionQuery.getQueryString(), species, reactionQuery.numParticipants)
	formattedData = FormattedData(reactionQuery.id)

	formattedData.reactionIDs = sabioResults.getFieldList(sabioResults.entryList, "reactionID")

	liftedQuery = None
	liftedSabioResults = None
	
	#get the km data. If data is not present, then it searches for similar reactions. 
	#proximLimit is the taxonomic limit (in nodes) after which the data is presumed to be irrelevant 
	formattedData.KmData = KmInfo(sabioResults)
	hasRelevantData = formattedData.KmData.hasRelevantData(proximLimit)
	if hasRelevantData == False:
		liftedQuery = reactionQuery.generateLiftedReactionQuery()
		queryString = liftedQuery.getQueryString()
		liftedSabioResults = getSabioData(queryString, species)
		formattedData.KmData = KmInfo(liftedSabioResults)
		formattedData.KmData.liftInfo = "Lifted From {}".format(liftedQuery.genericECNumber)


	formattedData.VmaxData = VmaxInfo(sabioResults)
	hasRelevantData = formattedData.VmaxData.hasRelevantData(proximLimit)
	if hasRelevantData == False:
		if liftedQuery == None:
			liftedQuery = reactionQuery.generateLiftedReactionQuery()
			queryString = liftedQuery.getQueryString()
			liftedSabioResults = getSabioData(queryString, species)
		formattedData.VmaxData = VmaxInfo(liftedSabioResults)
		formattedData.VmaxData.liftInfo = "Lifted From {}".format(liftedQuery.genericECNumber)


	return formattedData



def main(filename, species):
	
	#this is the endgame
	formattedDataList = []
	reactionQueries = ReactionQueries.generateReactionQueries(filename)
	for reactionQuery in reactionQueries:
		formattedData = createFormattedData(reactionQuery, species)
		formattedDataList.append(formattedData)

	#print fullResponse.__dict__
	
	for formattedData in formattedDataList:
		print formattedData.__dict__
		#print formattedData.KmData.__dict__
		#print formattedData.VmaxData.__dict__
	
		#print formattedData.VmaxData.__dict__

	#now I need to make a method to make this into an excel sheet
	createExcelSheet(formattedDataList, species)


 
if __name__ == '__main__':

	filename='SmilesStuff2.xlsx'
	species = 'mycoplasma pneumoniae'
	main(filename, species)

	#number = 1.0
	#number = float(number/2)+.1
	#print number
	#print np.around(number)

	"""
	query_dict = {
				#"Organism":'"Homo sapiens"',
				"Substrate": "AMP AND ADP", 
				"Product": "ADP",
				#"Enzymename":"Adk"
				#"EnzymeType":"wildtype"
				#"SabioReactionID" : "5305"

				#"Organism":'"Homo sapiens"',
				#"Substrate": "nad",
				#"Product": "nadh"
				#"Enzymename":"Adk"
				#"EnzymeType":"wildtyp
				}
	answer =  getSabioData(query_dict)
 	fullResponse = FormattedData("blue")
	#fullResponse.VmaxData = VmaxInfo(answer)
	fullResponse.KmData = KmInfo(answer)
	"""