from SabioInterface import getSabioData
import ReactionQueries
import numpy as np
from CreateExcelSheet import createExcelSheet

#to work on: filter by numparticipants
#add other search criteria
#make sure to use any data if less than (everything-1) was found. BUT!! I can use lift if all 
#smiles/inchi is present
#also, what should median be?
#add in default setting for temp, ph, and enzyme

#finish debugging the EC finder

class KineticInfo:
	def __init__(self, sabioResults, name = ""):
		self.name = name
		self.closestEntryIDs = []
		self.closestEntries = []
		self.closestValues = []
		self.medianEntry = None
		self.minEntry = None
		self.maxEntry = None
		self.liftInfo = "Lift Not Used"
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
		self.closestEntryIDs = sabioResults.getFieldList(orderedEntries, "entryID")
		self.closestValues = sabioResults.getFieldList(orderedEntries, self.name)

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


def createFormattedData(reactionQuery, species, defaultValues, proximLimit = 1000):
	#print reactionQuery.id
	#print reactionQuery.numParticipants
	searchString =reactionQuery.getQueryString()
	if len(searchString)>0:
		searchString = defaultValues+searchString
	print searchString
	sabioResults = getSabioData(searchString, species, reactionQuery.numParticipants)
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
		if len(liftedQuery.genericECNumber)>0:
			formattedData.KmData.liftInfo = "Lifted From {}".format(liftedQuery.genericECNumber)


	formattedData.VmaxData = VmaxInfo(sabioResults)
	hasRelevantData = formattedData.VmaxData.hasRelevantData(proximLimit)
	if hasRelevantData == False:
		if liftedQuery == None:
			liftedQuery = reactionQuery.generateLiftedReactionQuery()
			queryString = liftedQuery.getQueryString()
			liftedSabioResults = getSabioData(queryString, species)
		formattedData.VmaxData = VmaxInfo(liftedSabioResults)
		if len(liftedQuery.genericECNumber)>0:
			formattedData.VmaxData.liftInfo = "Lifted From {}".format(liftedQuery.genericECNumber)


	return formattedData



def main(filename, species, tempRange = [30, 40], enzymeType = "wildtype", phRange = [5,9], proximLimit=1000):

	defaultValues = "enzymeType:{} AND TemperatureRange:[{} TO {}] AND pHValueRange:[{} TO {}] AND ".format(enzymeType, tempRange[0], tempRange[1], phRange[0], phRange[1])


	file = open("Errors.txt", "w")
	file.write("")
	#this is the endgame
	formattedDataList = []

	reactionQueries = ReactionQueries.generateReactionQueries(filename)
	for reactionQuery in reactionQueries:
		try:
			formattedData = createFormattedData(reactionQuery, species, defaultValues, proximLimit)
			formattedDataList.append(formattedData)
		except:
			file = open("Errors.txt", "a")
			file.write("{}	{}".format(reactionQuery.id, reactionQuery.__dict__) + "\n")
	createExcelSheet(formattedDataList, species)


if __name__ == '__main__':

	filename='SmilesStuff.xlsx'
	species = 'mycoplasma pneumoniae'
	main(filename, species, proximLimit = 8)
