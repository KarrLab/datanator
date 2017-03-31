from SabioInterface import getSabioData
import ReactionQueries
import numpy as np
from CreateExcelSheet import createExcelSheet
import openpyxl
import TranslatorForSabio
import logging
import os
#add other search criteria
#make sure to only use lift if smiles/inchi is present!!!!
#also, what should median be?
#make sure i fix ecnumber in cases where we have (2)ADP or somethign like that


#finish debugging the EC finder

"""
.. autoclass:: Noodle
   :members: eat, slurp

   .. method:: boil(time=10)

      Boil the noodle *time* minutes
"""


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


	#this checks whether the current kinetic information is biologically useful. This is helpful to know
	#because if the data is not useful, then the use can try a more general search
	def hasRelevantData(self, proximLimit = 1000):
		hasRelevantData = False
		if len(self.closestEntries)>0 and (self.closestEntries[0].proximity <= proximLimit):
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
	#print "new"
	#for thing in reactionQuery.substrates:
	#	print thing.sabioNames
	#print (reactionQuery.id)
	#print (reactionQuery.numParticipants)
	searchString = TranslatorForSabio.getSubstrateProductQueryString(reactionQuery)
	logging.info("Datanator: Search String Found - {}".format(searchString))
	#print(searchString)
	#print searchString

	if len(searchString)>0:
		searchString = defaultValues+searchString
	logging.info("Datanator: About to check Sabio with regular search string")
	sabioResults = getSabioData(searchString, species, reactionQuery.numParticipants)
	#instantiate a formatted data object with the reactionID
	formattedData = FormattedData(reactionQuery.id)
	formattedData.reactionIDs = sabioResults.getFieldList(sabioResults.entryList, "reactionID")

	liftedSabioResults = None
	
	#get the km data. If data is not present, then it searches for similar reactions. 
	#proximLimit is the taxonomic limit (in nodes) after which the data is presumed to be irrelevant 
	formattedData.KmData = KmInfo(sabioResults)
	hasRelevantData = formattedData.KmData.hasRelevantData(proximLimit)
	logging.info("Datanator: Has relevant data for km - {}".format(hasRelevantData))
	if hasRelevantData == False:
		#we are about to set the genericECNumber. However, first we need ot make sure
		#the user didn't provide their own. If the user provided their own, that genericECNumber is used instead
		if len(reactionQuery.genericECNumber)==0:
			reactionQuery.setGenericECNumberFromEzymeAlgorithm()
		queryString = TranslatorForSabio.getGenericECQueryString(reactionQuery)
		if len(queryString)>0:
			queryString = defaultValues + queryString
		else:
			logging.error("Datanator: Tried to lift from EC, but did not work")
		liftedSabioResults = getSabioData(queryString, species)
		formattedData.KmData = KmInfo(liftedSabioResults)
		if len(reactionQuery.genericECNumber)>0:
			formattedData.KmData.liftInfo = "Lifted From {}".format(reactionQuery.genericECNumber)
			logging.info("Datanator: Lifted From {} for Km".format(reactionQuery.genericECNumber))


	formattedData.VmaxData = VmaxInfo(sabioResults)
	hasRelevantData = formattedData.VmaxData.hasRelevantData(proximLimit)
	logging.info("Datanator: Has relevant data for vmax - {}".format(hasRelevantData))
	if hasRelevantData == False:
		#check if liftesSabioResults has already been searched. If it has, no need to search again
		if liftedSabioResults == None:
			if len(reactionQuery.genericECNumber)==0:
				reactionQuery.setGenericECNumberFromEzymeAlgorithm()
			queryString = TranslatorForSabio.getGenericECQueryString(reactionQuery)
			if len(queryString)>0:
				queryString = defaultValues + queryString
			liftedSabioResults = getSabioData(queryString, species)
		formattedData.VmaxData = VmaxInfo(liftedSabioResults)
		if len(reactionQuery.genericECNumber)>0:
			formattedData.VmaxData.liftInfo = "Lifted From {}".format(reactionQuery.genericECNumber)
			logging.info("Datanator: Lifted From {} for Vmax".format(reactionQuery.genericECNumber))

	return formattedData



def getKineticData(inputFilename, outputFilename, species, tempRange = [15, 40], enzymeType = "wildtype", phRange = [5,9], proximLimit=1000):

	logging.basicConfig(filename=os.path.join('.', 'logging.log'), filemode='w', level=logging.INFO)
	logging.info('Datanator: Started')

	defaultValues = "enzymeType:{} AND TemperatureRange:[{} TO {}] AND pHValueRange:[{} TO {}] AND ".format(enzymeType, tempRange[0], tempRange[1], phRange[0], phRange[1])
	
	file = open("Errors.txt", "w")
	file.write("")
	file.close()

	
	#this is the endgame
	formattedDataList = []

	excelSheetObject = openpyxl.load_workbook(filename=inputFilename)
	reactionQueries = ReactionQueries.generateReactionQueries(excelSheetObject)
	for reactionQuery in reactionQueries:
		try:
			logging.info('\n')
			logging.info('Datanator: Started analyzing {}'.format(reactionQuery.id))
			formattedData = createFormattedData(reactionQuery, species, defaultValues, proximLimit)
			formattedDataList.append(formattedData)
		except:
			file = open("Errors.txt", "a")
			file.write("{}	{}".format(reactionQuery.id, reactionQuery.__dict__) + "\n")
			logging.error("Datanator: {} caused an error and did not work".format(reactionQuery.id))
	createExcelSheet(outputFilename, formattedDataList, species)
	logging.info('Datanator: Finished')

def getKineticDataFromDjango(inputFile, species, tempRange = [15, 40], enzymeType = "wildtype", phRange = [5,9], proximLimit=1000):

	logging.basicConfig(filename=os.path.join('.', 'logging.log'), filemode='w', level=logging.INFO)
	logging.info('Datanator: Started')

	defaultValues = "enzymeType:{} AND TemperatureRange:[{} TO {}] AND pHValueRange:[{} TO {}] AND ".format(enzymeType, tempRange[0], tempRange[1], phRange[0], phRange[1])
	
	file = open("Errors.txt", "w")
	file.write("")
	file.close()

	
	#this is the endgame
	formattedDataList = []

	excelSheetObject = inputFile
	reactionQueries = ReactionQueries.generateReactionQueries(excelSheetObject)
	for reactionQuery in reactionQueries:
		try:
			logging.info('\n')
			logging.info('Datanator: Started analyzing {}'.format(reactionQuery.id))
			formattedData = createFormattedData(reactionQuery, species, defaultValues, proximLimit)
			formattedDataList.append(formattedData)
		except:
			file = open("Errors.txt", "a")
			file.write("{}	{}".format(reactionQuery.id, reactionQuery.__dict__) + "\n")
			logging.error("Datanator: {} caused an error and did not work".format(reactionQuery.id))
	#createExcelSheet(outputFilename, formattedDataList, species)
	logging.info('Datanator: Finished')
	return formattedDataList
	
if __name__ == '__main__':

	inputFilename='SmilesStuff.xlsx'
	outputFilename = "THEDATATRYAGAIN.xlsx"
	species = 'mycoplasma pneumoniae'
	
	getKineticData(inputFilename, outputFilename, species, proximLimit = 8)
