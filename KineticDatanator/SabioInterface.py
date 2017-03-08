from TaxonFinder import getTaxonomicDistance
import requests 
import logging
from requests.packages.urllib3.exceptions import InsecureRequestWarning
requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

#to do later. Turn the data structures into a generic kinetic information gathering framework. 
#or perhaps I can use this for everything. Meaning, every database is going to be built on entries
#and upon totals. So what I could do is that I can have generic "database-dataStructures". This will be 
#Entries, and then total results which contain the results. 

class Entry:

	def __init__(self, textList):
		self.entryID = ""
		self.vmax = ""
		self.km = ""
		self.species = ""
		self.ECNumber = ""
		self.reactionID = ""
		self.proximity = ""
		self.reactionStoichiometry = ""

		line = textList[0].split("\t")

		self.species = line[0]
		self.reactionID = line[1]
		self.reactionStoichiometry = line[2]
		self.numParticipants = self.parseSabioReactionString(self.reactionStoichiometry)
		self.ECNumber = line[3]
		self.entryID = line[4]

		kineticInfo = []
		for line in textList:
			kineticInfo.append(line[line.find(self.entryID)+6:])

		for line in kineticInfo:
			if line[0:4] == "Vmax":
				parsed = line.split("\t")
				#print parsed
				#print self.entryID
				#if parsed[5] == 'mol*s^(-1)*g^(-1)' or parsed[5] == 'mol*g^(-1)*s^(-1)':
				self.vmax = parsed[2]
			if line[0:2]=="Km":
				parsed = line.split("\t")
				#print parsed
				#print self.entryID
				#if parsed[5] == "M":
				self.km = parsed[2]

	def parseSabioReactionString(self, reactionString):

		bothSides = reactionString.split(" = ")
		right = bothSides[0]
		left = bothSides[1]

		substrates = right.split(" + ")

		products = left.split(" + ")

		#remove the hydrogens to make the substrate-product count uniform
		for entry in substrates:
			if entry == "H+" or entry == "Hydrogen" or entry == "Hyrogen Ion" or entry == "H":
				substrates.remove(entry)

		for entry in products:
			if entry == "H+" or entry == "Hydrogen" or entry == "Hyrogen Ion" or entry == "H":
				products.remove(entry)

		balanced = []
		balanced.append(len(substrates))
		balanced.append(len(products))
		return balanced




class TotalResult:
	def __init__(self, sabioFile):
		self.entryList = []

		if len(sabioFile)>0:
			splitSabioFile = sabioFile.split("\n")
			entryTextArrays = []
			i = splitSabioFile[1].split("\t")[4]
			textArray = []

			#this code ensures that all the lines that relate to a single entryID get passed
			for line in splitSabioFile:
				if line[:8] != "Organism" and len(line)>0:
					if i == line.split("\t")[4]:
						textArray.append(line)
					else:
						entryTextArrays.append(textArray)
						textArray = []
						textArray.append(line)
						i = line.split("\t")[4]
			entryTextArrays.append(textArray)

			for entryData in entryTextArrays:
				entry = Entry(entryData)
				self.entryList.append(entry)



	def getFieldList(self, entryList, field):
		values = []
		for entry in entryList:
			values.append(entry.__dict__[field])
		orderedValues = sorted(set(values))
		return orderedValues

	def narrowByNums(self, numParticipants):
		i = 0
		for entry in self.entryList:
			if entry.numParticipants != numParticipants:
				self.entryList.remove(entry)


			#get the proximity for each entry
			#for entry in self.entryList:
			#	entry.proximity = getTaxonomicDistance('mycoplasma pneumoniae', entry.species)



#takes in a search dictionary or search string. Returns a TotalResult object if something is found. 
#If nothing is found, it returns "No results found for query"
def getSabioData(query_dict, baseSpecies, numParticipants = []):

	logging.info("Sabio Interface: Looking for query in Sabio")
	if len(query_dict)==0:
		logging.info("Sabio Interface: Sabio failed to respond - no search string entered")
		return TotalResult("") 

	ENTRYID_QUERY_URL = 'http://sabiork.h-its.org/sabioRestWebServices/searchKineticLaws/entryIDs' 
	#ENTRYID_QUERY_URL = 'http://sabiork.h-its.org/sabioRestWebServices/searchKineticLaws/sbml'
	PARAM_QUERY_URL = 'http://sabiork.h-its.org/entry/exportToExcelCustomizable' 

	# ask SABIO-RK for all EntryIDs matching a query 
	if isinstance(query_dict, dict):	
		query_string = ' AND '.join(['{}:{}'.format(k,v) for k,v in query_dict.items()]) 
	if isinstance(query_dict, str):
		query_string = query_dict
	query = {'format':'txt', 'q':query_string} 


	request = requests.get(ENTRYID_QUERY_URL, params = query) 
	#print request
	request.raise_for_status() # raise if 404 error 

	if request.text == "No results found for query":
		logging.info("Sabio Interface: Sabio failed to respond - No results found for query")
		return TotalResult("")
	# each entry is reported on a new line

	entryIDs = [int(x) for x in request.text.strip().split('\n')]
	print("{} matching entriess found".format(len(entryIDs)))
	logging.info("Sabio Interface: Sabio found {} results for query".format(len(entryIDs)))

	# encode next request, for parameter data given entry IDs 
	data_field = {'entryIDs[]': entryIDs} 
	query = {'format':'tsv', 'fields[]':['Organism',"SabioReactionID", 'Reaction','ECNumber','EntryID',  'Parameter']}


	request = requests.post(PARAM_QUERY_URL, params=query, data=data_field) 
	request.raise_for_status()
	textFile = request.text
	#print textFile
	resultObject = TotalResult(textFile)

	for entry in resultObject.entryList:
		entry.proximity = getTaxonomicDistance(baseSpecies, entry.species)

	if len(numParticipants)>0:
		resultObject.narrowByNums(numParticipants)


	return resultObject








if __name__ == '__main__':


	
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
	#answer =  getSabioData(query_dict)

	#print answer
	#blue = TotalResult(answer)
	
	searchString = """((Substrate:"Adenosine 3',5'-bisphosphate") AND (Substrate:"H2O" OR Substrate:"OH-")) AND ((Product:"AMP" OR Product:"Adenine-9-beta-D-arabinofuranoside 5'-monophosphate") AND (Product:"Dihydrogen phosphate" OR Product:"Phosphate"))"""
	searchString = """enzymeType:wildtype AND TemperatureRange:[30 TO 40] AND pHValueRange:[5 TO 9] AND ((Substrate:"Glyceraldehyde 3-phosphate" OR Substrate:"L-Glyceraldehyde 3-phosphate" OR Substrate:"Glycerone phosphate" OR Substrate:"D-Glyceraldehyde 3-phosphate") AND (Substrate:"D-Sedoheptulose 7-phosphate" OR Substrate:"Sedoheptulose 1-phosphate" OR Substrate:"Sedoheptulose 7-phosphate")) AND ((Product:"L-Xylulose 1-phosphate" OR Product:"D-Ribulose 5-phosphate" OR Product:"D-Xylose 5-phosphate" OR Product:"Ribose 5-phosphate" OR Product:"D-Arabinose 5-phosphate" OR Product:"D-Ribose 5-phosphate" OR Product:"D-Xylulose 1-phosphate" OR Product:"L-Xylulose 5-phosphate" OR Product:"Ribulose 5-phosphate" OR Product:"L-Ribulose 5-phosphate" OR Product:"Arabinose 5-phosphate" OR Product:"D-Xylulose 5-phosphate") AND (Product:"L-Xylulose 1-phosphate" OR Product:"D-Ribulose 5-phosphate" OR Product:"D-Xylose 5-phosphate" OR Product:"Ribose 5-phosphate" OR Product:"D-Arabinose 5-phosphate" OR Product:"D-Ribose 5-phosphate" OR Product:"D-Xylulose 1-phosphate" OR Product:"L-Xylulose 5-phosphate" OR Product:"Ribulose 5-phosphate" OR Product:"L-Ribulose 5-phosphate" OR Product:"Arabinose 5-phosphate" OR Product:"D-Xylulose 5-phosphate"))"""

	searchString = """ECNumber: ("1.2.7.0" OR "1.2.7.1" OR "1.2.7.2" OR "1.2.7.3" OR "1.2.7.4" OR "1.2.7.5" OR "1.2.7.6" OR "1.2.7.7" OR "1.2.7.8" OR "1.2.7.9" OR "1.2.7.10" OR "1.2.7.11" OR "1.2.7.12" OR "1.2.7.13" OR "1.2.7.14" OR "1.2.7.15" OR "1.2.7.16" OR "1.2.7.17" OR "1.2.7.18" OR "1.2.7.19" OR "1.2.7.20" OR "1.2.7.21" OR "1.2.7.22" OR "1.2.7.23" OR "1.2.7.24" OR "1.2.7.25" OR "1.2.7.26" OR "1.2.7.27" OR "1.2.7.28" OR "1.2.7.29" OR "1.2.7.30" OR "1.2.7.31" OR "1.2.7.32" OR "1.2.7.33" OR "1.2.7.34" OR "1.2.7.35" OR "1.2.7.36" OR "1.2.7.37" OR "1.2.7.38" OR "1.2.7.39" OR "1.2.7.40" OR "1.2.7.41" OR "1.2.7.42" OR "1.2.7.43" OR "1.2.7.44" OR "1.2.7.45" OR "1.2.7.46" OR "1.2.7.47" OR "1.2.7.48" OR "1.2.7.49" OR "1.2.7.50" OR "1.2.7.51" OR "1.2.7.52" OR "1.2.7.53" OR "1.2.7.54" OR "1.2.7.55" OR "1.2.7.56" OR "1.2.7.57" OR "1.2.7.58" OR "1.2.7.59" OR "1.2.7.60" OR "1.2.7.61" OR "1.2.7.62" OR "1.2.7.63" OR "1.2.7.64" OR "1.2.7.65" OR "1.2.7.66" OR "1.2.7.67" OR "1.2.7.68" OR "1.2.7.69" OR "1.2.7.70" OR "1.2.7.71" OR "1.2.7.72" OR "1.2.7.73" OR "1.2.7.74" OR "1.2.7.75" OR "1.2.7.76" OR "1.2.7.77" OR "1.2.7.78" OR "1.2.7.79" OR "1.2.7.80" OR "1.2.7.81" OR "1.2.7.82" OR "1.2.7.83" OR "1.2.7.84" OR "1.2.7.85" OR "1.2.7.86" OR "1.2.7.87" OR "1.2.7.88" OR "1.2.7.89" OR "1.2.7.90" OR "1.2.7.91" OR "1.2.7.92" OR "1.2.7.93" OR "1.2.7.94" OR "1.2.7.95" OR "1.2.7.96" OR "1.2.7.97" OR "1.2.7.98" OR "1.2.7.99" OR "1.2.7.100")"""
	baseSpecies = 'mycoplasma pneumoniae'
	results =  getSabioData(searchString, baseSpecies)
	print(len(results.entryList))
	
	for entry in results.entryList:
		print(entry.entryID)
		print(entry.km)
	
	"""
	array = [1,2,3]
	blue = []
	blue.append(array)
	print blue
	"""
