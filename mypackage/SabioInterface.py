from TaxonFinder import getTaxonomicDistance
import requests 

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

	if len(query_dict)==0:
		return TotalResult("") 

	ENTRYID_QUERY_URL = 'http://sabiork.h-its.org/sabioRestWebServices/searchKineticLaws/entryIDs' 
	#ENTRYID_QUERY_URL = 'http://sabiork.h-its.org/sabioRestWebServices/searchKineticLaws/sbml'
	PARAM_QUERY_URL = 'http://sabiork.h-its.org/entry/exportToExcelCustomizable' 

	# ask SABIO-RK for all EntryIDs matching a query 
	if isinstance(query_dict, dict):	
		query_string = ' AND '.join(['%s:%s' % (k,v) for k,v in query_dict.items()]) 
	if isinstance(query_dict, str):
		query_string = query_dict
	query = {'format':'txt', 'q':query_string} 


	request = requests.get(ENTRYID_QUERY_URL, params = query) 
	#print request
	request.raise_for_status() # raise if 404 error 

	if request.text == "No results found for query":
		return TotalResult("")
	# each entry is reported on a new line

	entryIDs = [int(x) for x in request.text.strip().split('\n')]
	print('%d matching entries found.' % len(entryIDs)) 

	# encode next request, for parameter data given entry IDs 
	data_field = {'entryIDs[]': entryIDs} 
	query = {'format':'tsv', 'fields[]':['Organism',"SabioReactionID", 'Reaction','ECNumber','EntryID',  'Parameter']}


	request = requests.post(PARAM_QUERY_URL, params=query, data=data_field) 
	request.raise_for_status()
	textFile = request.text
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
	answer =  getSabioData(query_dict)

	print answer
	#blue = TotalResult(answer)
	
	for entry in answer.entryList:
		print entry
	

