from SabioInterface import getSabioData
import ReactionQueries
import numpy as np


class KineticInfo:
	def __init__(self, response, name = ""):
		self.name = name
		self.closestEntriesIDs = []
		self.medianEntry = []
		self.minEntry = []
		self.maxEntry = []
		self.liftInfo = None

		#should maybe change these two to only include the values of of the reactions
		#that are in the lowest species
		self.SabioReactionIDs = self.getFieldList(response.entryList, "reactionID")
		self.ECNumbers = self.getFieldList(response.entryList, "ECNumber")

		proximNums = self.getFieldList(response.entryList, "proximity")

		n = 0
		narrowedEntries = []
		while n < len(proximNums): 
			for entry in response.entryList:
				if entry.proximity == proximNums[n] and len(entry.__dict__[name])>0:
					narrowedEntries.append(entry)
			if len(narrowedEntries)>0:
				break
			else:
				n = n+1

		orderedEntries = sorted(narrowedEntries, key=lambda entry: float(entry.__dict__[name]))#, reverse=True)
		#records closes entry Ids - this is currently not working
		self.closestEntriesIDs = self.getFieldList(orderedEntries, "entryID")

		#to work on: if there are an even number, what should the median be?
		if len(orderedEntries) >= 2:
			self.minEntry = orderedEntries[0]
			self.maxEntry = orderedEntries[len(orderedEntries)-1]

		#to work on: make sure this media is done right
		if len(orderedEntries)>2 or len(orderedEntries)==1:
			number = float(len(orderedEntries))
			number = int(np.around(number/2+.1))
			self.medianEntry = orderedEntries[number-1]

		
	def getFieldList(self, entryList, field):
		values = []
		for entry in entryList:
			values.append(entry.__dict__[field])
		orderedValues = sorted(set(values))
		return orderedValues

class VmaxInfo(KineticInfo):
	def __init__(self, response):
		KineticInfo.__init__(self, response, "vmax")


class KmInfo(KineticInfo):
	def __init__(self, response):
		KineticInfo.__init__(self, response, "km")



class FormattedData:
	def __init__(self, id):
		self.id = id
		self.KmData = None
		self.VmaxData = None


def main():
	
	#this is the endgame
	formattedDataList = []

	filename='SmilesStuff2.xlsx'
	reactionQueries = ReactionQueries.generateReactionQueries(filename)

	for reactionQuery in reactionQueries:
		queryString = reactionQuery.getQueryString()

		answer =  getSabioData(queryString)

		fullResponse = FormattedData(reactionQuery.id)
		fullResponse.VmaxData = VmaxInfo(answer)
		fullResponse.KmData = KmInfo(answer)
		formattedDataList.append(fullResponse)

	#print fullResponse.__dict__
	for formattedData in formattedDataList:
		print formattedData.KmData.__dict__
		print formattedData.VmaxData.__dict__

 
if __name__ == '__main__':
	main()

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