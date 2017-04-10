import requests
import os

def trimInchi(inchi):
	if "/h" in inchi:
		end = inchi.index("/h")
		inchi = inchi[:end]
	return inchi



file = open("InchiToCompoundPubChemSpringBoard2.txt", "w")
file.write("")

errorFile = open("datasetErrorsSpring.txt", "w")
errorFile.write("")

with open(os.path.join(os.path.join(os.path.dirname(os.path.realpath(__file__)), "InchiToCompoundPubChemSpringBoard.txt"))) as f:
	i = 0
	for line in f:

		newLine = line

		try:

			values = line.split(" - ")
			foundInchi = ""
			if values[1] == "No Inchi Found":
				compName = values[2][:-1]
				response = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/InChI/TXT".format(compName)).text[:-1]
				#print(request.text)
				if response[0:5].lower() == "inchi":
					#print(compName)
					#print(response)
					if response.count("InChI") == 1:
						foundInchi = response


					elif response.count("InChI") > 1:
						inchis = response.split("\n")
						allSame = True
						for inchi in inchis:
							if trimInchi(inchis[0])!=trimInchi(inchi):
								allSame = False


						if allSame==True:
							foundInchi = trimInchi(inchis[0])

			if len(foundInchi)>0:
				newLine = values[0] + " - " + """Inchi_0">""" + foundInchi + " - " + values[2]
				i = i+1
				print(i)
		except:
			errorFile.write(line)

		file.write(newLine)
	#return inchiCompoundTranslator

