
def getNameToInchiDict():
	inchiCompoundTranslator = {}
	with open("InchiToCompound.txt") as f:
		for line in f:
			values = line.split(" - ")
			#print values
			if values[1] != "No Inchi Found":
				#inchiCompoundTranslator[values[1][9:]] = values[2]#[9:]
				inchiCompoundTranslator[values[2][:-1]] = values[1][9:]
				#print values[1][9:]
			else:
				print "blue"
	return inchiCompoundTranslator
			