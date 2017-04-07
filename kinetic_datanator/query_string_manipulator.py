""" 
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-06
:Copyright: 2017, Karr Lab
:License: MIT
"""

#For get parsed reaction, in order to make it universal, I can make it take inputs where
#the user can specify what the middle piece is (equals, arrow, etc), and he can also specity
#what to take out (like [c])

#to work on: this currently takes out the hydrogens. there is probably a better way 
#to do this



def getQuerySearchString(subProdArray):
	#this creates the searchString from the substrate 

	substrates = []
	for array in subProdArray[0]:
		if len(array)>0:
			substrates.append(array)

	products = []
	for array in subProdArray[1]:
		if len(array)>0:
			products.append(array)


	searchString = ""

	if len(substrates)>0:
		searchString = searchString + "("
		for metabArray in substrates:
			subString = "("
			for entry in metabArray:
				subString = subString + "Substrate:" + '"' + entry + '"' + " OR "
			subString = subString[:-4] + ")"
			searchString = searchString + subString + " AND "
		searchString = searchString[:-5] + ")"

	if len(substrates)>0 and len(products)>0:
		searchString = searchString + " AND "
	
	if len(products)>0:
		searchString = searchString + "("
		for metabArray in products:
			subString = "("
			for entry in metabArray:
				subString = subString + "Product:" + '"' + entry + '"' + " OR "
			subString = subString[:-4] + ")"
			searchString = searchString + subString + " AND "
		searchString = searchString[:-5] + ")"
	


	return searchString


def getParsedReaction(reaction_string):
	#takes in a search string and parses it into two arrays.
	#one array for substrate IDs, and one array for product IDs

	balancedMetab = []
	substrates = []
	products = []
	
	if "<==>" in reaction_string:
		bothSides = reaction_string.split("<==>")
	elif "==>" in reaction_string:
		bothSides = reaction_string.split("==>")
	elif "-->" in reaction_string:
		bothSides = reaction_string.split("-->")
	elif "<->" in reaction_string:
		bothSides = reaction_string.split("<->")
	elif "<=>" in reaction_string:
		bothSides = reaction_string.split("<=>")
	elif "=" in reaction_string:
		bothSides = reaction_string.split("=")

	substrates = bothSides[0]
	products = bothSides[1]
	
	parsedSubstrates = []
	parsedProducts = []

	#IMPORTANT!!!!! we are getting rid of hydrogens here. 
	for entry in substrates.split(" "):
		if not (is_number(entry)) and entry != "[c]:" and entry != "H" and entry != "H[c]" and entry != "h_m" and entry != "h_x" and entry != "h_c" and entry != "H[e]" and entry != "[m]:" and entry !=  "[e]:"and entry != "(2)" and entry !=  "+" and entry !=  "==>" and entry !=  "":
			#get rid of the [c] tag on some molecules
			if entry.find("[") != -1:
				entry = entry[:entry.find("[")]
			parsedSubstrates.append(entry)
	for entry in products.split(" "):
		if not (is_number(entry)) and entry != "[c]:" and entry != "H" and entry != "H[c]" and entry != "h_m" and entry != "h_x" and entry != "h_c" and entry != "H[e]" and entry != "[m]:" and entry !=  "[e]:"and entry != "(2)" and entry !=  "+" and entry !=  "==>" and entry !=  "":
			#get rid of the [c] tag on some molecules
			if entry.find("[") != -1:
				entry = entry[:entry.find("[")]
			parsedProducts.append(entry)


	balancedMetab.append(parsedSubstrates)
	balancedMetab.append(parsedProducts)

	return balancedMetab



def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
