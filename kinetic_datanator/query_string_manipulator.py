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

