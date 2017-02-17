#For get parsed reaction, in order to make it universal, I can make it take inputs where
#the user can specify what the middle piece is (equals, arrow, etc), and he can also specity
#what to take out (like [c])

#to work on: this currently takes out the hydrogens. there is probably a better way 
#to do this



#this creates the searchString from the substrate 
def getQuerySearchString(subProdArray):
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


#takes in a search string and parses it into two arrays.
#one array for substrate IDs, and one array for product IDs


def getParsedReaction(reactionString):

	balancedMetab = []
	substrates = []
	products = []
	
	if "<==>" in reactionString:
		bothSides = reactionString.split("<==>")
	elif "==>" in reactionString:
		bothSides = reactionString.split("==>")
	elif "-->" in reactionString:
		bothSides = reactionString.split("-->")
	elif "<->" in reactionString:
		bothSides = reactionString.split("<->")
	elif "<=>" in reactionString:
		bothSides = reactionString.split("<=>")

	substrates = bothSides[0]
	products = bothSides[1]
	
	parsedSubstrates = []
	parsedProducts = []

	#IMPORTANT!!!!! we are getting rid of hydrogens here. 
	for entry in substrates.split(" "):
		if entry != "[c]:" and entry != "H" and entry != "H[c]" and entry != "h_m" and entry != "H[e]" and entry != "[m]:" and entry !=  "[e]:"and entry != "(2)" and entry !=  "+" and entry !=  "==>" and entry !=  "":
			#get rid of the [c] tag on some molecules
			if entry.find("[") != -1:
				entry = entry[:entry.find("[")]
			parsedSubstrates.append(entry)
	for entry in products.split(" "):
		if entry != "[c]:" and entry != "H" and entry != "H[c]" and entry != "H[e]" and entry != "[m]:" and entry !=  "[e]:"and entry != "(2)" and entry !=  "+" and entry !=  "==>" and entry !=  "":
			#get rid of the [c] tag on some molecules
			if entry.find("[") != -1:
				entry = entry[:entry.find("[")]
			parsedProducts.append(entry)


	balancedMetab.append(parsedSubstrates)
	balancedMetab.append(parsedProducts)
	return balancedMetab



if __name__ == '__main__':
	a = []#['dGMP', 'GMP', 'Lactose 6-phosphate', 'dGDP', "Orotidine 5'-phosphate", "Guanosine 3'-phosphate", "2',3'-Cyclic GMP", 'L-Arogenate', 'N-Acylneuraminate 9-phosphate', "Maltose 6'-phosphate", "5-Amino-6-(5'-phosphoribitylamino)uracil", '6-Phospho-beta-D-glucosyl-(1,4)-D-glucose', '2-Amino-4-hydroxy-6-(D-erythro-1,2,3-trihydroxypropyl)-7,8- dihydropteridine', '2-Amino-4-hydroxy-6-(erythro-1,2,3-trihydroxypropyl)dihydropteridine triphosphate', 'Dihydroneopterin phosphate', 'Ganciclovir', '8-Br-cGMP', "2'-Deoxyguanosine 3'-phosphate", "8-Azaguanosine-5'-monophosphate", '8-oxo-dGMP', 'Dihydroneopterin triphosphate', '8-oxo-dGTP', "2'-Deoxy-8-hydroxyguanosine"]
	b = ['dGMP', 'GMP', 'Lactose 6-phosphate', 'dGDP', "Orotidine 5'-phosphate", "Guanosine 3'-phosphate", "2',3'-Cyclic GMP", 'L-Arogenate', 'N-Acylneuraminate 9-phosphate', "Maltose 6'-phosphate", "5-Amino-6-(5'-phosphoribitylamino)uracil", '6-Phospho-beta-D-glucosyl-(1,4)-D-glucose', '2-Amino-4-hydroxy-6-(D-erythro-1,2,3-trihydroxypropyl)-7,8- dihydropteridine', '2-Amino-4-hydroxy-6-(erythro-1,2,3-trihydroxypropyl)dihydropteridine triphosphate', 'Dihydroneopterin phosphate', 'Ganciclovir', '8-Br-cGMP', "2'-Deoxyguanosine 3'-phosphate", "8-Azaguanosine-5'-monophosphate", '8-oxo-dGMP', 'Dihydroneopterin triphosphate', '8-oxo-dGTP', "2'-Deoxy-8-hydroxyguanosine"]
	c = []#["H2O"]
	d = ["phosphate"]
	stuff = []
	stuff.append([a, c])
	stuff.append([b, d])

	blank = [[], []]
	#print stuff


	string = 'h_m + 1a25dhvitd2_m + o2_m + nadph_m --> h2o_m + nadp_m + 1a2425thvitd2_m'
	print getParsedReaction(string)