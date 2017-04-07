""" 
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-06
:Copyright: 2017, Karr Lab
:License: MIT
"""

from . import ec_number_finder
from . import inchi_generator
from . import query_string_manipulator
import logging
import openpyxl


class Compound:
	def __init__(self, id, inchi_smiles = "", sabioNames = []):
		self.id = id
		self.inchi_smiles = inchi_smiles
		self.sabioNames = sabioNames

		#if len(self.sabioNames) == 0:
			#print(self.id)
			#print(self.inchi_smiles)
			#print(self.sabioNames)
			#print("\n")

class ReactionQuery:
	def __init__(self, id):
		self.id = id
		self.reaction_string = ""
		self.substrates = [] #this is a list of compound objects
		self.products = [] #this is a list of compound objects
		self.num_participants = [] #this is an array of two numbers, first number for substrates, second for products
		self.keggID = ""
		self.ec_number = ""
		self.generic_ec_number = ""

	def set_generic_ec_number_from_ezyme_algorithm(self):
		self.generic_ec_number = ec_number_finder.easy_find_ec_number(self.substrates, self.products)


def generateCompounds(wb):
	"""
	Args:
		wb (obj:`openpyxl.Workbook`): Takes an workbook that sheet with the name "Metabolites"
	
	Returns:
		:obj:`list` of `Compound`: list of compounds
	"""

	sabio_name_to_inchi_dict = inchi_generator.getSabioNameToInchiDict()
	ws = wb.get_sheet_by_name('Metabolites')
	
	#instantiate a compound object for each metabolite in the excel sheet
	compound_list = []
	for i in range(2, ws.max_row + 1):
		id = ws.cell(row=i, column=1).value

		smiles_or_inchi = ""
		structure = ws.cell(row=i, column=2).value
		sabio_names = []
		if structure != None:
			generic_inchi = inchi_generator.generateGenericInchi(structure)
			for name in sabio_name_to_inchi_dict:
				if sabio_name_to_inchi_dict[name] == generic_inchi:
					sabio_names.append(name)
			if len(sabio_names)==0:
				logging.info("reaction_queries: No Sabio Names found for {}".format(id))
			smiles_or_inchi = structure

		comp = Compound(id, smiles_or_inchi, sabio_names)
		compound_list.append(comp)

	return compound_list

def generate_reaction_queries(excelSheetObject):
	#Input an excel sheet data object (from openpyxl)
	#It outputs a list of reaction_query Objects

	logging.info("reaction_queries: Generating Reaction Queries")
	reaction_queries = [] #this is the only output of this method. Its a list of reaction_queries
	idToCompound = {} #this will be a dict used to find compound information about each metab id
	wb = excelSheetObject

	#this gets a dict that is used to correlate inchiString with Sabio's name for that compount
	sabioNameToInchiDict = inchi_generator.getSabioNameToInchiDict()

	#this instantiates a compound object for each metabolite in the excel sheet
	compoundList = generateCompounds(wb)
	#this creates a dictionary of ID to compound to make it quick to find the compound object 
	#with the compound ID. This works because each compound ID must be unique
	for comp in compoundList:
		idToCompound[comp.id] = comp


	#this creates a dictionary of reaction ID to reaction string
	#to work on: making it possible to only enter one row
	#to work on: make it possible not to have any reaction string data to begin with
	ws = wb.get_sheet_by_name('Reactions')
	idToReactionDict = {}
	for i in range(2, ws.max_row + 1):
		id = ws.cell(row=i, column=1).value
		reaction = ws.cell(row=i, column=2).value
		if reaction:
			idToReactionDict[id] = reaction
		else:
			idToReactionDict[id] = id


	#this is is where the ReactionQuery objects are instantiated
	#each reaction is parsed, each metabolite ID is converted to a Compound object
	#then the fields for each ReactionQuery object are set
	#finally the ReactionQuery object is added to the reaction_queries list
	for id in idToReactionDict:
		query = ReactionQuery(id)
		query.reaction_string = idToReactionDict[id]

		balancedMetab = query_string_manipulator.getParsedReaction(idToReactionDict[id])
		query.num_participants = [len(balancedMetab[0]), len(balancedMetab[1])]

		compoundBalancedMetab = []
		for side in balancedMetab:
			compoundSide = []
			for metab in side:
				try:
					compoundSide.append(idToCompound[metab])
				except:
					newUnrecognizedCompound = Compound(metab)
					compoundSide.append(newUnrecognizedCompound)
			compoundBalancedMetab.append(compoundSide)
		query.substrates = compoundBalancedMetab[0]
		query.products = compoundBalancedMetab[1]
		reaction_queries.append(query)

	return reaction_queries
