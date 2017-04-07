""" 
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-06
:Copyright: 2017, Karr Lab
:License: MIT
"""

from .taxon_finder import get_taxonomic_distance
from requests.packages.urllib3.exceptions import InsecureRequestWarning
import logging
import requests 
requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

#to do later. Turn the data structures into a generic kinetic information gathering framework. 
#or perhaps I can use this for everything. Meaning, every database is going to be built on entries
#and upon totals. So what I could do is that I can have generic "database-dataStructures". This will be 
#Entries, and then total results which contain the results. 

ENTRY_ID_QUERY_URL = 'http://sabiork.h-its.org/sabioRestWebServices/searchKineticLaws/entryIDs'
PARAM_QUERY_URL = 'http://sabiork.h-its.org/entry/exportToExcelCustomizable' 

class Entry:

	def __init__(self, text_list):
		self.entry_id = ""
		self.vmax = ""
		self.km = ""
		self.km_substrate = ""
		self.species = ""
		self.ec_number = ""
		self.reaction_id = ""
		self.proximity = ""
		self.reaction_stoichiometry = ""

		line = text_list[0].split("\t")

		self.species = line[0]
		self.reaction_id = line[1]
		self.reaction_stoichiometry = line[2]
		self.num_participants = self.parse_sabio_reaction_string(self.reaction_stoichiometry)
		self.ec_number = line[3]
		self.entry_id = line[4]

		kinetic_info = []
		for line in text_list:
			kinetic_info.append(line[line.find(self.entry_id)+6:])

		for line in kinetic_info:
			if line[0:4] == "Vmax":
				parsed = line.split("\t")
				#print(parsed)
				#print(self.entry_id)
				#if parsed[5] == 'mol*s^(-1)*g^(-1)' or parsed[5] == 'mol*g^(-1)*s^(-1)':
				self.vmax = parsed[2]
			if line[0:2]=="Km":
				parsed = line.split("\t")
				#print(parsed)
				#print(self.entry_id)
				#if parsed[5] == "M":
				self.km = parsed[2]
				self.km_substrate = parsed[1]

	def parse_sabio_reaction_string(self, reaction_string):

		right, left = reaction_string.split(" = ")
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
	def __init__(self, sabio_file):
		self.entry_list = []

		if len(sabio_file)>0:
			split_sabio_file = sabio_file.split("\n")
			entry_text_arrays = []
			i = split_sabio_file[1].split("\t")[4]
			text_array = []

			#this code ensures that all the lines that relate to a single entry_id get passed
			for line in split_sabio_file:
				if line[:8] != "Organism" and len(line)>0:
					if i == line.split("\t")[4]:
						text_array.append(line)
					else:
						entry_text_arrays.append(text_array)
						text_array = []
						text_array.append(line)
						i = line.split("\t")[4]
			entry_text_arrays.append(text_array)

			for entry_data in entry_text_arrays:
				#print(entry_data)
				#print("")
				entry = Entry(entry_data)
				self.entry_list.append(entry)



	def get_field_list(self, entry_list, field):
		values = []
		for entry in entry_list:
			values.append(entry.__dict__[field])
		ordered_values = sorted(set(values))
		return ordered_values

	def narrow_by_nums(self, num_participants):
		i = 0
		for entry in self.entry_list:
			if entry.num_participants != num_participants:
				self.entry_list.remove(entry)


			#get the proximity for each entry
			#for entry in self.entry_list:
			#	entry.proximity = get_taxonomic_distance('mycoplasma pneumoniae', entry.species)


def get_sabio_data(query, base_species, num_participants = []):
	#takes in a search dictionary or search string. Returns a TotalResult object if something is found. 
	#If nothing is found, it returns "No results found for query"

	logging.info("Sabio Interface: Looking for query in Sabio")
	if not query:
		logging.info("Sabio Interface: Sabio failed to respond - no search string entered")
		return TotalResult("")	

	# ask SABIO-RK for all entry_ids matching a query 
	if isinstance(query, dict):	
		query_string = ' AND '.join(['{}:{}'.format(k, v) for k, v in query.items()])
	if isinstance(query, str):
		query_string = query

	request = requests.get(ENTRY_ID_QUERY_URL, params={'format': 'txt', 'q': query_string}) 
	#print(request)
	request.raise_for_status() # raise if 404 error 

	if request.text == "No results found for query":
		logging.info("Sabio Interface: Sabio failed to respond - No results found for query")
		return TotalResult("")
	# each entry is reported on a new line

	entry_ids = [int(x) for x in request.text.strip().split('\n')]
	print("{} matching entries found".format(len(entry_ids)))
	logging.info("Sabio Interface: Sabio found {} results for query".format(len(entry_ids)))

	# encode next request, for parameter data given entry IDs 
	request = requests.post(PARAM_QUERY_URL, 
		params={
			'format':'tsv', 
			'fields[]': ['Organism', 'SabioReactionID', 'Reaction','ECNumber', 'EntryID', 'Parameter'],
		},
		data={
			'entryIDs[]': entry_ids
		})
	request.raise_for_status()
	text_file = request.text
	#print(text_file)
	result = TotalResult(text_file)

	for entry in result.entry_list:
		entry.proximity = get_taxonomic_distance(base_species, entry.species)

	if len(num_participants)>0:
		result.narrow_by_nums(num_participants)


	return result
