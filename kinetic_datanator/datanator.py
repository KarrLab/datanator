""" 
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-06
:Copyright: 2017, Karr Lab
:License: MIT
"""

from . import reaction_queries
from . import translator_for_sabio
from .create_excel_sheet import create_excel_sheet
from .sabio_interface import get_sabio_data
import logging
import numpy as np
import openpyxl
import os
#add other search criteria
#make sure to only use lift if smiles/inchi is present!!!!
#also, what should median be?
#make sure i fix ec_number in cases where we have (2)ADP or somethign like that


#finish debugging the EC finder

"""
.. autoclass:: Noodle
   :members: eat, slurp

   .. method:: boil(time=10)

      Boil the noodle *time* minutes
"""


class KineticInfo:
	def __init__(self, sabio_results, name = ""):

		self.name = name
		self.closest_entry_ids = []
		self.closest_entries = []
		self.closest_values = []
		self.median_entry = None
		self.min_entry = None
		self.max_entry = None
		self.lift_info = "Lift Not Used"
		self.reaction_list = sabio_results.get_field_list(sabio_results.entry_list, "reaction_stoichiometry")

		#should maybe change these two to only include the values of of the reactions
		#that are in the lowest species
		self.sabio_reaction_ids = sabio_results.get_field_list(sabio_results.entry_list, "reaction_id")
		self.ec_numbers = sabio_results.get_field_list(sabio_results.entry_list, "ec_number")


		proxim_nums = sabio_results.get_field_list(sabio_results.entry_list, "proximity")

		n = 0
		narrowed_entries = []
		while n < len(proxim_nums): 
			for entry in sabio_results.entry_list:
				if entry.proximity == proxim_nums[n] and len(entry.__dict__[name])>0:
					narrowed_entries.append(entry)
			if len(narrowed_entries)>0:
				break
			else:
				n = n+1
		self.closest_entries = narrowed_entries

		ordered_entries = sorted(narrowed_entries, key=lambda entry: float(entry.__dict__[name]))#, reverse=True)
		#records closest entry Ids - this is currently not working
		self.closest_entry_ids = sabio_results.get_field_list(ordered_entries, "entry_id")
		self.closest_values = sabio_results.get_field_list(ordered_entries, self.name)

		#to work on: if there are an even number, what should the median be?
		if len(ordered_entries) >= 2:
			self.min_entry = ordered_entries[0]
			self.max_entry = ordered_entries[len(ordered_entries)-1]

		#to work on: make sure this media is done right
		if len(ordered_entries)>2 or len(ordered_entries)==1:
			number = float(len(ordered_entries))
			number = int(np.around(number/2+.1))
			self.median_entry = ordered_entries[number-1]


	#this checks whether the current kinetic information is biologically useful. This is helpful to know
	#because if the data is not useful, then the use can try a more general search
	def has_relevant_data(self, proxim_limit = 1000):
		has_relevant_data = False
		if len(self.closest_entries)>0 and (self.closest_entries[0].proximity <= proxim_limit):
			has_relevant_data = True
		return has_relevant_data


class VmaxInfo(KineticInfo):
	def __init__(self, sabio_results):
		KineticInfo.__init__(self, sabio_results, "vmax")


class KmInfo(KineticInfo):
	def __init__(self, sabio_results):
		KineticInfo.__init__(self, sabio_results, "km")


class FormattedData:
	def __init__(self, id):#, sabio_results):
		self.id = id
		self.reaction_ids = []
		self.km_data = None
		self.vmax_data = None


def create_formatted_data(reaction_query, species, defaultValues, proxim_limit = 1000):
	#print("new")
	#for thing in reaction_query.substrates:
	#	print(thing.sabioNames)
	#print(reaction_query.id)
	#print(reaction_query.num_participants)
	searchString = translator_for_sabio.getSubstrateProductQueryString(reaction_query)
	logging.info("Datanator: Search String Found - {}".format(searchString))
	#print(searchString)

	if len(searchString)>0:
		searchString = defaultValues+searchString
	logging.info("Datanator: About to check Sabio with regular search string")
	sabio_results = get_sabio_data(searchString, species, reaction_query.num_participants)
	#instantiate a formatted data object with the reaction_id
	rxn = FormattedData(reaction_query.id)
	rxn.reaction_ids = sabio_results.get_field_list(sabio_results.entry_list, "reaction_id")

	lifted_sabio_results = None
	
	#get the km data. If data is not present, then it searches for similar reactions. 
	#proxim_limit is the taxonomic limit (in nodes) after which the data is presumed to be irrelevant 
	rxn.km_data = KmInfo(sabio_results)
	has_relevant_data = rxn.km_data.has_relevant_data(proxim_limit)
	logging.info("Datanator: Has relevant data for km - {}".format(has_relevant_data))
	if has_relevant_data == False:
		#we are about to set the generic_ec_number. However, first we need ot make sure
		#the user didn't provide their own. If the user provided their own, that generic_ec_number is used instead
		if len(reaction_query.generic_ec_number)==0:
			reaction_query.set_generic_ec_number_from_ezyme_algorithm()
		queryString = translator_for_sabio.getGenericECQueryString(reaction_query)
		if len(queryString)>0:
			queryString = defaultValues + queryString
		else:
			logging.error("Datanator: Tried to lift from EC, but did not work")
		lifted_sabio_results = get_sabio_data(queryString, species)
		rxn.km_data = KmInfo(lifted_sabio_results)
		if len(reaction_query.generic_ec_number)>0:
			rxn.km_data.lift_info = "Lifted From {}".format(reaction_query.generic_ec_number)
			logging.info("Datanator: Lifted From {} for Km".format(reaction_query.generic_ec_number))


	rxn.vmax_data = VmaxInfo(sabio_results)
	has_relevant_data = rxn.vmax_data.has_relevant_data(proxim_limit)
	logging.info("Datanator: Has relevant data for vmax - {}".format(has_relevant_data))
	if has_relevant_data == False:
		#check if liftessabio_results has already been searched. If it has, no need to search again
		if lifted_sabio_results == None:
			if len(reaction_query.generic_ec_number)==0:
				reaction_query.set_generic_ec_number_from_ezyme_algorithm()
			queryString = translator_for_sabio.getGenericECQueryString(reaction_query)
			if len(queryString)>0:
				queryString = defaultValues + queryString
			lifted_sabio_results = get_sabio_data(queryString, species)
		rxn.vmax_data = VmaxInfo(lifted_sabio_results)
		if len(reaction_query.generic_ec_number)>0:
			rxn.vmax_data.lift_info = "Lifted From {}".format(reaction_query.generic_ec_number)
			logging.info("Datanator: Lifted From {} for Vmax".format(reaction_query.generic_ec_number))

	return rxn

def get_kinetic_data(input_filename, output_filename, species, temp_range = [15, 40], enzyme_type = "wildtype", ph_range = [5,9], proxim_limit=1000):
	default_values = "enzymeType:{} AND TemperatureRange:[{} TO {}] AND pHValueRange:[{} TO {}] AND ".format(enzyme_type, temp_range[0], temp_range[1], ph_range[0], ph_range[1])
	
	rxns = []
	wb = openpyxl.load_workbook(filename=input_filename)
	rxn_qs = reaction_queries.generate_reaction_queries(wb)
	for rxn_q in rxn_qs:
		rxn = create_formatted_data(rxn_q, species, default_values, proxim_limit)
		rxns.append(rxn)
	create_excel_sheet(species, rxns, output_filename)
	return rxns

def get_kinetic_data_from_django(input_file, species, temp_range = [15, 40], enzyme_type = "wildtype", ph_range = [5,9], proxim_limit=1000):
	default_values = "enzymeType:{} AND TemperatureRange:[{} TO {}] AND pHValueRange:[{} TO {}] AND ".format(enzyme_type, temp_range[0], temp_range[1], ph_range[0], ph_range[1])
	
	rxns = []
	rxn_qs = reaction_queries.generate_reaction_queries(input_file)
	for rxn_q in rxn_qs:
		rxn = create_formatted_data(rxn_q, species, default_values, proxim_limit)
		rxns.append(rxn)
	#create_excel_sheet(species, rxns, outputFilename)
	return rxns
