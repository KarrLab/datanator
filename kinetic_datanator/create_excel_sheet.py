""" 
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-06
:Copyright: 2017, Karr Lab
:License: MIT
"""

from . import taxon_finder
from openpyxl import Workbook
import os

def create_excel_sheet(species, formatted_data_list, output_filename):
	""" Takes in FormattedData as an input. Saves a file 

	Args:
		species
		formatted_data_list
		output_filename
	"""
	
	# create a workbook
	wb = Workbook()
	
	# write kinetic data to a worksheet
	ws = wb.active
	ws.title = "Kinetics"

	headers = [
		'Name', 'Sabio Reaction IDs', 
		'Vmax Median Value', 'Vmax Median Entry', 'Vmax Min Entry', 'Vmax Max Entry', 'Vmax Proximity', 'Vmax Lift Info' 'Closest Vmax', 'Sabio Entry IDs', 'Closest Vmax Values',
		'Km Median Value', 'Km Median Entry', 'Km Min Enry', 'Km Max Entry', 'Km Proximity', 'Km Lift Info', 'Closest Km Sabio Entry IDs', 'Closest Km Values'
		]
	for i_col, value in enumerate(headers):
		ws.cell(row=1, column = i_col + 1).value = value

	for i_row, formatted_data in enumerate(formatted_data_list):
		ws.cell(row=i_row + 2, column=1).value = formatted_data.id
		ws.cell(row=i_row + 2, column=2).value = toCSV(formatted_data.reaction_ids)

		# vmax information
		ws.cell(row=i_row + 2, column=3).value = formatValue(formatted_data.vmax_data.median_entry, "vmax")
		ws.cell(row=i_row + 2, column=4).value = formatValue(formatted_data.vmax_data.median_entry)
		ws.cell(row=i_row + 2, column=5).value = formatValue(formatted_data.vmax_data.min_entry)
		ws.cell(row=i_row + 2, column=6).value = formatValue(formatted_data.vmax_data.max_entry)
		ws.cell(row=i_row + 2, column=7).value = formatProxim(formatted_data.vmax_data.closest_entries)
		ws.cell(row=i_row + 2, column=8).value = formatValue(formatted_data.vmax_data, "lift_info")
		ws.cell(row=i_row + 2, column=9).value = toCSV(formatted_data.vmax_data.closest_entry_ids)
		ws.cell(row=i_row + 2, column=10).value = toCSV(formatted_data.vmax_data.closest_values)

		# km information
		ws.cell(row=i_row + 2, column=11).value = formatValue(formatted_data.km_data.median_entry, "km")
		ws.cell(row=i_row + 2, column=12).value = formatValue(formatted_data.km_data.median_entry)
		ws.cell(row=i_row + 2, column=13).value = formatValue(formatted_data.km_data.min_entry)
		ws.cell(row=i_row + 2, column=14).value = formatValue(formatted_data.km_data.max_entry)
		ws.cell(row=i_row + 2, column=15).value = formatProxim(formatted_data.km_data.closest_entries)
		ws.cell(row=i_row + 2, column=16).value = formatValue(formatted_data.km_data, "lift_info")
		ws.cell(row=i_row + 2, column=17).value = toCSV(formatted_data.km_data.closest_entry_ids)
		ws.cell(row=i_row + 2, column=19).value = toCSV(formatted_data.km_data.closest_values)

	# write taxonomic information to a sheet
	ws = wb.create_sheet(title="Taxonomic information")
	ws.cell(row=1, column=1, value='Taxonomy')
	ws.cell(row=1, column=2, value='Proximity')
	for i_row, value in enumerate(taxon_finder.get_taxonomic_lineage(species)):
		ws.cell(row=i_row+2, column=1, value=value)
		ws.cell(row=i_row+2, column=2, value=i_row + 1)

	# save workbook to a file
	wb.save(output_filename)

def formatValue(outerField, innerField = ""):
	value = ""
	if outerField != None:
		if innerField:
			value =  "{}".format(outerField.__dict__[innerField])
		else:
			value = "{}".format(outerField.__dict__)
	if outerField == None:
		value =  "No Data"
	return value

def formatProxim(closest_entries):
	if len(closest_entries)>0:
		return "{}".format(closest_entries[0].proximity)
	else:
		return "No Data"

def toCSV(array):

	csv = ""
	for item in array:
		csv = csv + "{}".format(item) + ", "

	if len(csv)>0:
		csv = csv[:-2]
	return csv
