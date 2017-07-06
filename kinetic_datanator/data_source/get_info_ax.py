import requests
import json
import datetime
import sys  

reload(sys)  
sys.setdefaultencoding('utf8')


def get_text_file(a_list, title):
	file = open("{}.txt".format(title), "w")
	file.write("")
	file.close()
	file = open("{}.txt".format(title), "a")
	for thing in a_list:
		if type(thing) is list:
			file.write(str(thing))
		else:
			file.write(thing)
			file.write('\n')



unique_keys = []
unique_variables = []
unique_exp_types = []
unique_characteristics = []
list_accession_nums = []



for year in range(2001,2018):
	#open("AX_Json_Objects/{}.txt".format(year), 'r').read().encode('utf8')
	metadata = json.loads(open("AX_Json_Objects/{}.txt".format(year), 'r').read().encode('utf8'))
	all_keys = []
	for entry in metadata['experiments']['experiment']:
		for key in entry:
			all_keys.append(key)
	unique_keys = list(set(unique_keys + (all_keys)))



	all_variables = []
	for entry in metadata['experiments']['experiment']:
		if 'experimentalvariable' in entry:
			for var in entry['experimentalvariable']:
				all_variables.append(str(var['name']))

	unique_variables = list(set(unique_variables + (all_variables)))

	all_exp_types = []
	for entry in metadata['experiments']['experiment']:
		if 'experimenttype' in entry:
			all_exp_types.append(str(entry['experimenttype'][0]))
	unique_exp_types = list(set(unique_exp_types + (all_exp_types)))



	all_characteristics = []
	for entry in metadata['experiments']['experiment']:
		if 'samplecharacteristic' in entry:
			for char in entry['samplecharacteristic']:
				all_characteristics.append(str(char['category']))
	unique_characteristics = list(set(unique_characteristics + (all_characteristics)))


	#list_accession_nums = []
	for entry in metadata['experiments']['experiment']:
		list_accession_nums.append(entry['accession'])
		if year == 2017:
			print entry['accession']


get_text_file(unique_keys, "All Experiment Fields")
get_text_file(unique_variables, "All Experimental Variables")
get_text_file(unique_exp_types, "All Experiment Types")
get_text_file(unique_characteristics, "All Sample Charateristics")
get_text_file(list_accession_nums, "All Accession Numbers")
