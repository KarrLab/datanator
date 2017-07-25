import requests
import os
import sys
import datetime
import demjson

#reload(sys)
#sys.setdefaultencoding('utf8')

ENDPOINT = 'https://www.ebi.ac.uk/arrayexpress/json/v3/experiments'

class DownloadExperiments(object):	

	def download_single_year(self, year):
		"""
		Gets a JSON of all queries for a single year. Saves it into a directory called
		"AllSamples". Creates this directory if it doesn't exist.

		Args:
			year (obj:'int':) the year that is being collected
		"""
		directory = 'AllExperiments'
		if not os.path.exists(directory):
			os.makedirs(directory)
		response = requests.get(ENDPOINT + '?date=[{}-01-01+{}-12-31]'.format(year, year))
		response.raise_for_status()
		with open(os.path.join(directory, '{}.txt'.format(year)), 'wb') as file:
			file.write(response.content)

	def download_all_experiments(self, start_year, end_year):
		"""
		Downloads all experiments by iterating through the years and calling "donload_single_year"
		"""
		for year in range(start_year, end_year+1):
			print(year)
			self.download_single_year(year)
		

class DownloadSamples():

	def download_single_sample(self, ax_num):
		"""
		Gets a JSON of all samples in a single experiment. Saves it into a directory called
		"AllSamples". Creates this directory if it doesn't exist.

		Args:
			ax_num (obj:'str':) accession number of the experiment
		"""

		directory = 'AllSamples'
		if not os.path.exists(directory):
			os.makedirs(directory)
		response = requests.get(ENDPOINT + "/{}/samples".format(ax_num))
		response.raise_for_status()
		with open(os.path.join(directory, '{}.txt'.format(ax_num)), 'wb') as file:
			file.write(response.content)


def download_all_metadata(start_year=2001, end_year=datetime.datetime.now().year):
	"""
	Downloads all medatata from array exrpess on their samples and experiments. The metadata
	is saved as the text file. Within the text files, the data is stored as a json object. 
	"""
	DownloadExperiments().download_all_experiments(2001, end_year)

	all_ax_nums = []
	for year in range(start_year, end_year+1):
		metadata = demjson.decode_file(os.path.join('AllExperiments', '{}.txt'.format(year)), encoding='utf-8')
		for entry in metadata['experiments']['experiment']:
			all_ax_nums.append(entry['accession'])
	for num in all_ax_nums:
		DownloadSamples().download_single_sample(num)

