import requests
import os
import sys
import datetime
import json

reload(sys)  
sys.setdefaultencoding('utf8')

class DownloadExperiments():

	def download_single_year(self, year):
		"""
		Gets a JSON of all queries for a single year. Saves it into a directory called
		"AllSamples". Creates this directory if it doesn't exist.

		Args:
			year (obj:'int':) the year that is being collected
		"""
		directory = os.path.join('.', 'AllExperiments')
		if not os.path.exists(directory):
			os.makedirs(directory)
		file = open('./AllExperiments/{}.txt'.format(year), 'w')

		file.write(requests.get("https://www.ebi.ac.uk/arrayexpress/json/v3/experiments?date=[{}-01-01+{}-12-31]"
			.format(year, year)).text)

	def download_all_experiments(self, download_range=datetime.datetime.now().year+1):
		"""
		Downloads all experiments by iterating through the years and calling "donload_single_year"
		"""
		for year in range(2001, download_range):
			self.download_single_year(year)
		

class DownloadSamples():

	def download_single_sample(self, ax_num):
		"""
		Gets a JSON of all samples in a single experiment. Saves it into a directory called
		"AllSamples". Creates this directory if it doesn't exist.

		Args:
			ax_num (obj:'str':) accession number of the experiment
		"""

		directory = os.path.join('.', 'AllSamples')
		if not os.path.exists(directory):
			os.makedirs(directory)
		file = open('./AllSamples/{}.txt'.format(ax_num), 'w')
		file.write(requests.get("https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}/samples".format(ax_num)).text)


def download_all_metadata(download_range=datetime.datetime.now().year+1):
	"""
	Downloads all medatata from array exrpess on their samples and experiments. The metadata
	is saved as the text file. Within the text files, the data is stored as a json object. 
	"""
	DownloadExperiments().download_all_experiments(download_range)

	all_ax_nums = []
	for year in range(2001, download_range):
		metadata = json.loads(open("./AllExperiments/{}.txt".format(year), 'r').read().encode('utf8'))
		for entry in metadata['experiments']['experiment']:
			all_ax_nums.append(entry['accession'])
	for num in all_ax_nums:
		DownloadSamples().download_single_sample(num)

