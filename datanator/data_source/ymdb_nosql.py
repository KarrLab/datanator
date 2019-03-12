import io
import os
import json
import requests
import requests.exceptions
import warnings
import zipfile
import xmltodict

class YmdbNoSQL():
	def __init__(self, verbose):
		self.verbose = verbose

	domain = 'http://ymdb.ca'
	compound_index = domain + '/system/downloads/current/ymdb.json.zip' # list of metabolites
	compound_url = domain + '/compounds/{}.xml'
	collection_dir = './cache/ymdb/'

	# each compound's collection of imformation is written into a json file for mongodb
	def write_to_json(self):

		if self.verbose:
			print('Download list of all compounds: ...')
		
		response = requests.get(self.compound_index)
		response.raise_for_status()

		if self.verbose:
			print('... Done!')
		if self.verbose:
			print('Unzipping and parsing compound list ...')

		with zipfile.ZipFile(io.BytesIO(response.content), 'r') as zip_file:
			with zip_file.open('ymdb.json', 'r') as json_file:
				entries = json.load(json_file)

		if self.verbose:
			print ('  found {} compounds'.format(len(entries)))

		entries.sort(key=lambda e: e['ymdb_id'])


		if self.verbose:
			print('Downloading {} compounds ...'.format(len(entries)))

		for i_entry, entry in enumerate(entries):

			if self.verbose and (i_entry % 10 == 0):
				print('  Downloading compound {} of {}'.format(i_entry + 1, len(entries)))

			# actual xml file for each metabolite
			response = requests.get(self.compound_url.format(entry['ymdb_id']))
			try:
				response.raise_for_status()
			except requests.exceptions.HTTPError:
				warnings.warn('Unable to download data for compound {}'.format(entry['ymdb_id']))
				continue

			doc = xmltodict.parse(response.text)

			new_doc = doc['compound']  # delete key "compound" but keep key's value
			
			os.makedirs(os.path.dirname(self.collection_dir), exist_ok=True)
			file_name = os.path.join(self.collection_dir + entry['ymdb_id'] + '.json')
			with open(file_name, "w") as f:
				f.write(json.dumps(new_doc))

		

if __name__ == '__main__':
	r = YmdbNoSQL(True)
	r.write_to_json()



