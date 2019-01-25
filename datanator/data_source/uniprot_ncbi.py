"""
Downloads and parse the latest release for 
ortholog -> uniprot mappings from omabrowser.org
"""

import panda as pd
from time import time
from datetime import datetime
from sqlalchemy import Column, Integer, Float, Date
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import urllib.request

file_name = ["oma-uniprot", "oma-gi"]
file_type = ".txt.gz"
root_folder = './cache/'

def download_data(file_name, file_type, root_folder):
	print('Beginning file download with urllib...')
	BASE_URL = "https://omabrowser.org/All/"
	#read in gzip files
	for file in file_name:
		download_url = BASE_URL + file + file_type
		storage_location = root_folder + file + file_type
		response = urllib.request.urlretrieve(download_url, storage_location)

def load_data(file_name, file_tpye, root_folder):
	for file in file_name:
		pd_file = root_folder + file + file_type
		df = pd.read_csv(pd_file, compression='gzip', header = 0, sep='\t', comment='#')


