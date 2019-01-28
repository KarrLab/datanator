"""
Downloads and parse the latest release for 
ortholog -> uniprot mappings from omabrowser.org
"""

import pandas as pd
from time import time
from datetime import datetime
import sqlalchemy
from sqlalchemy import Column, Integer, Float, Date
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
import urllib.request

Base = sqlalchemy.ext.declarative.declarative_base()

def download_data(file_name, file_type, root_folder):
	print('Beginning file download with urllib: Downloading  ' + file_name)
	BASE_URL = "https://omabrowser.org/All/"
	#read in gzip files
	# for file in file_name:
	download_url = BASE_URL + file_name + file_type
	storage_location = root_folder + file_name + file_type
	response = urllib.request.urlretrieve(download_url, storage_location)

	return response

#load files into separate dataframes
def load_data(file_name, file_tpye, root_folder):
	# engine = create_engine('sqlite://', echo=False)
	# for file in file_name:
	pd_file = root_folder + file_name + file_type
	df = pd.read_csv(pd_file, compression='gzip', header = None, names = [x.upper() for x in file_name.split('-')], sep='\t', comment='#')
		# df.to_sql(file_name, con=engine)
	return df 

file_name = ["oma-uniprot", "oma-ncbi"]
file_type = ".txt.gz"
root_folder = 'cache/'

uniprot = download_data(file_name[0],file_type,root_folder)
oma_ncbi = download_data(file_name[1],file_type,root_folder)

uniprot_df = load_data(file_name[0], file_type, root_folder)
ncbi_df = load_data(file_name[1], file_type, root_folder)

print('merging dataframes...')
ncbi_uniprot = pd.merge(uniprot_df, ncbi_df, how='outer', on=file_name[0].split('-')[0].upper())

# toStore = NcbiUniprot(oma=ncbi_uniprot[0], uniprot=ncbi_uniprot['1_x'], ncbi=ncbi_uniprot['1_y'])
print('store merged dataframes in sql file...')
engine = sqlalchemy.create_engine('sqlite:///' + root_folder + 'ncbi_uniprot.sqlite')
# db_session = sqlalchemy.orm.sessionmaker(bind=engine)()
ncbi_uniprot.to_sql('ncbi_uniprot', con=engine)
