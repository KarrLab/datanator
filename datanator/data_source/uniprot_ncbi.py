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
def load_data(file_name, file_type, root_folder):
	# engine = create_engine('sqlite://', echo=False)
	# for file in file_name:
	pd_file = root_folder + file_name + file_type
	df = pd.read_csv(pd_file, compression='gzip', header = None, names = [x.upper() for x in file_name.split('-')], sep='\t', comment='#')
		# df.to_sql(file_name, con=engine)
	return df 

def sql_data(df1,df2,file_name,root_folder,db_type):

	print('merging dataframes...')
	#1 uniprot_df 2 ncbi_df
	merged_df = pd.merge(df1, df2, how='outer', on=file_name[0].split('-')[0].upper())
	name1 = file_name[0].split('-')[1]
	name2 = file_name[1].split('-')[1]
	db_name = name1 + '_' + name2
	print('store merged dataframes in ' + db_type + ' file...')
	engine = sqlalchemy.create_engine('sqlite:///' + root_folder + db_name + db_type)
	merged_df.to_sql(db_name, con=engine)

if __name__ == '__main__':

	file_name = ["oma-uniprot", "oma-ncbi"]
	file_type = ".txt.gz"
	root_folder = 'cache/'
	db_type = '.sqlite'

	uniprot = download_data(file_name[0],file_type,root_folder)
	oma_ncbi = download_data(file_name[1],file_type,root_folder)

	uniprot_df = load_data(file_name[0], file_type, root_folder)
	ncbi_df = load_data(file_name[1], file_type, root_folder)

	sql_data(uniprot_df, ncbi_df,file_name,root_folder,db_type)