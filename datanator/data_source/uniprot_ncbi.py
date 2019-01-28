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
from sqlalchemy.orm import sessionmaker
import urllib.request

Base = sqlalchemy.ext.declarative.declarative_base()

class Uniprot(Base):
    """ Represents the oma-uniprot association table

    Attributes:
        oma (:obj:`str`): oma id
        uniprot (:obj:`list`): uniprot id
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    oma = sqlalchemy.Column(sqlalchemy.String(), unique=True, index=True)

    uniprot = sqlalchemy.Column(sqlalchemy.String(), unique=True, index=True)

    __tablename__ = 'uniprot'

class Ncbi(Base):
    """ Represents the oma-ncbi association table

    Attributes:
        oma (:obj:`str`): oma id
        ncbi (:obj:`list'): ncbi id
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    oma = sqlalchemy.Column(sqlalchemy.String(), unique=True, index=True)

    ncbi = sqlalchemy.Column(sqlalchemy.String(), unique=True, index=True)

    __tablename__ = 'ncbi'


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
	df = pd.read_csv(pd_file, compression='gzip', header = 0, sep='\t', comment='#')
		# df.to_sql(file_name, con=engine)
	return df 

file_name = ["oma-uniprot", "oma-ncbi"]
file_type = ".txt.gz"
root_folder = './cache/'

uniprot = download_data(file_name[0],file_type,root_folder)
oma_gi = download_data(file_name[1],file_type,root_folder)

uniprot_df = load_data(file_name[0], file_type, root_folder)
oma_ncbi_df = load_data(file_name[1], file_type, root_folder)

