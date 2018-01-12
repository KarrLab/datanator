#import numpy as np
#from ete3 import NCBITaxa
#import requests, sys
from six.moves.urllib.request import urlretrieve
#import ftplib
#from ftplib import FTP
#import requests
import os
import shutil


def run(sample, top_dir):
	#build kallisto index file
	DIRNAME = "{}/CDNA_FILES".format(top_dir)
	if not os.path.isdir(DIRNAME):
		os.makedirs(DIRNAME)
	spec_name = sample.ensembl_info[0].organism_strain
	file_name = "{}/{}.cdna.all.fa.gz".format(DIRNAME, spec_name)
	url = sample.ensembl_info[0].url
	if not os.path.isfile(file_name):
		file = urlretrieve(url, '{}/{}.cdna.all.fa.gz'.format(top_dir, spec_name))
		shutil.move('{}/{}.cdna.all.fa.gz'.format(top_dir, spec_name), DIRNAME)
	cur_dir = "{}/explore".format(os.path.dirname(os.path.abspath(__file__)))
	os.chdir(top_dir)
	KALLISTO_DIR = "{}/kallisto_index_files".format(top_dir)
	if not os.path.isdir(KALLISTO_DIR):
		os.makedirs(KALLISTO_DIR)
	if not os.path.isfile("{}/{}.idx".format(KALLISTO_DIR, spec_name)):
		os.system("kallisto index -i {}.idx {}".format(spec_name, file_name))
		shutil.move("{}/{}.idx".format(top_dir, spec_name), KALLISTO_DIR)
		
# I might as well create the special kallisto files here too. 

