from . import download_cDNA
import requests
from six.moves.urllib.request import urlretrieve
import os
import pandas as pd
import numpy as np


class Sample():
	def __init__(self, name, organism, url):
		self.name = name
		self.organism = organism
		self.url = url
class Experiment():
	def __init__(self, name, samples):
		self.name = name
		self.samples = samples


def get_processed_data(experiment, top_dir):
	TOP_DIR=top_dir
	EXP_DIRNAME = "{}/{}".format(TOP_DIR, experiment.id)
	if not os.path.isdir(EXP_DIRNAME):
		os.makedirs(EXP_DIRNAME)

	for sample in experiment.samples:
		#print(organism)
		download_cDNA.run(sample, top_dir)
		spec_name = sample.ensembl_info[0].organism_strain
		SAMP_DIR = "{}/{}".format(EXP_DIRNAME, sample.name)
		if not os.path.isdir(SAMP_DIR):
			os.makedirs(SAMP_DIR)
		sample_url = sample.fastq_urls[0].url
		print(sample_url)
		fastq_files = """"""
		for num, url in enumerate(sample.fastq_urls):
			file = urlretrieve(url.url, '{}/{}_{}.fastq.gz'.format(SAMP_DIR, sample.name, num))#there used to be a spacae after "gz". I removed it
			fastq_files = fastq_files + """"{}/{}_{}.fastq.gz" """.format(SAMP_DIR, sample.name, num)

		print("read type: " + sample.experiment.read_type)
		if sample.experiment.read_type == "single":
			print("single muffins")
			os.system("""kallisto quant -i "{}/kallisto_index_files/{}.idx" -o "{}/output" --single -l 180 -s 20 {} """.format(TOP_DIR, spec_name, SAMP_DIR, fastq_files))
		elif sample.experiment.read_type == "paired":
			print(fastq_files)
			print("""kallisto quant -i "{}/kallisto_index_files/{}.idx" -o "{}/output" {} """.format(TOP_DIR, spec_name, SAMP_DIR, fastq_files))
			os.system("""kallisto quant -i "{}/kallisto_index_files/{}.idx" -o "{}/output" {} """.format(TOP_DIR, spec_name, SAMP_DIR, fastq_files))
	
		new_pandas = pd.read_csv('{}/output/abundance.tsv'.format(SAMP_DIR), sep='\t').set_index("target_id")
		total_tpm = 0
		for num in new_pandas['tpm']: 
			total_tpm = total_tpm + num
		new_column = [num/total_tpm for num in new_pandas['tpm']]
		new_pandas['percent total'] = new_column
		new_pandas.to_pickle("{}/{}_abundances_binary".format(SAMP_DIR, sample.name))
		new_pandas.to_csv("{}/{}_abundances_csv".format(SAMP_DIR, sample.name))
