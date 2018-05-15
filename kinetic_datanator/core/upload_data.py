from kinetic_datanator.data_source import refseq
import os
from os import path
from Bio import SeqIO


CACHE_DIRNAME = os.path.join(os.path.dirname(__file__), '..', 'data_source', 'cache')


class UploadData():
	def __init__(self, cache_dirname=CACHE_DIRNAME):
		self.cache_dirname = cache_dirname


	def upload_reference_genome(self, path_to_annotation_file):
		#bio_seqio_object = SeqIO.parse(path_to_annotation_file, "genbank")
		refseq.Refseq(cache_dirname=self.cache_dirname).load_content([SeqIO.parse(path_to_annotation_file, "genbank")])



	def upload_processed_data(sample_name, csv_file, top_directory=path.dirname(__file__)):
	    #this takes in a csv file, and puts the file in the proper place so that it can be queried. 
	    new_pandas = pd.read_csv(csv_file, sep=',').set_index("gene_locus")
	    new_pandas.to_pickle("{}/{}".format(top_directory, sample_name))

