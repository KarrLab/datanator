from . import download_cdna
from kinetic_datanator.util import rna_seq_util
#from six.moves.urllib.request import urlretrieve
import urllib
#from six.moves import urllib 
import numpy as np
import os
import pandas as pd
import requests
import shutil
<<<<<<< HEAD
import ftplib
import gzip

=======
>>>>>>> a6a94d8301672297d9b8fa4d2f015f34ea8033a3

def get_processed_data(experiment, top_dirname):
    """ Download FASTQ files and CDNA for all samples in an expirement, and processes the data. 
            All the data is downloaded into the top directory. 

        Args:
            experiment(:obj:`array_express.Experiment`): the array express experiment
            top_dirname(:obj:`str`): the name of the directory where the overall data is being stored

    """
    samples = experiment.samples
    get_processed_data(samples, top_dirname)
<<<<<<< HEAD



def get_processed_data_samples(samples, top_dirname):

    for sample in samples:
        if sample.ensembl_info and sample.fastq_urls:
            download_cdna.run(sample.ensembl_info[0], top_dirname)

            exp_dirname = "{}/{}".format(top_dirname, sample.experiment_id)
            if not os.path.isdir(exp_dirname):
                os.makedirs(exp_dirname)
            species_name = sample.ensembl_info[0].organism_strain
            sample_name = sample.name.replace(" ", "_")
            sample_dirname = "{}/{}".format(exp_dirname, sample_name)
            print("here: {}".format(sample_dirname))
            if not os.path.isdir(sample_dirname):
                os.makedirs(sample_dirname)
            sample_url = sample.fastq_urls[0].url
            fastq_filenames = []
            for num, url in enumerate(sample.fastq_urls):
                print("starting {}".format(num))
                print(fastq_filenames)
                file_name = '{}/{}_{}.fastq.gz'.format(sample_dirname, sample_name, num)
                fastq_filenames.append(file_name)
                file_must_be_downloaded = False
                if os.path.isfile(file_name):
                    try:
                        with gzip.open(file_name, 'rb') as f:
                            file_content = f.read()
                        file_must_be_downloaded = False
                    except:
                        file_must_be_downloaded = True
                else:
                    pass
                if file_must_be_downloaded:
                    file = urllib.urlretrieve(url.url, file_name)  # there used to be a space after "gz". I removed it
                    urllib.urlcleanup()
                else:
                    pass
                print("done with {}".format(num))


            index_filename = '{}/kallisto_index_files/{}.idx'.format(top_dirname, species_name)
            output_dirname = '{}/output'.format(sample_dirname)
            if sample.experiment.read_type == "single":
                rna_seq_util.Kallisto().quant(fastq_filenames, index_filename=index_filename, output_dirname=output_dirname,
                             single_end_reads=True, fragment_length=180, fragment_length_std=20)
            elif sample.experiment.read_type == "paired":
                rna_seq_util.Kallisto().quant(fastq_filenames, index_filename=index_filename, output_dirname=output_dirname)

=======



def get_processed_data_samples(samples, top_dirname):

    for sample in samples:
        if sample.ensembl_info and sample.fastq_urls:

            exp_dirname = "{}/{}".format(top_dirname, sample.experiment_id)
            if not os.path.isdir(exp_dirname):
                os.makedirs(exp_dirname)
            species_name = sample.ensembl_info[0].organism_strain
            sample_name = sample.name.replace(" ", "_")
            sample_dirname = "{}/{}".format(exp_dirname, sample_name)
            print("here: {}".format(sample_dirname))
            if not os.path.isdir(sample_dirname):
                os.makedirs(sample_dirname)
            sample_url = sample.fastq_urls[0].url
            fastq_filenames = []
            for num, url in enumerate(sample.fastq_urls):
                file = urlretrieve(url.url, '{}/{}_{}.fastq.gz'.format(sample_dirname, sample_name, num)
                                   )  # there used to be a space after "gz". I removed it
                fastq_filenames.append('{}/{}_{}.fastq.gz'.format(sample_dirname, sample_name, num))

            index_filename = '{}/kallisto_index_files/{}.idx'.format(top_dirname, species_name)
            output_dirname = '{}/output'.format(sample_dirname)
            if sample.experiment.read_type == "single":
                rna_seq_util.Kallisto().quant(fastq_filenames, index_filename=index_filename, output_dirname=output_dirname,
                             single_end_reads=True, fragment_length=180, fragment_length_std=20)
            elif sample.experiment.read_type == "paired":
                rna_seq_util.Kallisto().quant(fastq_filenames, index_filename=index_filename, output_dirname=output_dirname)

>>>>>>> a6a94d8301672297d9b8fa4d2f015f34ea8033a3
            new_pandas = pd.read_csv('{}/output/abundance.tsv'.format(sample_dirname), sep='\t').set_index("target_id")
            total_tpm = 0
            for num in new_pandas['tpm']:
                total_tpm = total_tpm + num
            new_column = [num/total_tpm for num in new_pandas['tpm']]
            new_pandas['target_id'] = new_column
            new_pandas.to_pickle("{}/{}_abundances_binary".format(sample_dirname, sample_name))
            new_pandas.to_csv("{}/{}_abundances_csv".format(sample_dirname, sample_name))
        else:
<<<<<<< HEAD
            print("No Ensembl Info or Fastq files")
=======
            print("No Ensembl Info or Fastq files")
>>>>>>> a6a94d8301672297d9b8fa4d2f015f34ea8033a3
