#from . import download_cdna
from kinetic_datanator.util import rna_seq_util
#from six.moves.urllib.request import urlretrieve
import urllib
#from six.moves import urllib 
import numpy as np
import os
import pandas as pd
import requests
import shutil
import ftplib
import gzip



def get_processed_data_samples(samples, output_directory, temp_directory):
    for sample in samples:
        if sample.ensembl_info and sample.fastq_urls:
            download_cdna(sample.ensembl_info[0], temp_directory)
            download_all_fastq([sample], temp_directory)
            process_cdna(sample.ensembl_info[0], output_directory, temp_directory)
            process_fastq([sample], output_directory, temp_directory)
            delete_cdna_files(sample.ensembl_info[0], temp_directory)
            delete_fastq_files([sample], temp_directory)
        else:
            print("No FASTQ or no Ensembl for {}_{}".format(sample.experiment_id, sample.name))


def download_cdna(ensembl_info, temp_directory):
    CDNA_DIR = "{}/CDNA_FILES".format(temp_directory)
    if not os.path.isdir(CDNA_DIR):
        os.makedirs(CDNA_DIR)
    file_name = "{}/{}.cdna.all.fa.gz".format(CDNA_DIR, ensembl_info.organism_strain)
    if not os.path.isfile(file_name):
        file = urlretrieve(ensembl_info.url, file_name)
    os.chdir(temp_directory)


def download_all_fastq(samples, temp_directory):
    for sample in samples:
        for num, url in enumerate(sample.fastq_urls):
            print("starting {}".format(num))
            file_name = '{}/FASTQ_Files/{}__{}__{}.fastq.gz'.format(temp_directory, sample.experiment_id, sample.name, num)
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

def process_cdna(ensembl_info, output_directory, temp_directory):

    CDNA_DIR = "{}/CDNA_FILES".format(temp_directory)
    cdna_file = "{}/{}.cdna.all.fa.gz".format(CDNA_DIR, ensembl_info.organism_strain)
    KALLISTO_DIR = "{}/kallisto_index_files".format(output_directory)
    if not os.path.isdir(KALLISTO_DIR):
        os.makedirs(KALLISTO_DIR)
    kallisto_file = "{}/{}.idx".format(KALLISTO_DIR, ensembl_info.organism_strain)
    if not os.path.isfile(kallisto_file):
        os.system("kallisto index -i {} {}".format(kallisto_file, cdna_file))


def process_fastq(samples, output_directory, temp_directory):
    for sample in samples:
        if sample.ensembl_info and sample.fastq_urls:

            exp_dirname = "{}/{}".format(output_directory, sample.experiment_id)
            if not os.path.isdir(exp_dirname):
                os.makedirs(exp_dirname)

            sample_dirname = "{}/{}".format(exp_dirname, sample.name)
            if not os.path.isdir(sample_dirname):
                os.makedirs(sample_dirname)

            fastq_filenames = []
            for num in range(len(sample.fastq_urls)):
                fastq_filenames.append("{}/FASTQ_Files/{}__{}__{}.fastq.gz".format(temp_directory, sample.experiment_id, sample.name, num))
            index_filename = '{}/kallisto_index_files/{}.idx'.format(output_directory, sample.ensembl_info[0].organism_strain)
            output_dirname = '{}/output'.format(sample_dirname)
            if sample.experiment.read_type == "single":
                    rna_seq_util.Kallisto().quant(fastq_filenames, index_filename=index_filename, output_dirname=output_dirname,
                                 single_end_reads=True, fragment_length=180, fragment_length_std=20)
            elif sample.experiment.read_type == "paired":
                rna_seq_util.Kallisto().quant(fastq_filenames, index_filename=index_filename, output_dirname=output_dirname)

            new_pandas = pd.read_csv('{}/output/abundance.tsv'.format(sample_dirname), sep='\t').set_index("target_id")
            total_tpm = 0
            for num in new_pandas['tpm']:
                total_tpm = total_tpm + num
            new_column = [num/total_tpm for num in new_pandas['tpm']]
            new_pandas['target_id'] = new_column
            new_pandas.to_pickle("{}/{}_abundances_binary".format(sample_dirname, sample.name))
            #new_pandas.to_csv("{}/{}_abundances_csv".format(sample_dirname, sample.name))
        else:
            print("No Ensembl Info or Fastq files")



def delete_cdna_files(ensembl_info, temp_directory):
    file_name = "{}/CDNA_FILES/{}.cdna.all.fa.gz".format(temp_directory, ensembl_info.organism_strain)
    if os.path.isfile(file_name):
        os.remove(file_name)


def delete_fastq_files(samples, temp_directory):
    for sample in samples:
        for num in range(len(sample.fastq_urls)):
            file_name = "{}/FASTQ_Files/{}__{}__{}.fastq.gz".format(temp_directory, sample.experiment_id, sample.name, num)
            if os.path.isfile(file_name):
                os.remove(file_name)

