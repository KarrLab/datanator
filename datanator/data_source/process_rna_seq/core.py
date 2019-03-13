#from . import download_cdna
from datanator.util import rna_seq_util
from six.moves.urllib.request import urlretrieve
import urllib
#from six.moves import urllib 
import numpy as np
import os
import pandas as pd
import requests
import shutil
import ftplib
import gzip
from Bio import SeqIO





def get_processed_data_samples(samples, output_directory, temp_directory):
    for sample in samples:
        if sample.ensembl_info and sample.fastq_urls:
            download_cdna(sample.ensembl_info[0].ref_genome, sample.ensembl_info[0].organism_strain, sample.ensembl_info[0].url, temp_directory)
            fastq_urls = ""
            for url in sample.fastq_urls:
                fastq_urls = fastq_urls + url.url + " "
            fastq_urls = fastq_urls[:-1]
            download_fastq(sample.experiment_id, sample.name, temp_directory, fastq_urls)
            process_cdna(sample.ensembl_info[0].organism_strain, output_directory, temp_directory)
            process_fastq(sample.experiment_id, sample.name, sample.ensembl_info[0].organism_strain, len(sample.fastq_urls), sample.experiment.read_type, output_directory, temp_directory)
            #delete_cdna_files(sample.ensembl_info[0].organism_strain, temp_directory)
            delete_fastq_files(sample.experiment_id, sample.name, temp_directory)
        else:
            print("No FASTQ or no Ensembl for {}_{}".format(sample.experiment_id, sample.name))


def download_cdna(ref_genome, strain_name, url, temp_directory):
    CDNA_DIR = "{}/CDNA_FILES".format(temp_directory)
    if not os.path.isdir(CDNA_DIR):
        os.makedirs(CDNA_DIR)
    file_name = "{}/{}.cdna.all.fa.gz".format(CDNA_DIR, strain_name)
    if not os.path.isfile(file_name):
        if ref_genome == "ensembl":
            if not os.path.isfile(file_name):
                print(url)
                file = urlretrieve(url, file_name)
            os.chdir(temp_directory)
        elif ref_genome == "genbank":
            new_cdna_file = gzip.open(file_name, 'wb')
            temp_file_name = "{}/{}.temp.all.fa.gz".format(CDNA_DIR, strain_name)
            file = urlretrieve(url, temp_file_name)
            list_locus_tags = []
            with gzip.open(temp_file_name, "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    locus_tag = record.description[record.description.find("[locus_tag=")+11:record.description.find("]", record.description.find("[locus_tag="))]
                    if locus_tag in list_locus_tags:
                        locus_tag = "{}_variation".format(locus_tag)
                    list_locus_tags.append(locus_tag)
                    new_cdna_file.write(bytes(">{}\n".format(locus_tag), 'utf-8'))
                    new_cdna_file.write(bytes("{}\n".format(record.seq), 'utf-8'))
    #print(locus_tag)
            os.remove(temp_file_name)



#def download_fastq(sample, temp_directory):
def download_fastq(experiment_name,  sample_name, temp_directory, fastq_urls):
    for num, url in enumerate(fastq_urls.split(" ")):
        print("starting {}".format(num))
        file_name = '{}/FASTQ_Files/{}__{}__{}.fastq.gz'.format(temp_directory, experiment_name, sample_name, num)
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

def process_cdna(strain_name, output_directory, temp_directory):

    CDNA_DIR = "{}/CDNA_FILES".format(temp_directory)
    cdna_file = "{}/{}.cdna.all.fa.gz".format(CDNA_DIR, strain_name)
    KALLISTO_DIR = "{}/kallisto_index_files".format(output_directory)
    if not os.path.isdir(KALLISTO_DIR):
        os.makedirs(KALLISTO_DIR)
    kallisto_file = "{}/{}.idx".format(KALLISTO_DIR, strain_name)
    if not os.path.isfile(kallisto_file):
        os.system("kallisto index -i {} {}".format(kallisto_file, cdna_file))


def process_fastq(experiment_name, sample_name, strain_name, num_fastq_files, read_type, output_directory, temp_directory):

    exp_dirname = "{}/{}".format(output_directory, experiment_name)
    if not os.path.isdir(exp_dirname):
        os.makedirs(exp_dirname)

    sample_dirname = "{}/{}".format(exp_dirname, sample_name)
    if not os.path.isdir(sample_dirname):
        os.makedirs(sample_dirname)

    fastq_filenames = []
    for num in range(num_fastq_files):
        fastq_filenames.append("{}/FASTQ_Files/{}__{}__{}.fastq.gz".format(temp_directory, experiment_name, sample_name, num))
    index_filename = '{}/kallisto_index_files/{}.idx'.format(output_directory, strain_name)
    output_dirname = '{}/output'.format(sample_dirname)
    if read_type == "single":
            rna_seq_util.Kallisto().quant(fastq_filenames, index_filename=index_filename, output_dirname=output_dirname,
                         single_end_reads=True, fragment_length=180, fragment_length_std=20)
    elif read_type == "paired":
        rna_seq_util.Kallisto().quant(fastq_filenames, index_filename=index_filename, output_dirname=output_dirname)

    new_pandas = pd.read_csv('{}/output/abundance.tsv'.format(sample_dirname), sep='\t').set_index("target_id")
    total_tpm = 0
    for num in new_pandas['tpm']:
        total_tpm = total_tpm + num
    new_column = [num/total_tpm for num in new_pandas['tpm']]
    new_pandas['target_id'] = new_column
    new_pandas.to_pickle("{}/{}_abundances_binary".format(sample_dirname, sample_name))
    #new_pandas.to_csv("{}/{}_abundances_csv".format(sample_dirname, sample.name))



def delete_cdna_files(strain_name, temp_directory):
    file_name = "{}/CDNA_FILES/{}.cdna.all.fa.gz".format(temp_directory, strain_name)
    if os.path.isfile(file_name):
        os.remove(file_name)


def delete_fastq_files(experiment_name, sample_name, temp_directory):
    for file in os.listdir("{}/FASTQ_Files".format(temp_directory)):
        if file.startswith("{}__{}".format(experiment_name, sample_name)):
            os.remove("{}/FASTQ_Files/{}".format(temp_directory, file))
