from . import download_cdna
from kinetic_datanator.util import rna_seq_util
from six.moves.urllib.request import urlretrieve
import numpy as np
import os
import pandas as pd
import requests

def get_processed_data(experiment, top_dirname):
    """ Download FASTQ files and CDNA for all samples in an expirement, and processes the data. 
            All the data is downloaded into the top directory. 

        Args:
            experiment(:obj:`array_express.Experiment`): the array express experiment
            top_dirname(:obj:`str`): the name of the directory where the overall data is being stored

        """
    exp_dirname = "{}/{}".format(top_dirname, experiment.id)
    if not os.path.isdir(exp_dirname):
        os.makedirs(exp_dirname)

    for sample in experiment.samples:
        download_cdna.run(sample, top_dirname)
        species_name = sample.ensembl_info[0].organism_strain
        sample_dirname = "{}/{}".format(exp_dirname, sample.name)
        if not os.path.isdir(sample_dirname):
            os.makedirs(sample_dirname)
        sample_url = sample.fastq_urls[0].url
        fastq_filenames = []
        for num, url in enumerate(sample.fastq_urls):
            file = urlretrieve(url.url, '{}/{}_{}.fastq.gz'.format(sample_dirname, sample.name, num)
                               )  # there used to be a space after "gz". I removed it
            fastq_filenames.append('{}/{}_{}.fastq.gz'.format(sample_dirname, sample.name, num))

        index_filename = '{}/kallisto_index_files/{}.idx'.format(top_dirname, species_name)
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
        new_pandas['percent total'] = new_column
        new_pandas.to_pickle("{}/{}_abundances_binary".format(sample_dirname, sample.name))
        new_pandas.to_csv("{}/{}_abundances_csv".format(sample_dirname, sample.name))

