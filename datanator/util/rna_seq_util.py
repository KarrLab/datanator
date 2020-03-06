""" Utilities for RNA-seq data

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2018-01-15
:Copyright: 2018, Karr Lab
:License: MIT
"""

import subprocess

class Kallisto(object):
    """ Python interface to `kallisto <https://pachterlab.github.io/kallisto>`_. """

    def index(self, fasta_filenames, index_filename=None, kmer_size=31, make_unique=False):
        """ Generate index from FASTA files 

        Args:
            fastq_filenames (:obj:`list` of :obj:`str`): paths to FASTA files
            index_filename (:obj:`str`, optional): path to the kallisto index file to be created
            kmer_size (:obj:`int`, optional): k-mer length
            make_unique (:obj:`bool`, optional): if :obj:`True`, replace repeated target names with unique names            
        """
        # process options
        options = []

        if index_filename is not None:
            options.append('--index=' + index_filename)

        if kmer_size is not None:
            options.append('--kmer-size=' + str(kmer_size))

        if make_unique:
            options.append('--make-unique')

        # run kallisto
        self._run('index', options + fasta_filenames)

    def quant(self, fastq_filenames, index_filename=None, output_dirname=None,
              bias=False, bootstrap_samples=0, seed=42, plaintext=False, fusion=False,
              single_end_reads=False, forward_stranded=False, reverse_stranded=False,
              fragment_length=None, fragment_length_std=None,
              threads=1, pseudobam=False):
        """ Process RNA-seq FASTQ files 

        Args:
            fastq_filenames (:obj:`list` of :obj:`str`): paths to FASTQ files                
            index_filename (:obj:`str`, optional): path to the kallisto index file to be used for quantification 
            output_dirname (:obj:`str`, optional): path to the output directory

            single_end_reads (:obj:`bool`, optional): if :obj:`True`, quantify single-end reads
            fragment_length (:obj:`float`, optional): estimated average fragment length
            fragment_length_std (:obj:`float`, optional): estimated standard deviation of fragment length
        """
        # process options
        options = []

        if index_filename is not None:
            options.append('--index=' + index_filename)

        if output_dirname is not None:
            options.append('--output-dir=' + output_dirname)

        if bias:
            options.append('--bias')

        if bootstrap_samples is not None:
            options.append('--bootstrap-samples=' + str(bootstrap_samples))

        if seed is not None:
            options.append('--seed=' + str(seed))

        if plaintext:
            options.append('--plaintext')

        if fusion:
            options.append('--fusion')

        if single_end_reads:
            options.append('--single')

        if forward_stranded:
            options.append('--fr-stranded')

        if reverse_stranded:
            options.append('--rf-stranded')

        if fragment_length is not None:
            options.append('--fragment-length=' + str(fragment_length))

        if fragment_length_std is not None:
            options.append('--sd=' + str(fragment_length_std))

        if threads is not None:
            options.append('--threads=' + str(threads))

        if pseudobam:
            options.append('--pseudobam')

        # run kallisto
        self._run('quant', options + fastq_filenames)

    def _run(self, cmd, args, verbose=False):
        """ Run a kallisto command

        Args:
            cmd (:obj:`str`): kallisto command
            args (:obj:`list` of :obj:`str`): arguments to the kallisto command
            verbose (:obj:`bool`, optional): if :obj:`True`, write status information to stdout
        """
        if verbose:
            stdout = None
            stderr = None
        else:
            stdout = subprocess.DEVNULL
            stderr = subprocess.DEVNULL
        subprocess.check_call(['kallisto', cmd] + args, stdout=stdout, stderr=stderr)
        