""" Test for RNA-seq utilities

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-01-15
:Copyright: 2018, Karr Lab
:License: MIT
"""

from datanator.util import rna_seq_util
from six.moves import urllib
import capturer
import os
import shutil
import tempfile
import unittest


class TestKallisto(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test(self):
        fasta_filename = os.path.join(self.temp_dir, 'transcripts.fasta.gz')
        fastq_filenames = [
            os.path.join(self.temp_dir, 'reads_1.fastq.gz'),
            os.path.join(self.temp_dir, 'reads_2.fastq.gz'),
        ]

        # download files
        urllib.request.urlretrieve('https://github.com/pachterlab/kallisto/raw/master/test/transcripts.fasta.gz',
                                   fasta_filename)
        urllib.request.urlretrieve('https://github.com/pachterlab/kallisto/raw/master/test/reads_1.fastq.gz',
                                   fastq_filenames[0])
        urllib.request.urlretrieve('https://github.com/pachterlab/kallisto/raw/master/test/reads_2.fastq.gz',
                                   fastq_filenames[1])

        # run kallisto on test files
        index_filename = os.path.join(self.temp_dir, 'index.idx')
        rna_seq_util.Kallisto().index([fasta_filename], index_filename=index_filename)
        self.assertTrue(os.path.isfile(index_filename))

        output_dirname = os.path.join(self.temp_dir, 'out')
        rna_seq_util.Kallisto().quant(fastq_filenames, index_filename=index_filename, output_dirname=output_dirname)
        self.assertTrue(os.path.isdir(output_dirname))
        self.assertTrue(os.path.isfile(os.path.join(output_dirname, 'abundance.tsv')))
        self.assertTrue(os.path.isfile(os.path.join(output_dirname, 'abundance.h5')))
        self.assertTrue(os.path.isfile(os.path.join(output_dirname, 'run_info.json')))

    def test_error(self):
        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with self.assertRaises(Exception):
                rna_seq_util.Kallisto()._run('index', ['__undefined__.fasta'], verbose=True)
            self.assertNotEqual(captured.stdout.get_text(), '')
            self.assertNotEqual(captured.stderr.get_text(), '')

        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with self.assertRaises(Exception):
                rna_seq_util.Kallisto()._run('index', ['__undefined__.fasta'], verbose=False)
            self.assertEqual(captured.stdout.get_text(), '')
            self.assertEqual(captured.stderr.get_text(), '')
