""" Test the ArrayExpress database
:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-08-16
:Copyright: 2017, Karr Lab
:License: MIT
"""

from datanator.data_source import array_express
from datanator.data_source.array_express_tools import ensembl_tools
from datanator.data_source.process_rna_seq import core
from datanator.data_source.process_rna_seq import download_cdna
from six.moves.urllib.request import urlretrieve
import datetime
import os
import pandas
import shutil
import tempfile
import unittest

#@unittest.skip('skip')
class QuickTest(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.src = array_express.ArrayExpress(cache_dirname=self.cache_dirname, download_backups=False, load_content=False)

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_load_experiment_metadata(self):
        src = self.src
        session = src.session

        self.assertEqual.__self__.maxDiff = None

        src.load_experiment_metadata(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/?date=[2001-01-01+2002-12-31]")
        q = session.query(array_express.Experiment)
        self.assertGreater(q.count(), 3)
        self.assertEqual(session.query(array_express.Experiment).filter_by(id='E-GEOD-6').count(), 1)
        self.assertEqual(session.query(array_express.Experiment).filter_by(id='E-GEOD-8').count(), 1)
        self.assertEqual(session.query(array_express.Experiment).filter_by(id='E-GEOD-10').count(), 1)

        experiment = session.query(array_express.Experiment).filter_by(id='E-GEOD-10').first()
        self.assertEqual(experiment.name, 'eye-SAGE')
        self.assertEqual(experiment.name_2, None)
        self.assertEqual(experiment.organisms, [session.query(array_express.Organism).filter_by(name='Homo sapiens').first()])
        self.assertTrue(experiment.description.startswith('Human retinal and RPE SAGE libraries.'))
        self.assertTrue(experiment.description.endswith('. Keywords: other'))
        self.assertEqual(experiment.types, [session.query(array_express.ExperimentType).filter_by(
            name='transcription profiling by SAGE').first()])
        self.assertEqual(experiment.designs, [])
        self.assertEqual(sorted([df.name for df in experiment.data_formats]), ['normalization', 'processedData'])
        self.assertEqual(sorted([df.name for df in experiment.data_formats]), sorted([df.name for df in [
            session.query(array_express.DataFormat).filter_by(name='normalization').first(),
            session.query(array_express.DataFormat).filter_by(name='processedData').first(),
        ]]))

        self.assertEqual(experiment.submission_date, datetime.date(2001, 10, 2))
        self.assertEqual(experiment.release_date, datetime.date(2001, 10, 3))

        experiment = session.query(array_express.Experiment).filter_by(id='E-SNGR-7').first()
        self.assertEqual(experiment.designs, [session.query(array_express.ExperimentDesign).filter_by(name='time series').first()])

    @unittest.skip("need to use an example of RNA-SEQ sample")
    def test_load_experiment_samples(self):

        src = self.src
        session = src.session

        src.load_experiment_metadata(src.load_experiment_metadata(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/?date=[2001-01-01+2002-12-31]"))

        # E-GEOD-10
        experiment = session.query(array_express.Experiment).filter_by(id='E-GEOD-10').first()
        src.load_experiment_samples(experiment)

        self.assertEqual(len(experiment.samples), 4)
        self.assertEqual(set([sample.name for sample in experiment.samples]), set(['GSM571 1', 'GSM572 1', 'GSM573 1', 'GSM574 1']))

        sample = next(sample for sample in experiment.samples if sample.name == 'GSM574 1')
        self.assertEqual(sample.extracts, [session.query(array_express.Extract).filter_by(name='GSM574 extract 1').first()])
        self.assertEqual(sample.assay, 'GSM574')
        self.assertEqual(sample.characteristics, [session.query(
            array_express.Characteristic).filter_by(category='Organism', value='Homo sapiens').first()])
        self.assertEqual(sample.variables, [])

        # E-SNGR-6
        experiment = session.query(array_express.Experiment).filter_by(id='E-SNGR-6').first()
        src.load_experiment_samples(experiment)

        self.assertEqual(len(experiment.samples), 1)

        sample = next(sample for sample in experiment.samples if sample.index == 0)
        self.assertEqual(sample.name, 'Schizosaccharomyces pombe, cultured with nitrogen source NH4Cl')
        self.assertEqual(sample.extracts, [
            session.query(array_express.Extract).filter_by(name='1h, added nitrogen').first()
        ])
        self.assertEqual(sample.assay, 'jm_1h_135-15')

        self.assertEqual(sorted(sample.characteristics, key=lambda characteristic: characteristic.category), [
            session.query(array_express.Characteristic).filter_by(category='Genotype', value='h+/h+ ade6-M210/ade6-M216').first(),
            session.query(array_express.Characteristic).filter_by(category='Organism', value='Schizosaccharomyces pombe').first(),
        ])
        self.assertEqual(sample.variables, [
            session.query(array_express.Variable).filter_by(name='sampling interval', value='1', unit=None).first()
        ])

    def test_load_experiment_protocols(self):
        src = self.src
        session = src.session

        src.load_experiment_metadata(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/?date=[2001-01-01+2002-12-31]")

        # E-GEOD-10
        experiment = session.query(array_express.Experiment).filter_by(id='E-GEOD-6').first()
        src.load_experiment_protocols(experiment)
        self.assertEqual(len(experiment.protocols), 2)
        self.assertEqual(set([protocol.protocol_accession for protocol in experiment.protocols]), set(['P-GSE6-1', 'P-GSE6-2', ]))

        q = session.query(array_express.Protocol)
        protocol = session.query(array_express.Protocol).filter_by(protocol_accession='P-GSE6-1').first()
        text = "ID_REF = \nCH1B_MEAN = mean channel 1 background\nCH1B_MEDIAN = median channel 1 background\nCH1D_MEAN = difference between CH1I_MEAN and CH1B_MEAN\nCH2I_MEAN = mean channel 2 signal\nCH2D_MEAN = difference of CH2I_MEAN and CH2B_MEAN\nCH2B_MEAN = mean of channel 2 background\nCH2B_MEDIAN = median of channel 2 background\nCH2BN_MEDIAN = normalized CH2B_MEDIAN\nCH2DN_MEAN = normalized CH2D_MEAN\nCH2IN_MEAN = normalized CH2I_MEAN\nCORR = simple correlation coefficient of channel 1 and channel 2 pixels\nFLAG = user defined flag (default=0)\nVALUE = same as UNF_VALUE but with flagged values removed\nPIX_RAT2_MEDIAN =  \nPERGTBCH1I_1SD = standard deviation of fraction of pixels greater than CH1B_MEAN\nPERGTBCH2I_1SD = standard deviation of fraction of pixels greater than CH2B_MEAN\nRAT1_MEAN = ratio of CH1D_MEAN to CH2D_MEAN\nRAT1N_MEAN = ratio of CH1DN_MEAN to CH2DN_MEAN\nRAT2_MEAN = ratio of CH2D_MEAN to CH1D_MEAN\nRAT2N_MEAN = ratio of CH2DN_MEAN and CH1DN_MEAN\nREGR = slope of the simple regression line of channel 2 to channel 1 pixels\nTOT_BPIX = total number of background pixels\nTOT_SPIX = total number of signal (spot) pixels\nTOP = vertical, short-axis spot ellipse minimum pixel coordinate\nBOT = vertical, short-axis spot ellipse maximum pixel coordinate\nLEFT = horizontal, long-axis spot ellipse minimum pixel coordinate\nRIGHT = horizontal, long-axis spot ellipse maximum pixel coordinate\nUNF_VALUE = aka LOG_RAT2N_MEAN, or log2 of ratio of CH2DN_MEAN and CH1DN_MEAN"
        self.assertEqual(protocol.text, text)

        test_experiment = array_express.Experiment(id='E-MTAB-5281')
        src.load_experiment_protocols(test_experiment)
        protocol = session.query(array_express.Protocol).filter_by(protocol_accession='P-MTAB-42502').first()
        self.assertEqual(protocol.protocol_type, 'labelling')
        self.assertEqual(protocol.text, 'Labelling was performed by using enzymatic attachment of nucleotides coupled to biotin.')
        self.assertEqual(protocol.performer, 'Stephanie Boue')
        self.assertEqual(protocol.hardware, 'Affymetrix GeneChip Scanner 3000 7G')
        self.assertEqual(protocol.software, 'Affymetrix AGCC')

    def test_load_content(self):
        src = self.src
        session = src.session
        # print(len(session.query(array_express.Experiment).all()))
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/E-MTAB-4262")
        exp = session.query(array_express.Experiment).filter_by(id="E-MTAB-4262").first()
        self.assertEqual(exp.id, "E-MTAB-4262")
        self.assertEqual(exp.name, "Expression analysis of Arabidopsis seedlings 3h and 2h after transfer to sucrose")
        self.assertEqual(exp.description, "Seedlings are grown on a mesh covering MS media without carbon source, and subsequently transferred at 9 DAS to control media without sucrose or medium supplemented with 15 mM sucrose. 3 hours and 24 hours after transfer, seedlings were harvested and the third leaf micro-dissected for RNA extraction.")
        self.assertEqual(exp.organisms[0].name, "Arabidopsis thaliana")
        self.assertEqual(list(set([t.name for t in exp.types])), list(set(["RNA-seq of coding RNA"])))
        self.assertEqual(list(set([t.name for t in exp.designs])), list(set(["compound treatment design", "time series design"])))
        self.assertEqual(len(exp.protocols), 8)
        self.assertTrue(exp.has_fastq_files)
        self.assertEqual(exp.read_type, "paired")
        self.assertEqual(len(exp.samples), 12)
        a_sample = session.query(array_express.Sample).filter_by(experiment=exp).filter_by(name="R1_0mM_24h").first()
        self.assertEqual(a_sample.read_type, "paired")
        self.assertEqual(list(set([
            "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR173/002/ERR1736192/ERR1736192_1.fastq.gz",
            "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR173/002/ERR1736192/ERR1736192_2.fastq.gz",
        ])), list(set([u.url for u in a_sample.fastq_urls])))

@unittest.skip('skip')
class TestLoadFASTQ_Url(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.src = array_express.ArrayExpress(cache_dirname=self.cache_dirname, download_backups=False, load_content=False)

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_single_read(self):
        session = self.src.session
        self.src.load_experiment_samples(array_express.Experiment(id="E-GEOD-41373"))
        exp = session.query(array_express.Experiment).filter_by(id="E-GEOD-41373").first()
        self.assertTrue(exp.has_fastq_files)
        self.assertEqual(exp.read_type, "single")
        self.assertEqual(len(exp.samples), 12)
        a_sample = session.query(array_express.Sample).filter_by(experiment=exp).filter_by(name="GSM1015798_1").first()
        self.assertEqual(a_sample.read_type, "single")
        self.assertEqual(a_sample.fastq_urls[0].url, "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR580/SRR580556/SRR580556.fastq.gz")


@unittest.skip('skip')
class TestDownloadCDNA(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.src = array_express.ArrayExpress(cache_dirname=self.cache_dirname, download_backups=False, load_content=False)

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_get_cdna_prokaryotes(self):
        # the strain is added to the organism and found
        src = self.src
        session = src.session
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/E-MTAB-6099")
        exp = session.query(array_express.Experiment).filter_by(id="E-MTAB-6099").first()
        ensembl_info = exp.samples[0].ensembl_info[0]
        download_cdna.run(ensembl_info, self.cache_dirname)
        self.assertTrue(os.path.isfile('{}/CDNA_FILES/burkholderia_cenocepacia_j2315.cdna.all.fa.gz'.format(self.cache_dirname)))
        self.assertTrue(os.path.isfile('{}/kallisto_index_files/burkholderia_cenocepacia_j2315.idx'.format(self.cache_dirname)))

@unittest.skip('skip')
class TestEnsemblTools(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.src = array_express.ArrayExpress(cache_dirname=self.cache_dirname, download_backups=False, load_content=False)

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_prokaryote_1(self):
        src = self.src
        session = src.session
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/E-MTAB-5971")
        exp = session.query(array_express.Experiment).filter_by(id="E-MTAB-5971").first()
        sample = exp.samples[0]
        self.assertEqual(sample.ensembl_info[0].organism_strain, "escherichia_coli_k_12_mg1655")
        self.assertTrue(sample.full_strain_specificity)
        self.assertEqual(sample.ensembl_info[0].url, "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz")

    def test_prokaryote_2(self):
        src = self.src
        session = src.session
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/E-GEOD-68277")
        exp = session.query(array_express.Experiment).filter_by(id="E-GEOD-68277").first()
        sample = exp.samples[0]
        self.assertEqual(sample.ensembl_info[0].organism_strain, "streptococcus_pyogenes_m1_gas_(serotype_m1)")
        self.assertFalse(sample.full_strain_specificity)
        self.assertEqual(sample.ensembl_info[0].url, "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/785/GCF_000006785.2_ASM678v2/GCF_000006785.2_ASM678v2_genomic.gbff.gz")

    def test_eukaryote(self):
        src = self.src
        session = src.session
        ax = "E-MTAB-5149"
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}".format(ax))
        exp = session.query(array_express.Experiment).filter_by(id=ax).first()
        sample = exp.samples[0]
        self.assertEqual(sample.ensembl_info[0].organism_strain, "saccharomyces_cerevisiae")
        self.assertTrue(sample.full_strain_specificity)
        self.assertEqual(sample.ensembl_info[0].url,
                         "ftp://ftp.ensembl.org/pub/current_fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz")

    def test_plant(self):
        src = self.src
        session = src.session
        ax = "E-GEOD-61293"
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}".format(ax))
        exp = session.query(array_express.Experiment).filter_by(id=ax).first()
        sample = exp.samples[0]
        self.assertEqual(sample.ensembl_info[0].organism_strain, "arabidopsis_thaliana")
        self.assertTrue(sample.full_strain_specificity)
        self.assertEqual(sample.ensembl_info[
                         0].url, "ftp://ftp.ensemblgenomes.org/pub/current/plants/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz")

    def test_error_multiple_organisms(self):
        session = self.src.session
        self.src.load_experiment_samples(array_express.Experiment(id="E-GEOD-58388"))
        exp = session.query(array_express.Experiment).filter_by(id="E-GEOD-58388").first()
        a_sample = session.query(array_express.Sample).filter_by(experiment=exp).filter_by(name="GSM1409710_1").first()
        with self.assertRaises(LookupError):
            ensembl_tools.get_strain_info(a_sample)

    def test_error(self):
        src = self.src
        session = src.session
        ax = "E-GEOD-58388"
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}".format(ax))
        exp = session.query(array_express.Experiment).filter_by(id=ax).first()
        sample = exp.samples[0]
        self.assertEqual(sample.ensembl_info, [])

@unittest.skip('handling of files not compatible with circleci')
class TestProcessData(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = '{}/test_processing'.format(os.path.dirname(os.path.realpath(__file__)))
        self.src = array_express.ArrayExpress(cache_dirname=self.cache_dirname, download_backups=False, load_content=False)

    def tearDown(self):
        os.unlink('{}/ArrayExpress.sqlite'.format(self.cache_dirname))
        os.unlink('{}/E-MTAB-6099/Control_2/Control_2_abundances_binary'.format(self.cache_dirname))
        shutil.rmtree("{}/kallisto_index_files".format(self.cache_dirname))
        #shutil.copy("{}/backup_temporary_files/CDNA_FILES/burkholderia_cenocepacia_j2315.cdna.all.fa.gz".format(self.cache_dirname), "{}/temporary_files/CDNA_FILES/burkholderia_cenocepacia_j2315.cdna.all.fa.gz".format(self.cache_dirname))
        shutil.copy("{}/backup_temporary_files/FASTQ_Files/E-MTAB-6099__Control_2__0.fastq.gz".format(self.cache_dirname), "{}/temporary_files/FASTQ_Files/E-MTAB-6099__Control_2__0.fastq.gz".format(self.cache_dirname))

    def test_process_data(self):
        src = self.src
        session = src.session
        ax = "E-MTAB-6099"
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}".format(ax))
        sample_name = 'Control_2'
        sample = session.query(array_express.Sample).filter_by(name=sample_name).first()
        core.get_processed_data_samples([sample], self.cache_dirname, "{}/temporary_files".format(self.cache_dirname))
        self.assertTrue(os.path.isfile("{}/E-MTAB-6099/Control_2/Control_2_abundances_binary".format(self.cache_dirname)))
        data = pandas.read_pickle("{}/E-MTAB-6099/Control_2/Control_2_abundances_binary".format(self.cache_dirname))
        self.assertEqual(data.loc["BCAL0002", "est_counts"], 0)
        self.assertEqual(data.loc["BCAL0022", "tpm"], 0)
        self.assertFalse(os.path.isfile("{}/temporary_files/FASTQ_Files/E-MTAB-6099__Control_2__0.fastq.gz".format(self.cache_dirname)))
        #self.assertFalse(os.path.isfile("{}/temporary_files/CDNA_FILES/burkholderia_cenocepacia_j2315.cdna.all.fa.gz".format(self.cache_dirname)))


@unittest.skip('skip for other tests')
class TestCommandLineProcessing(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = '{}/test_processing'.format(os.path.dirname(os.path.realpath(__file__)))
        self.src = array_express.ArrayExpress(cache_dirname=self.cache_dirname, download_backups=False, load_content=False)

    def tearDown(self):
        os.unlink('{}/ArrayExpress.sqlite'.format(self.cache_dirname))
        os.unlink('{}/E-MTAB-6099/Control_2/Control_2_abundances_binary'.format(self.cache_dirname))
        shutil.rmtree("{}/kallisto_index_files".format(self.cache_dirname))
        shutil.copy("{}/backup_temporary_files/CDNA_FILES/burkholderia_cenocepacia_j2315.cdna.all.fa.gz".format(self.cache_dirname), "{}/temporary_files/CDNA_FILES/burkholderia_cenocepacia_j2315.cdna.all.fa.gz".format(self.cache_dirname))
        shutil.copy("{}/backup_temporary_files/FASTQ_Files/E-MTAB-6099__Control_2__0.fastq.gz".format(self.cache_dirname), "{}/temporary_files/FASTQ_Files/E-MTAB-6099__Control_2__0.fastq.gz".format(self.cache_dirname))

    @unittest.skip("this test is slow")
    def test_process_data_single(self):
        src = self.src
        session = src.session
        ax = "E-MTAB-6099"
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}".format(ax))
        sample_name = 'Control_2'
        sample = session.query(array_express.Sample).filter_by(name=sample_name).first()
        #core.get_processed_data_samples([sample], self.cache_dirname, "{}/temporary_files".format(self.cache_dirname))
        python_file = "python datanator/data_source/process_rna_seq/command_line_core.py"
        output_directory = self.cache_dirname
        temp_directory = "{}/temporary_files".format(self.cache_dirname)

        os.system("{} download-cdna {} {} {}".format(python_file, sample.ensembl_info[0].organism_strain, sample.ensembl_info[0].url, temp_directory))

        fastq_urls = ""
        for url in sample.fastq_urls:
            fastq_urls = fastq_urls + url.url + " "
        fastq_urls = fastq_urls[:-1]

        os.system("""{} download-fastq {} {} {} "{}" """.format(python_file, sample.experiment_id, sample.name, temp_directory, fastq_urls))

        os.system("{} process-cdna {} {} {}".format(python_file, sample.ensembl_info[0].organism_strain, output_directory, temp_directory))

        os.system("{} process-fastq {} {} {} {} {} {} {}".format(python_file, sample.experiment_id, sample.name, sample.ensembl_info[0].organism_strain, len(sample.fastq_urls), sample.experiment.read_type, output_directory, temp_directory))

        os.system("{} delete-cdna {} {}".format(python_file, sample.ensembl_info[0].organism_strain, temp_directory))

        os.system("{} delete-fastq {} {} {}".format(python_file, sample.experiment_id, sample.name, temp_directory))



        self.assertTrue(os.path.isfile("{}/E-MTAB-6099/Control_2/Control_2_abundances_binary".format(self.cache_dirname)))
        data = pandas.read_pickle("{}/E-MTAB-6099/Control_2/Control_2_abundances_binary".format(self.cache_dirname))
        self.assertEqual(data.loc["CAO00538", "est_counts"], 2)
        self.assertEqual(data.loc["CAO00538", "tpm"], 361319)
        self.assertFalse(os.path.isfile("{}/temporary_files/FASTQ_Files/E-MTAB-6099__Control_2__0.fastq.gz".format(self.cache_dirname)))
        self.assertFalse(os.path.isfile("{}/temporary_files/CDNA_FILES/burkholderia_cenocepacia_j2315.cdna.all.fa.gz".format(self.cache_dirname)))

    def test_process_data(self):
        src = self.src
        session = src.session
        ax = "E-MTAB-6099"
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}".format(ax))
        sample_name = 'Control_2'
        sample = session.query(array_express.Sample).filter_by(name=sample_name).first()
        #core.get_processed_data_samples([sample], self.cache_dirname, "{}/temporary_files".format(self.cache_dirname))
        python_file = "python datanator/data_source/process_rna_seq/command_line_core.py"
        output_directory = self.cache_dirname
        temp_directory = "{}/temporary_files".format(self.cache_dirname)

        os.system("{} download-cdna {} {} {}".format(python_file, sample.ensembl_info[0].organism_strain, sample.ensembl_info[0].url, temp_directory))

        fastq_urls = ""
        for url in sample.fastq_urls:
            fastq_urls = fastq_urls + url.url + " "
        fastq_urls = fastq_urls[:-1]

        os.system("""{} download-fastq {} {} {} "{}" """.format(python_file, sample.experiment_id, sample.name, temp_directory, fastq_urls))

        os.system("{} process-cdna {} {} {}".format(python_file, sample.ensembl_info[0].organism_strain, output_directory, temp_directory))

        os.system("{} process-fastq {} {} {} {} {} {} {}".format(python_file, sample.experiment_id, sample.name, sample.ensembl_info[0].organism_strain, len(sample.fastq_urls), sample.experiment.read_type, output_directory, temp_directory))

        os.system("{} delete-cdna {} {}".format(python_file, sample.ensembl_info[0].organism_strain, temp_directory))

        os.system("{} delete-fastq {} {} {}".format(python_file, sample.experiment_id, sample.name, temp_directory))



        self.assertTrue(os.path.isfile("{}/E-MTAB-6099/Control_2/Control_2_abundances_binary".format(self.cache_dirname)))
        data = pandas.read_pickle("{}/E-MTAB-6099/Control_2/Control_2_abundances_binary".format(self.cache_dirname))
        self.assertEqual(data.loc["CAO00538", "est_counts"], 2)
        self.assertEqual(data.loc["CAO00538", "tpm"], 361319)
        self.assertFalse(os.path.isfile("{}/temporary_files/FASTQ_Files/E-MTAB-6099__Control_2__0.fastq.gz".format(self.cache_dirname)))
        self.assertFalse(os.path.isfile("{}/temporary_files/CDNA_FILES/burkholderia_cenocepacia_j2315.cdna.all.fa.gz".format(self.cache_dirname)))


    #def test_weird(self):
    #    src = self.src
    #    session = src.session
    #    ax = "E-GEOD-32228"
    #    ax = "E-GEOD-50870"
    #    src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}".format(ax))
    #    exp = session.query(array_express.Experiment).filter_by(id=ax).first()
    #    sample = exp.samples[0]
        #self.assertEqual(sample.ensembl_info[0].organism_strain, "saccharomyces_cerevisiae")
        #self.assertTrue(sample.full_strain_specificity)
        #self.assertEqual(sample.ensembl_info[0].url,
        #                 "ftp://ftp.ensembl.org/pub/current_fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz")


#@unittest.skip("need to use an example of RNA-SEQ sample")
@unittest.skip('skip')
class TestHTTPError(unittest.TestCase):
    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.src = array_express.ArrayExpress(cache_dirname=self.cache_dirname, download_backups=False, load_content=False)

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_error(self):
        src = self.src
        session = src.session
        ax = "E-MTAB-5530"
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}".format(ax))
        print("blue")


@unittest.skip('skip')
class TestMinervaProcessing(unittest.TestCase):

    #def setUp(self):
    #    self.cache_dirname = '{}/test_for_minerva'.format(os.path.dirname(os.path.realpath(__file__)))
    #    self.output_directory= '{}/test_for_minerva/output_directory'.format(os.path.dirname(os.path.realpath(__file__)))
    #    self.temporary_directory= '{}/test_for_minerva/temporary_directory'.format(os.path.dirname(os.path.realpath(__file__)))

    def setUp(self):
        #self.cache_dirname = tempfile.mkdtemp()
        self.cache_dirname = '{}/test_for_minerva'.format(os.path.dirname(os.path.realpath(__file__)))
        self.output_directory= '{}/test_for_minerva/output_directory'.format(os.path.dirname(os.path.realpath(__file__)))
        self.temporary_directory= '{}/test_for_minerva/temporary_directory'.format(os.path.dirname(os.path.realpath(__file__)))
        self.src = array_express.ArrayExpress(cache_dirname=self.cache_dirname, download_backups=False, load_content=False)

    def tearDown(self):
        pass
        #os.unlink('{}/ArrayExpress.sqlite'.format(self.cache_dirname))

    def test_minerva(self):
        a = self.src
        #a.load_content("https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/E-MTAB-6099")
        session = a.session
        print("beginning processing")
        exp = session.query(array_express.Experiment).filter_by(id="E-MTAB-6099").first()
        samples = exp.samples
        core.get_processed_data_samples(samples, self.output_directory, self.temporary_directory)
