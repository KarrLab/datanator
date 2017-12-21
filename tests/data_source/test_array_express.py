""" Test the ArrayExpress database
:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-08-16
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.data_source import array_express
import shutil
import tempfile
import unittest
import datetime


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

        src.load_experiment_metadata(2001, 2002)
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

    def test_load_experiment_samples(self):
        src = self.src
        session = src.session

        src.load_experiment_metadata(2001, 2002)

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

        src.load_experiment_metadata(2001, 2002)

        # E-GEOD-10
        experiment = session.query(array_express.Experiment).filter_by(id='E-GEOD-6').first()
        src.load_experiment_protocols(experiment)
        self.assertEqual(len(experiment.protocols), 2)
        self.assertEqual(set([protocol.protocol_accession for protocol in experiment.protocols]), set(['P-GSE6-1', 'P-GSE6-2',]))

        q = session.query(array_express.Protocol)
        protocol = session.query(array_express.Protocol).filter_by(protocol_accession='P-GSE6-1').first()
        text = "ID_REF = \nCH1B_MEAN = mean channel 1 background\nCH1B_MEDIAN = median channel 1 background\nCH1D_MEAN = difference between CH1I_MEAN and CH1B_MEAN\nCH2I_MEAN = mean channel 2 signal\nCH2D_MEAN = difference of CH2I_MEAN and CH2B_MEAN\nCH2B_MEAN = mean of channel 2 background\nCH2B_MEDIAN = median of channel 2 background\nCH2BN_MEDIAN = normalized CH2B_MEDIAN\nCH2DN_MEAN = normalized CH2D_MEAN\nCH2IN_MEAN = normalized CH2I_MEAN\nCORR = simple correlation coefficient of channel 1 and channel 2 pixels\nFLAG = user defined flag (default=0)\nVALUE = same as UNF_VALUE but with flagged values removed\nPIX_RAT2_MEDIAN =  \nPERGTBCH1I_1SD = standard deviation of fraction of pixels greater than CH1B_MEAN\nPERGTBCH2I_1SD = standard deviation of fraction of pixels greater than CH2B_MEAN\nRAT1_MEAN = ratio of CH1D_MEAN to CH2D_MEAN\nRAT1N_MEAN = ratio of CH1DN_MEAN to CH2DN_MEAN\nRAT2_MEAN = ratio of CH2D_MEAN to CH1D_MEAN\nRAT2N_MEAN = ratio of CH2DN_MEAN and CH1DN_MEAN\nREGR = slope of the simple regression line of channel 2 to channel 1 pixels\nTOT_BPIX = total number of background pixels\nTOT_SPIX = total number of signal (spot) pixels\nTOP = vertical, short-axis spot ellipse minimum pixel coordinate\nBOT = vertical, short-axis spot ellipse maximum pixel coordinate\nLEFT = horizontal, long-axis spot ellipse minimum pixel coordinate\nRIGHT = horizontal, long-axis spot ellipse maximum pixel coordinate\nUNF_VALUE = aka LOG_RAT2N_MEAN, or log2 of ratio of CH2DN_MEAN and CH1DN_MEAN"
        self.assertEqual(protocol.text, text)

        test_experiment = array_express.Experiment(id = 'E-MTAB-5281')
        src.load_experiment_protocols(test_experiment)
        protocol = session.query(array_express.Protocol).filter_by(protocol_accession='P-MTAB-42502').first()
        self.assertEqual(protocol.protocol_type,'labelling')
        self.assertEqual(protocol.text,'Labelling was performed by using enzymatic attachment of nucleotides coupled to biotin.')
        self.assertEqual(protocol.performer,'Stephanie Boue')
        self.assertEqual(protocol.hardware,'Affymetrix GeneChip Scanner 3000 7G')
        self.assertEqual(protocol.software,'Affymetrix AGCC')

    def test_load_content(self):
        src = self.src
        session = src.session
        #print(len(session.query(array_express.Experiment).all()))
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
        a_sample = session.query(array_express.Sample).filter_by(experiment=exp).filter_by(name="GSM1015798 1").first()
        self.assertEqual(a_sample.read_type, "single")
        self.assertEqual(a_sample.fastq_urls[0].url, "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR580/SRR580556/SRR580556.fastq.gz")

    def test_paired_read(self):
        session = self.src.session
        self.src.load_experiment_samples(array_express.Experiment(id="E-MTAB-4262"))
        exp = session.query(array_express.Experiment).filter_by(id="E-MTAB-4262").first()
        self.assertTrue(exp.has_fastq_files)
        self.assertEqual(exp.read_type, "paired")
        self.assertEqual(len(exp.samples), 12)
        a_sample = session.query(array_express.Sample).filter_by(experiment=exp).filter_by(name="R1_0mM_24h").first()
        self.assertEqual(a_sample.read_type, "paired")
        self.assertEqual(list(set([
            "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR173/002/ERR1736192/ERR1736192_1.fastq.gz",
            "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR173/002/ERR1736192/ERR1736192_2.fastq.gz",
            ])), list(set([u.url for u in a_sample.fastq_urls])))


    def test_multiple_reads(self):
        session = self.src.session
        self.src.load_experiment_samples(array_express.Experiment(id="E-GEOD-58388"))
        exp = session.query(array_express.Experiment).filter_by(id="E-GEOD-58388").first()
        self.assertEqual(exp.read_type, "single")
        self.assertEqual(len(exp.samples), 12)
        a_sample = session.query(array_express.Sample).filter_by(experiment=exp).filter_by(name="GSM1409710 1").first()
        self.assertEqual(a_sample.read_type, "single")
        self.assertEqual(len(a_sample.fastq_urls), 1) #this sample has one entry
        a_sample = session.query(array_express.Sample).filter_by(experiment=exp).filter_by(name="GSM1409720 1").first()
        self.assertEqual(len(a_sample.fastq_urls), 1) #this sample has two entries, but the url is the same
        self.assertEqual(a_sample.fastq_urls[0].url, "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/000/SRR1377310/SRR1377310.fastq.gz")
        
        self.src.load_experiment_samples(array_express.Experiment(id="E-GEOD-45474"))
        exp = session.query(array_express.Experiment).filter_by(id="E-GEOD-45474").first()
        self.assertEqual(exp.read_type, "single")
        self.assertEqual(len(exp.samples), 33)
        a_sample = session.query(array_express.Sample).filter_by(experiment=exp).filter_by(name="GSM1105145 1").first()
        self.assertEqual(a_sample.read_type, "single")
        self.assertEqual(len(a_sample.fastq_urls), 2) #this sample has two entries, but the urls are different, so they are both added

    def test_without_fastq_files(self):
        session = self.src.session
        self.src.load_experiment_samples(array_express.Experiment(id="E-GEOD-32190"))
        exp = session.query(array_express.Experiment).filter_by(id="E-GEOD-32190").first()
        self.assertFalse(exp.has_fastq_files)
        self.assertEqual(exp.read_type, "paired") #even though the fastq files aren't available, it is still important to record the read type for the user
