""" Test the ArrayExpress database
:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-08-16
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.data_source import array_express
from kinetic_datanator.data_source.array_express_tools import ensembl_tools
from kinetic_datanator.data_source.process_rna_seq import core
from kinetic_datanator.data_source.process_rna_seq import download_cdna
from six.moves.urllib.request import urlretrieve
import datetime
import os
import pandas
import shutil
import tempfile
import unittest

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

    @unittest.skip("need to use an example of RNA-SEQ sample")
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

@unittest.skip("too long")
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
        self.assertEqual(len(a_sample.fastq_urls), 1)  # this sample has one entry
        a_sample = session.query(array_express.Sample).filter_by(experiment=exp).filter_by(name="GSM1409720 1").first()
        self.assertEqual(len(a_sample.fastq_urls), 1)  # this sample has two entries, but the url is the same
        self.assertEqual(a_sample.fastq_urls[0].url, "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/000/SRR1377310/SRR1377310.fastq.gz")

        self.src.load_experiment_samples(array_express.Experiment(id="E-GEOD-45474"))
        exp = session.query(array_express.Experiment).filter_by(id="E-GEOD-45474").first()
        self.assertEqual(exp.read_type, "single")
        self.assertEqual(len(exp.samples), 33)
        a_sample = session.query(array_express.Sample).filter_by(experiment=exp).filter_by(name="GSM1105145 1").first()
        self.assertEqual(a_sample.read_type, "single")
        self.assertEqual(len(a_sample.fastq_urls), 2)  # this sample has two entries, but the urls are different, so they are both added

    def test_without_fastq_files(self):
        session = self.src.session
        self.src.load_experiment_samples(array_express.Experiment(id="E-GEOD-32190"))
        exp = session.query(array_express.Experiment).filter_by(id="E-GEOD-32190").first()
        self.assertFalse(exp.has_fastq_files)
        # even though the fastq files aren't available, it is still important to record the read type for the user
        self.assertEqual(exp.read_type, "paired")

@unittest.skip("too long")
class TestProcessRnaseq(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.src = array_express.ArrayExpress(cache_dirname=self.cache_dirname, download_backups=False, load_content=False)

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)


    def test_get_processed_data_prokaryote_paired(self):
        src = self.src
        session = src.session
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/E-MTAB-3381")
        exp = session.query(array_express.Experiment).filter_by(id="E-MTAB-3381").first()
        core.get_processed_data(exp, self.cache_dirname)

    def test_total_eukaryote_single_end(self):
        src = self.src
        session = src.session
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/E-MTAB-3252")
        exp = session.query(array_express.Experiment).filter_by(id="E-MTAB-3252").first()
        core.get_processed_data(exp, self.cache_dirname)
        self.assertTrue(os.path.isfile("{}/E-MTAB-3252/Sample 1/Sample 1_abundances_binary".format(self.cache_dirname)))
        data = pandas.read_pickle("""{}/E-MTAB-3252/Sample 1/Sample 1_abundances_binary""".format(self.cache_dirname))
        self.assertEqual(data.loc["CCP42723", "percent total"], 0.00044199100205569983)

@unittest.skip("too long")
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
        sample = exp.samples[0]
        download_cdna.run(sample, self.cache_dirname)
        self.assertTrue(os.path.isfile('{}/CDNA_FILES/burkholderia_cenocepacia_j2315.cdna.all.fa.gz'.format(self.cache_dirname)))
        self.assertTrue(os.path.isfile('{}/kallisto_index_files/burkholderia_cenocepacia_j2315.idx'.format(self.cache_dirname)))

    def test_get_cdna_plants(self):
        # Arabidopsis thaliana
        src = self.src
        session = src.session
        ax = "E-GEOD-61293"
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}".format(ax))
        exp = session.query(array_express.Experiment).filter_by(id=ax).first()
        sample = exp.samples[0]
        download_cdna.run(sample, self.cache_dirname)
        self.assertTrue(os.path.isfile('{}/CDNA_FILES/arabidopsis_thaliana.cdna.all.fa.gz'.format(self.cache_dirname)))
        self.assertTrue(os.path.isfile('{}/kallisto_index_files/arabidopsis_thaliana.idx'.format(self.cache_dirname)))

    def test_get_cdna_eukaryotes(self):
        # Arabidopsis thaliana
        src = self.src
        session = src.session
        ax = "E-MTAB-5149"
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}".format(ax))
        exp = session.query(array_express.Experiment).filter_by(id=ax).first()
        sample = exp.samples[0]
        download_cdna.run(sample, self.cache_dirname)
        self.assertTrue(os.path.isfile('{}/CDNA_FILES/saccharomyces_cerevisiae.cdna.all.fa.gz'.format(self.cache_dirname)))
        self.assertTrue(os.path.isfile('{}/kallisto_index_files/saccharomyces_cerevisiae.idx'.format(self.cache_dirname)))


@unittest.skip("too long")
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
        ensembl_info = ensembl_tools.get_ensembl_info(sample)
        self.assertEqual(ensembl_info.organism_strain, "escherichia_coli_k_12_mg1655")
        self.assertTrue(ensembl_info.full_strain_specificity)
        self.assertEqual(ensembl_info.download_url, "ftp://ftp.ensemblgenomes.org/pub/bacteria/current/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cdna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa.gz")
        #"escherichia_coli_str_k_12_substr_mg1655"

    def test_prokaryote_2(self):
        src = self.src
        session = src.session
        ax = "E-MTAB-3396"
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}".format(ax))
        exp = session.query(array_express.Experiment).filter_by(id=ax).first()
        sample = exp.samples[0]
        ensembl_info = ensembl_tools.get_ensembl_info(sample)
        self.assertEqual(ensembl_info.organism_strain, "escherichia_coli_k_12")
        self.assertFalse(ensembl_info.full_strain_specificity)

    def test_eukaryote(self):
        src = self.src
        session = src.session
        ax = "E-MTAB-5149"
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}".format(ax))
        exp = session.query(array_express.Experiment).filter_by(id=ax).first()
        sample = exp.samples[0]
        ensembl_info = ensembl_tools.get_ensembl_info(sample)
        self.assertEqual(ensembl_info.organism_strain, "saccharomyces_cerevisiae")
        self.assertTrue(ensembl_info.full_strain_specificity)
        self.assertEqual(ensembl_info.download_url,
                         "ftp://ftp.ensembl.org/pub/current_fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz")

    def test_plant(self):
        src = self.src
        session = src.session
        ax = "E-GEOD-61293"
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}".format(ax))
        exp = session.query(array_express.Experiment).filter_by(id=ax).first()
        sample = exp.samples[0]
        ensembl_info = ensembl_tools.get_ensembl_info(sample)
        self.assertEqual(sample.ensembl_info[0].organism_strain, "arabidopsis_thaliana")
        self.assertTrue(sample.full_strain_specificity)
        self.assertEqual(sample.ensembl_info[
                         0].url, "ftp://ftp.ensemblgenomes.org/pub/current/plants/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz")


    def test_substrain_3(self):
        # in this sample, the strain ("avian pathogen") is not recognized, and
        # therefore only the organism (escherichia coli k12) is downloaded
        src = self.src
        session = src.session
        ax = "E-MTAB-3396"
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}".format(ax))
        exp = session.query(array_express.Experiment).filter_by(id=ax).first()
        sample = exp.samples[0]
        ensembl_info = ensembl_tools.get_ensembl_info(sample)
        self.assertEqual(sample.ensembl_info[0].organism_strain, "escherichia_coli_k_12")


    def test_substrain_4(self):
        src = self.src
        session = src.session
        # ax = "E-GEOD-65711" #This one is a problem because strain is listed as "K-12 MG1655" instead of "str. K-12 subtr MG1655"
        # E-MTAB-3699
        ax = "E-GEOD-65711"
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}".format(ax))
        exp = session.query(array_express.Experiment).filter_by(id=ax).first()
        sample = exp.samples[0]
        ensembl_info = ensembl_tools.get_ensembl_info(sample)
        self.assertEqual(sample.ensembl_info[0].organism_strain, "escherichia_coli_k_12_mg1655")

    def test_substrain_5(self):
        # test baccillus subtilis
        src = self.src
        session = src.session
        ax = "E-GEOD-35641"
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}".format(ax))
        exp = session.query(array_express.Experiment).filter_by(id=ax).first()
        sample = []
        for sam in exp.samples:
            if sam.name == 'GSM872395 1':
                sample = sam
        ensembl_info = ensembl_tools.get_ensembl_info(sample)
        self.assertEqual(sample.ensembl_info[0].organism_strain, "bacillus_subtilis_subtilis_168")


    def test_substrain_6(self):
        # Staphylococcus aureus subsp. aureus USA300
        src = self.src
        session = src.session
        ax = "E-MTAB-5666"
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}".format(ax))
        exp = session.query(array_express.Experiment).filter_by(id=ax).first()
        sample = exp.samples[0]
        ensembl_info = ensembl_tools.get_ensembl_info(sample)
        self.assertEqual(sample.ensembl_info[0].organism_strain, "staphylococcus_aureus_aureus")


    def test_substrain_8(self):
        # Listeria monocytogenes 10403S
        src = self.src
        session = src.session
        ax = "E-GEOD-70265"
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}".format(ax))
        exp = session.query(array_express.Experiment).filter_by(id=ax).first()
        sample = exp.samples[0]
        ensembl_info = ensembl_tools.get_ensembl_info(sample)
        self.assertEqual(sample.ensembl_info[0].organism_strain, "listeria_monocytogenes_10403s")

    def test_substrain_9(self):
        # Desulfovibrio vulgaris str. Hildenborough
        src = self.src
        session = src.session
        ax = "E-GEOD-79022"
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}".format(ax))
        exp = session.query(array_express.Experiment).filter_by(id=ax).first()
        sample = exp.samples[0]
        ensembl_info = ensembl_tools.get_ensembl_info(sample)
        self.assertEqual(sample.ensembl_info[0].organism_strain, "desulfovibrio_vulgaris_hildenborough")


    def test_substrain_11(self):
        # Listeria monocytogenes EGD
        src = self.src
        session = src.session
        ax = "E-GEOD-65558"
        src.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}".format(ax))
        exp = session.query(array_express.Experiment).filter_by(id=ax).first()
        sample = exp.samples[0]
        ensembl_info = ensembl_tools.get_ensembl_info(sample)
        self.assertEqual(sample.ensembl_info[0].organism_strain, "listeria_monocytogenes_egd")
