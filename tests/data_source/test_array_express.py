""" Test the ArrayExpress database
:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-08-16
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.data_source import array_express
import datetime
import shutil
import tempfile
import unittest
import os.path
import zipfile
import datetime
import dateutil.parser


class QuickTest(unittest.TestCase):


    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.src = array_express.ArrayExpress(cache_dirname=self.cache_dirname, download_backup=False, load_content=False)

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
        # print "yello"
        self.assertEqual(sorted([df.name for df in experiment.data_formats]), sorted([df.name for df in [
            session.query(array_express.DataFormat).filter_by(name='normalization').first(),
            session.query(array_express.DataFormat).filter_by(name='processedData').first(),
        ]]))

        #self.assertEqual(sorted(experiment.data_formats, key=lambda format: format.name), [
        #    session.query(array_express.DataFormat).filter_by(name='normalization').first(),
        #    session.query(array_express.DataFormat).filter_by(name='processedData').first(),
        #])
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

        self.assertEqual(len(experiment.samples), 20)

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



    def test_load_processed_data(self):
        test_experiment = array_express.Experiment(id = 'E-MTAB-4549', organisms=[array_express.Organism(name="Mus musculus")])
        src = self.src
        src.load_processed_data(test_experiment)
        self.assertTrue(os.path.isfile("{}/array_express_processed_data/E-MTAB-4549.processed.1.zip".format(src.cache_dirname)))
        zip_ref = zipfile.ZipFile("{}/array_express_processed_data/E-MTAB-4549.processed.1.zip".format(src.cache_dirname), 'r')
        zip_ref.extractall(src.cache_dirname)
        zip_ref.close()
        filename = "{}/raw_countsHSC_bulk.txt".format(src.cache_dirname)
        file = open(filename, "r")
        self.assertEqual(file.read()[:9], 'dHSC dHSC')





    def test_load_content(self):
        src = self.src
        session = src.session

        src.load_content(start_year=2001, end_year=2001)

        experiment = session.query(array_express.Experiment).filter_by(id='E-GEOD-10').first()
        self.assertEqual(set([sample.name for sample in experiment.samples]), set(['GSM571 1', 'GSM572 1', 'GSM573 1', 'GSM574 1']))

        sample = next(sample for sample in experiment.samples if sample.name == 'GSM574 1')
        self.assertEqual(sample.name, 'GSM574 1')
        self.assertEqual(sample.extracts, [session.query(array_express.Extract).filter_by(name='GSM574 extract 1').first()])
        self.assertEqual(sample.assay, 'GSM574')
        self.assertEqual(sample.characteristics, [session.query(
            array_express.Characteristic).filter_by(category='Organism', value='Homo sapiens').first()])
        self.assertEqual(sample.variables, [])


class LongTest(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.src = array_express.ArrayExpress(cache_dirname=self.cache_dirname, download_backup=False, load_content=False)

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_load_content(self):
        src = self.src
        session = src.session

        src.load_content(start_year=2001, end_year=2002)

        # E-GEOD-10
        experiment = session.query(array_express.Experiment).filter_by(id='E-GEOD-10').first()
        self.assertEqual(set([sample.name for sample in experiment.samples]), set(['GSM571 1', 'GSM572 1', 'GSM573 1', 'GSM574 1']))

        sample = next(sample for sample in experiment.samples if sample.name == 'GSM574 1')
        self.assertEqual(sample.name, 'GSM574 1')
        self.assertEqual(sample.extracts, [session.query(array_express.Extract).filter_by(name='GSM574 extract 1').first()])
        self.assertEqual(sample.assay, 'GSM574')
        self.assertEqual(sample.characteristics, [session.query(
            array_express.Characteristic).filter_by(category='Organism', value='Homo sapiens').first()])
        self.assertEqual(sample.variables, [])

        # E-SNGR-6
        experiment = session.query(array_express.Experiment).filter_by(id='E-SNGR-6').first()
        self.assertEqual(len(experiment.samples), 20)

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


class TestDownloadProcessedData(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.src = array_express.ArrayExpress(cache_dirname=self.cache_dirname, download_backup=False, load_content=False, verbose=True)

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_download_processed_data(self):
        src = self.src

        response = src.requests_session.get(src.ENDPOINT_DOMAINS['array_express'] + '?date=[2008-12-01+2008-12-20]')
        response.raise_for_status()
        for expt_json in response.json()['experiments']['experiment']:
            id = expt_json['accession']
            experiment = src.get_or_create_object(array_express.Experiment, id=id)
            if isinstance(expt_json['name'], list):
                experiment.name = expt_json['name'][0]
                experiment.name_2 = expt_json['name'][1]
            else:
                experiment.name = expt_json['name']
            if 'organism' in expt_json:
                for organism_name in expt_json['organism']:
                    experiment.organisms.append(src.get_or_create_object(array_express.Organism, name=organism_name))
            if 'description' in expt_json:
                experiment.description = expt_json['description'][0]['text']
            if 'experimenttype' in expt_json:
                entries = expt_json['experimenttype']
                for entry in entries:
                    experiment.types.append(src.get_or_create_object(array_express.ExperimentType, name=entry))
            if 'experimentdesign' in expt_json:
                entries = expt_json['experimentdesign']
                for entry in entries:
                    experiment.designs.append(src.get_or_create_object(array_express.ExperimentDesign, name=entry))
            if 'bioassaydatagroup' in expt_json:
                for entry in expt_json['bioassaydatagroup']:
                    experiment.data_formats.append(src.get_or_create_object(array_express.DataFormat, name=entry['dataformat']))
            if 'submissiondate' in expt_json:
                experiment.submission_date = dateutil.parser.parse(expt_json['submissiondate']).date()
            if 'releasedate' in expt_json:
                experiment.release_date = dateutil.parser.parse(expt_json['releasedate']).date()
            src.session.add(experiment)

        for i_experiment, experiment in enumerate(src.session.query(array_express.Experiment).all()):

            if ('processedData' in [d.name for d in experiment.data_formats]) and ("RNA-seq of coding RNA" in [d.name for d in experiment.types]):
                src.load_processed_data(experiment)

        self.assertTrue(os.path.isfile("{}/array_express_processed_data/E-GEOD-13518.processed.1.zip".format(src.cache_dirname)))

        #the following should not download because there is not processed file
        self.assertFalse(os.path.isfile("{}/array_express_processed_data/E-MTAB-71.processed.1.zip".format(src.cache_dirname)))

        #the following should not download becuase the experiment is not rna seq of coding rna
        self.assertFalse(os.path.isfile("{}/array_express_processed_data/E-MEXP-1386.processed.1.zip".format(src.cache_dirname)))


class TestIncorporateProcessedData(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.src = array_express.ArrayExpress(cache_dirname=self.cache_dirname, download_backup=False, load_content=False, verbose=True)

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_download_processed_data(self):
        src = self.src
        session = src.session

        #src.load_experiment_metadata(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments?date=[2016-02-02+2016-02-02]")
        src.load_experiment_metadata(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments?accession=E-MTAB-4063")
        experiments = session.query(array_express.Experiment)

        for i_experiment, experiment in enumerate(session.query(array_express.Experiment).all()):
            src.load_experiment_samples(experiment)
            src.load_experiment_protocols(experiment)
            single_dimensional_processed_RNA_seq_data = False
            # print [d.name for d in experiment.types]
            if ("RNA-seq of coding RNA" in [d.name for d in experiment.types]):
                for df in experiment.data_formats:
                    # print df.name
                    # print df.bio_assay_data_cubes
                    # print ""
                    if df.name == 'processedData' and df.bio_assay_data_cubes == 1:
                        single_dimensional_processed_RNA_seq_data = True
                        break

            #if ('processedData' in [d.name for d in experiment.data_formats]) and ("RNA-seq of coding RNA" in [d.name for d in experiment.types]):
            if single_dimensional_processed_RNA_seq_data:
                src.load_processed_data(experiment)

        self.assertTrue(os.path.isfile("{}/array_express_processed_data/E-MTAB-4063.processed.1.zip".format(src.cache_dirname)))
        self.assertTrue(os.path.isfile("{}/array_express_processed_data/E-MTAB-4063/NJMU24_normalized_counts.txt".format(src.cache_dirname)))
        experiment = session.query(array_express.Experiment).filter_by(id='E-MTAB-4063').first()
        self.assertTrue(experiment.whole_genome)
        self.assertTrue(len(experiment.genes), 57820)
        some_gene = session.query(array_express.Gene).filter_by(name='ENSG00000001084.6').first()
        self.assertEqual(some_gene.experiments[0].id, 'E-MTAB-4063')

    def test_download_processed_data_long(self):
        src = self.src
        session = src.session

        src.load_experiment_metadata(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments?date=[2016-02-02+2016-02-04]")
        #src.load_experiment_metadata(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments?accession=E-MTAB-4063")
        experiments = session.query(array_express.Experiment)

        for i_experiment, experiment in enumerate(session.query(array_express.Experiment).all()):
            # print experiment.id
            src.load_experiment_samples(experiment)
            src.load_experiment_protocols(experiment)
            single_dimensional_processed_RNA_seq_data = False
            if ("RNA-seq of coding RNA" in [d.name for d in experiment.types]):
                for df in experiment.data_formats:
                    if df.name == 'processedData' and df.bio_assay_data_cubes == 1:
                        single_dimensional_processed_RNA_seq_data = True
                        break

            #if ('processedData' in [d.name for d in experiment.data_formats]) and ("RNA-seq of coding RNA" in [d.name for d in experiment.types]):
            if single_dimensional_processed_RNA_seq_data:
                src.load_processed_data(experiment)

        self.assertEqual(len(session.query(array_express.Gene).all()), 97481)
        incorporated_exps = list(set([g.experiment_id for g in session.query(array_express.Gene).all()]))
        self.assertEqual(incorporated_exps, list(set([u'E-MTAB-3965', u'E-MTAB-4063'])))
        total_exps = set([exp.id for exp in session.query(array_express.Experiment).all()])
        unincorporated_exps = list(total_exps - set(incorporated_exps))
        self.assertEqual(set(unincorporated_exps), set([u'E-MTAB-4133',
            u'E-GEOD-77512',
            u'E-MTAB-3941', #unincorporated because no processed data
            u'E-GEOD-77428',
            u'E-MTAB-3680', #unincorporated because no processed data
            u'E-GEOD-69220',
            u'E-GEOD-77479',
            u'E-GEOD-66821',
            u'E-GEOD-72210',
            u'E-GEOD-65859',
            u'E-MTAB-4231', #unincorporated because not an RNA seq experiment
            u'E-MTAB-3301', #unincorporated because has too many bio_assay_data_cubes
            u'E-GEOD-71319',
            u'E-GEOD-77224'
            ]))
