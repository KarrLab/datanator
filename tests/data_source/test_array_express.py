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


class QuickTest(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.src = array_express.ArrayExpress(cache_dirname=self.cache_dirname, download_backup=False, load_content=False)

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_load_experiment_metadata(self):
        src = self.src
        session = src.session

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
        self.assertEqual(sorted(experiment.data_formats, key=lambda format: format.name), [
            session.query(array_express.DataFormat).filter_by(name='normalization').first(),
            session.query(array_express.DataFormat).filter_by(name='processedData').first(),
        ])
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

        src.load_experiment_metadata(2001, 2001)

        # E-GEOD-10
        experiment = session.query(array_express.Experiment).filter_by(id='E-GEOD-6').first()
        src.load_experiment_protocols(experiment)
        self.assertEqual(len(experiment.protocols), 2)
        self.assertEqual(set([protocol.protocol_accession for protocol in experiment.protocols]), set(['P-GSE6-1', 'P-GSE6-2', ]))

        q = session.query(array_express.Protocol)
        protocol = session.query(array_express.Protocol).filter_by(protocol_accession='P-GSE6-1').first()
        start_text = "ID_REF = \nCH1B_MEAN = mean channel 1 background\n"
        end_text = "\nUNF_VALUE = aka LOG_RAT2N_MEAN, or log2 of ratio of CH2DN_MEAN and CH1DN_MEAN"
        self.assertEqual(protocol.text[0:len(start_text)], start_text)
        self.assertEqual(protocol.text[-len(end_text):], end_text)

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
