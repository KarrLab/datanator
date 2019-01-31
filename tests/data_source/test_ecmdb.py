""" Tests of ecmdb
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-05
:Copyright: 2017, Karr Lab
:License: MIT
"""

from datanator.data_source import ecmdb
import datetime
import dateutil
import os
import shutil
import tempfile
import unittest



class TestEcmdbFromRemote(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        src = ecmdb.Ecmdb(cache_dirname=cls.cache_dirname, download_backups=False, load_content=False, verbose=True, max_entries=12)
        src.load_content()
        src.session.close()
        src.engine.dispose()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    def setUp(self):
        self.src = ecmdb.Ecmdb(cache_dirname=self.cache_dirname, download_backups=False, load_content=False, verbose=True, max_entries=12)

    def tearDown(self):
        src = self.src
        src.session.close()
        src.engine.dispose()

    def test_download_from_remote(self):
        src = self.src
        session = src.session
        q = session.query(ecmdb.Compound)

        # check compound ids
        self.assertEqual(q.count(), src.max_entries)
        self.assertEqual([c.id for c in q.all()][0:5], [
            'M2MDB000001',
            'M2MDB000002',
            'M2MDB000003',
            'M2MDB000004',
            'M2MDB000005',
        ])

        # check one compound in depth
        compound = q.first()
        self.assertEqual(compound.name, '2-Ketobutyric acid')

        self.assertEqual([s.name for s in compound.synonyms[0:4]], [
            '2-oxobutanoic acid',
            '&alpha;-ketobutyrate',
            '&alpha;-ketobutyric acid',
            '&alpha;-oxobutyrate',
        ])
        self.assertTrue(compound.description.startswith('2-Ketobutyric acid (alpha ketobutyric acid)'))
        self.assertEqual(compound.structure, 'InChI=1S/C4H6O3/c1-2-3(5)4(6)7/h2H2,1H3,(H,6,7)')
        self.assertEqual(compound._structure_formula_connectivity, 'C4H6O3/c1-2-3(5)4(6)7')
        self.assertEqual([c.name for c in compound.compartments], ['Cytosol'])

        self.assertEqual(compound.concentrations, [])

        self.assertEqual(set([(xr.namespace, xr.id) for xr in compound.cross_references]), set([
            ('biocyc', '2-OXOBUTANOATE'),
            ('cas', '600-18-0'),
            ('hmdb', 'HMDB00005'),
            ('chemspider', '57'),
            ('ligandexpo', '2KT'),
            ('wikipedia.en', 'Alpha-ketobutyric_acid'),
            ('kegg.compound', 'C00109'),
            ('chebi', 'CHEBI:16763'),
            ('pubchem.compound', '58'),
        ]))
        self.assertEqual(compound.comment, None)

        self.assertEqual(compound.created, dateutil.parser.parse('2011-05-29 15:47:17 UTC').replace(tzinfo=None))
        #self.assertEqual(compound.updated, dateutil.parser.parse('2015-06-03 15:00:41 -0600').replace(tzinfo=None))
        self.assertLess((datetime.datetime.utcnow() - compound.downloaded).total_seconds(), 3000)

        # compound with multiple compartments
        compound = session.query(ecmdb.Compound).filter_by(id='M2MDB000002').first()
        self.assertEqual([c.name for c in compound.compartments], ['Cytosol', 'Extra-organism', 'Periplasm'])

        # compound with one concentrations
        compound = session.query(ecmdb.Compound).filter_by(id='M2MDB000012').first()
        self.assertEqual(len(compound.concentrations), 1)
        self.assertEqual(compound.concentrations[0].value, 658.0)
        self.assertEqual(compound.concentrations[0].error, 25.)
        self.assertEqual(compound.concentrations[0].strain, 'BL21 DE3')
        self.assertEqual(compound.concentrations[0].growth_status, 'Stationary phase cultures (overnight culture)')
        self.assertEqual(compound.concentrations[0].media, 'Luria-Bertani (LB) media')
        self.assertEqual(compound.concentrations[0].temperature, 37)
        self.assertEqual(compound.concentrations[0].growth_system, 'Shake flask')
        self.assertEqual(len(compound.concentrations[0].references), 1)
        self.assertEqual(compound.concentrations[0].references[0].namespace, 'pubmed')
        self.assertEqual(compound.concentrations[0].references[0].id, '17535911')

        # compound with multiple concentrations
        compound = session.query(ecmdb.Compound).filter_by(id='M2MDB000010').first()
        self.assertEqual(len(compound.concentrations), 5)
        self.assertEqual(compound.concentrations[0].value, 1.47)
        self.assertEqual(compound.concentrations[0].error, 0.0)
        self.assertEqual(compound.concentrations[0].strain, 'K12 NCM3722')
        self.assertEqual(compound.concentrations[0].growth_status, 'Mid-Log Phase')
        self.assertTrue(compound.concentrations[0].media.startswith('Gutnick minimal complete medium'))
        self.assertEqual(compound.concentrations[0].temperature, 37)
        self.assertEqual(compound.concentrations[0].growth_system, 'Shake flask and filter culture')
        self.assertEqual(len(compound.concentrations[0].references), 1)
        self.assertEqual(compound.concentrations[0].references[0].namespace, 'pubmed')
        self.assertEqual(compound.concentrations[0].references[0].id, '19561621')

        # compound with comment
        # compound = session.query(ecmdb.Compound).filter_by(id='M2MDB000734').first()

    def test_cascade_concentrations(self):
        src = self.src
        session = src.session

        n_compound = session.query(ecmdb.Compound).count()
        n_concentration = session.query(ecmdb.Concentration).count()

        compound = session.query(ecmdb.Compound).filter_by(id='M2MDB000010').first()
        n_concentration_compound = len(compound.concentrations)
        session.delete(compound)

        self.assertEqual(session.query(ecmdb.Compound).count(), n_compound - 1)
        self.assertEqual(session.query(ecmdb.Concentration).count(), n_concentration - n_concentration_compound)

    def test_cascade_compounds_1(self):
        src = self.src
        session = src.session

        n_compound = session.query(ecmdb.Compound).count()
        n_concentration = session.query(ecmdb.Concentration).count()

        compound = session.query(ecmdb.Compound).filter_by(id='M2MDB000010').first()
        n_concentration_compound = len(compound.concentrations)
        compound.concentrations[0].compound = None

        self.assertEqual(len(compound.concentrations), n_concentration_compound - 1)
        self.assertEqual(session.query(ecmdb.Compound).count(), n_compound)
        self.assertEqual(session.query(ecmdb.Concentration).count(), n_concentration - 1)

    def test_cascade_compounds_2(self):
        src = self.src
        session = src.session

        n_compound = session.query(ecmdb.Compound).count()
        n_concentration = session.query(ecmdb.Concentration).count()

        compound = session.query(ecmdb.Compound).filter_by(id='M2MDB000010').first()
        n_concentration_compound = len(compound.concentrations)
        concentration = compound.concentrations[0]
        concentration.compound = None
        session.delete(concentration)

        self.assertEqual(len(compound.concentrations), n_concentration_compound - 1)
        self.assertEqual(session.query(ecmdb.Compound).count(), n_compound)
        self.assertEqual(session.query(ecmdb.Concentration).count(), n_concentration - 1)

    def test_cascade_compounds_3(self):
        src = self.src
        session = src.session

        n_compound = session.query(ecmdb.Compound).count()
        n_concentration = session.query(ecmdb.Concentration).count()

        concentration = session.query(ecmdb.Concentration).first()
        session.delete(concentration)

        self.assertEqual(session.query(ecmdb.Compound).count(), n_compound)
        self.assertEqual(session.query(ecmdb.Concentration).count(), n_concentration - 1)


# class TestEcmdbFromBlank(unittest.TestCase):

#     def setUp(self):
#         self.cache_dirname = tempfile.mkdtemp()

#     def tearDown(self):
#         shutil.rmtree(self.cache_dirname)

#     def test_commit_intermediate_results(self):
#         src = ecmdb.Ecmdb(cache_dirname=self.cache_dirname, download_backups=False, load_content=False, verbose=True, max_entries=5,
#                           commit_intermediate_results=True)
#         src.load_content()

#     @unittest.skip('Skip because this is a long test')
#     def test_download_all(self):
#         src = ecmdb.Ecmdb(download_backups=False, load_content=True, clear_content=True, clear_requests_cache=False, verbose=True)
#         src.upload_backups()

#         session = src.session
#         self.assertGreater(session.query(ecmdb.Compound).count(), 3500)

#     def test_download_from_backup(self):
#         src = ecmdb.Ecmdb(cache_dirname=self.cache_dirname)
#         session = src.session

#         self.assertGreater(session.query(ecmdb.Compound).count(), 3500)

# class TestEcmdbFromCache(unittest.TestCase):
#     """
#     Quick test to ensure Ecmdb on Karr Lab Server is correct

#     """

#     @classmethod
#     def setUpClass(cls):
#         cls.cache_dirname = tempfile.mkdtemp()
#         cls.src = ecmdb.Ecmdb(cache_dirname=cls.cache_dirname, download_backups=True, load_content=False, verbose=True)

#     @classmethod
#     def tearDownClass(cls):
#         shutil.rmtree(cls.cache_dirname)

#     def test_proper_loading(self):
#         session = self.src.session

#         q = session.query(ecmdb.Compound).filter_by(id = 'M2MDB000090')
#         self.assertEqual(q.count(), 1)

#         self.assertEqual(q.first().name, 'Orotidylic acid')
#         self.assertEqual(q.first().structure, 'InChI=1S/C10H13N2O11P/c13-5-1-3(9(16)17)12(10(18)11-5)8-7(15)6(14)4(23-8)2-22-24(19,20)21/h1,4,6-8,14-15H,2H2,(H,16,17)(H,11,13,18)(H2,19,20,21)/t4-,6-,7-,8-/m1/s1')

#         q = session.query(ecmdb.Concentration).get(123)

#         self.assertEqual(q.value, 207.0)
#         self.assertEqual(q.error, 0.0)
#         self.assertEqual(q.growth_status, 'Stationary Phase, glucose limited')
