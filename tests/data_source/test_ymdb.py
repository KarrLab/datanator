""" Tests of ecmdb
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-05
:Copyright: 2017, Karr Lab
:License: MIT
"""

from datanator.data_source import ymdb
import datetime
import dateutil
import os
import shutil
import tempfile
import unittest



class TestYmdbFromRemote(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        src = ymdb.Ymdb(cache_dirname=cls.cache_dirname, download_backups=False, load_content=False, verbose=True, max_entries=12)
        src.load_content()
        src.session.close()
        src.engine.dispose()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    def setUp(self):
        self.src = ymdb.Ymdb(cache_dirname=self.cache_dirname, download_backups=False, load_content=False, verbose=True, max_entries=12)

    def tearDown(self):
        src = self.src
        src.session.close()
        src.engine.dispose()

    def test_download_from_remote(self):
        src = self.src
        session = src.session
        q = session.query(ymdb.Compound)

        # check compound ids
        self.assertEqual(q.count(), src.max_entries)
        self.assertEqual([c.id for c in q.all()][0:5], [
            'YMDB000001',
            'YMDB000002',
            'YMDB000003',
            'YMDB000004',
            'YMDB000005',
        ])

        # check one compound in depth
        compound = q.first()
        self.assertEqual(compound.name, "5'-Phosphoribosyl-N-formylglycinamide")

        self.assertEqual([s.name for s in compound.synonyms[0:4]], [
            "5'-P-ribosyl-N-formylglycineamide",
            "5'-phosphoribosyl-formylglycinamide",
            "5'-phosphoribosyl-N-formylglycinamide",
            "5'-phosphoribosyl-N-formylglycineamide",
        ])
        self.assertTrue(compound.description.startswith("5'-Phosphoribosyl-N-formylglycineamide (FGAR or N-Formyl-GAR)"))
        self.assertEqual(compound.structure, 'InChI=1S/C8H15N2O9P/c11-3-9-1-5(12)10-8-7(14)6(13)4(19-8)2-18-20(15,16)17/h3-4,6-8,13-14H,1-2H2,(H,9,11)(H,10,12)(H2,15,16,17)/t4-,6-,7-,8-/m1/s1')
        self.assertEqual(compound._structure_formula_connectivity, 'C8H15N2O9P')
        self.assertEqual([c.name for c in compound.compartments], ['cytoplasm'])

        self.assertEqual(compound.concentrations, [])

        self.assertEqual(set([(xr.namespace, xr.id) for xr in compound.cross_references]), set([
            ('biocyc', '5-P-RIBOSYL-N-FORMYLGLYCINEAMIDE'),
            ('cas', '349-34-8'),
            ('chebi', 'CHEBI:18272'),
            ('chemspider', '115687'),
            ('hmdb', 'HMDB01308'),
            ('kegg.compound', 'C04376'),
            ('pubchem.compound', '151'),
        ]))
        self.assertEqual(compound.comment, None)

        self.assertEqual(compound.created, dateutil.parser.parse('2011-05-29 15:47:17 UTC').replace(tzinfo=None))
        #self.assertEqual(compound.updated, dateutil.parser.parse('2015-06-03 15:00:41 -0600').replace(tzinfo=None))
        self.assertLess((datetime.datetime.utcnow() - compound.downloaded).total_seconds(), 3000)

        # compound with multiple compartments
        compound = session.query(ymdb.Compound).filter_by(id='YMDB000002').first()
        self.assertEqual([c.name for c in compound.compartments], ['extracellular', 'nucleus', 'vacuole', 'cytoplasm'])

        # compound with multiple concentrations
        compound = session.query(ymdb.Compound).filter_by(id='YMDB000022').first()
        self.assertEqual(len(compound.concentrations), 2)
        self.assertEqual(compound.concentrations[0].value, 9080.0)
        self.assertEqual(compound.concentrations[0].error, 0.0)
        self.assertEqual(compound.concentrations[0].strain, '')
        self.assertEqual(compound.concentrations[0].growth_status, 'Stationary phase cultures (overnight culture)')
        self.assertEqual(compound.concentrations[0].media.startswith('glucose~(140 g/L)'))
        self.assertEqual(len(compound.concentrations[0].references), 1)
        self.assertEqual(compound.concentrations[0].references[0].namespace, 'pubmed')
        self.assertEqual(compound.concentrations[0].references[0].id, '18609643')

        # # compound with multiple concentrations
        # compound = session.query(ecmdb.Compound).filter_by(id='M2MDB000010').first()
        # self.assertEqual(len(compound.concentrations), 5)
        # self.assertEqual(compound.concentrations[0].value, 1.47)
        # self.assertEqual(compound.concentrations[0].error, 0.0)
        # self.assertEqual(compound.concentrations[0].strain, 'K12 NCM3722')
        # self.assertEqual(compound.concentrations[0].growth_status, 'Mid-Log Phase')
        # self.assertTrue(compound.concentrations[0].media.startswith('Gutnick minimal complete medium'))
        # self.assertEqual(compound.concentrations[0].temperature, 37)
        # self.assertEqual(compound.concentrations[0].growth_system, 'Shake flask and filter culture')
        # self.assertEqual(len(compound.concentrations[0].references), 1)
        # self.assertEqual(compound.concentrations[0].references[0].namespace, 'pubmed')
        # self.assertEqual(compound.concentrations[0].references[0].id, '19561621')

        # compound with comment
        # compound = session.query(ecmdb.Compound).filter_by(id='M2MDB000734').first()

    def test_cascade_concentrations(self):
        src = self.src
        session = src.session

        n_compound = session.query(ymdb.Compound).count()
        n_concentration = session.query(ymdb.Concentration).count()

        compound = session.query(ymdb.Compound).filter_by(id='YMDB000010').first()
        n_concentration_compound = len(compound.concentrations)
        session.delete(compound)

        self.assertEqual(session.query(ymdb.Compound).count(), n_compound - 1)
        self.assertEqual(session.query(ymdb.Concentration).count(), n_concentration - n_concentration_compound)

    def test_cascade_compounds_1(self):
        src = self.src
        session = src.session

        n_compound = session.query(ymdb.Compound).count()
        n_concentration = session.query(ymdb.Concentration).count()

        compound = session.query(ymdb.Compound).filter_by(id='YMDB000010').first()
        n_concentration_compound = len(compound.concentrations)
        compound.concentrations[0].compound = None

        self.assertEqual(len(compound.concentrations), n_concentration_compound - 1)
        self.assertEqual(session.query(ymdb.Compound).count(), n_compound)
        self.assertEqual(session.query(ymdb.Concentration).count(), n_concentration - 1)

    def test_cascade_compounds_2(self):
        src = self.src
        session = src.session

        n_compound = session.query(ymdb.Compound).count()
        n_concentration = session.query(ymdb.Concentration).count()

        compound = session.query(ymdb.Compound).filter_by(id='YMDB000010').first()
        n_concentration_compound = len(compound.concentrations)
        concentration = compound.concentrations[0]
        concentration.compound = None
        session.delete(concentration)

        self.assertEqual(len(compound.concentrations), n_concentration_compound - 1)
        self.assertEqual(session.query(ymdb.Compound).count(), n_compound)
        self.assertEqual(session.query(ymdb.Concentration).count(), n_concentration - 1)

    def test_cascade_compounds_3(self):
        src = self.src
        session = src.session

        n_compound = session.query(ymdb.Compound).count()
        n_concentration = session.query(ymdb.Concentration).count()

        concentration = session.query(ymdb.Concentration).first()
        session.delete(concentration)

        self.assertEqual(session.query(ymdb.Compound).count(), n_compound)
        self.assertEqual(session.query(ymdb.Concentration).count(), n_concentration - 1)


# class TestYmdbFromBlank(unittest.TestCase):

#     def setUp(self):
#         self.cache_dirname = tempfile.mkdtemp()

#     def tearDown(self):
#         shutil.rmtree(self.cache_dirname)

#     def test_commit_intermediate_results(self):
#         src = ymdb.Ymdb(cache_dirname=self.cache_dirname, download_backups=False, load_content=False, verbose=True, max_entries=5,
#                           commit_intermediate_results=True)
#         src.load_content()

#     @unittest.skip('Skip because this is a long test')
#     def test_download_all(self):
#         src = ymdb.Ymdb(download_backups=False, load_content=True, clear_content=True, clear_requests_cache=False, verbose=True)
#         src.upload_backups()

#         session = src.session
#         self.assertGreater(session.query(ymdb.Compound).count(), 3500)

#     def test_download_from_backup(self):
#         src = ymdb.Ymdb(cache_dirname=self.cache_dirname)
#         session = src.session

#         self.assertGreater(session.query(ymdb.Compound).count(), 3500)

# class TestEcmdbFromCache(unittest.TestCase):
#     """
#     Quick test to ensure Ecmdb on Karr Lab Server is correct

#     """

#     @classmethod
#     def setUpClass(cls):
#         cls.cache_dirname = tempfile.mkdtemp()
#         cls.src = ymdb.Ymdb(cache_dirname=cls.cache_dirname, download_backups=True, load_content=False, verbose=True)

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
