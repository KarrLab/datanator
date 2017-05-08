""" Tests of ecmdb

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-05
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.data_source import ecmdb
import datetime
import dateutil
import os
import shutil
import tempfile
import unittest


class TestEcmdb(unittest.TestCase):

    def setUp(self):
        self.cache_dir = tempfile.mkdtemp()
        self.db_filename = os.path.join(self.cache_dir, 'db.sqlite')
        self.req_filename = os.path.join(self.cache_dir, 'req.sqlite')

    def tearDown(self):
        shutil.rmtree(self.cache_dir)

    def test_download_from_remote(self):
        engine = ecmdb.get_engine(filename=self.db_filename)
        session = ecmdb.get_session(engine=engine, auto_download=False)
        downloader = ecmdb.Downloader(session, requests_cache_filename=self.req_filename, max_entries=12, verbose=True)
        downloader.download()
        session.close()
        engine.dispose()

        engine = ecmdb.get_engine(filename=self.db_filename)
        session = ecmdb.get_session(engine=engine, auto_download=False)
        q = session.query(ecmdb.Compound)

        # check compound ids
        self.assertEqual(q.count(), downloader.max_entries)
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
        self.assertEqual(compound._structure_formula_connectivity, 'C4O3/c1-2-3(5)4(6)7')
        self.assertEqual([c.name for c in compound.compartments], ['Cytosol'])

        self.assertEqual(compound.concentrations, [])

        self.assertEqual(set([(xr.namespace, xr.id) for xr in compound.cross_references]), set([
            ('biocyc', '2-OXOBUTANOATE'),
            ('cas', '600-18-0'),
            ('chebi', 'CHEBI:16763'),
            ('chemspider', '57'),
            ('ligandexpo', '2KT'),
            ('hmdb', 'HMDB00005'),
            ('kegg.compound', 'C00109'),
            ('pubchem.compound', '58'),
            ('wikipedia.en', 'Alpha-ketobutyric_acid'),
        ]))
        self.assertEqual(compound.comment, None)

        self.assertEqual(compound.created, dateutil.parser.parse('2012-05-31 09:55:11 -0600').replace(tzinfo=None))
        self.assertEqual(compound.updated, dateutil.parser.parse('2015-06-03 15:00:41 -0600').replace(tzinfo=None))
        self.assertLess((datetime.datetime.utcnow() - compound.downloaded).total_seconds(), 60)

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

    @unittest.skip('Skip because this is a long test')
    def test_download_all(self):
        engine = ecmdb.get_engine()
        session = ecmdb.get_session(engine=engine)
        downloader = ecmdb.Downloader(session, verbose=True)
        downloader.download(update_database=True, update_requests=False)
        session.close()
        engine.dispose()

        session = ecmdb.get_session()
        self.assertGreater(session.query(ecmdb.Compound).count(), 3500)

    def test_download_from_cache(self):
        engine = ecmdb.get_engine(filename=self.db_filename)
        session = ecmdb.get_session(engine=engine)

        self.assertGreater(session.query(ecmdb.Compound).count(), 3500)
