# -*- coding: utf-8 -*-

""" Tests of sabio_rk

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-26
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.data_source import sabio_rk
from kinetic_datanator.data_source.sabio_rk import (Entry, Compartment, Compound, Enzyme, Reaction,
                                                    ReactionParticipant, KineticLaw, Parameter, Resource, Downloader)
import datetime
import math
import os
import six
import sqlalchemy
import tempfile
import unittest


if six.PY3:
    from test.support import EnvironmentVarGuard
else:
    from test.test_support import EnvironmentVarGuard


class TestDownloader(unittest.TestCase):

    def setUp(self):
        _, self.engine_filename = tempfile.mkstemp()
        os.remove(self.engine_filename)

        _, self.requests_cache_name = tempfile.mkstemp()
        os.remove(self.requests_cache_name)

    def tearDown(self):
        if os.path.isfile(self.engine_filename):
            os.remove(self.engine_filename)

        if os.path.isfile(self.requests_cache_name + '.sqlite'):
            os.remove(self.requests_cache_name + '.sqlite')

    def test_download_kinetic_law_ids(self):
        engine = sabio_rk.get_engine(filename=self.engine_filename)
        session = sabio_rk.get_session(engine=engine, auto_download=False)
        downloader = sabio_rk.Downloader(session, requests_cache_name=self.requests_cache_name, max_laws=10,
                                         index_batch_size=5, webservice_batch_size=5, excel_batch_size=5, compound_batch_size=5)

        ids = downloader.download_kinetic_law_ids()
        self.assertEqual(ids, [1, 10, 100, 1000, 10000, 10001, 10002, 10003, 10004, 10005])

    def test_download_kinetic_laws_and_compounds(self):
        engine = sabio_rk.get_engine(filename=self.engine_filename)
        session = sabio_rk.get_session(engine=engine, auto_download=False)
        downloader = sabio_rk.Downloader(session, requests_cache_name=self.requests_cache_name)
        downloader.download_kinetic_laws([1])

        """ compartments """
        q = session.query(Compartment).filter_by(name='Cell')
        self.assertEqual(q.count(), 0)

        """ compounds """
        c = session.query(Compound).filter_by(id=40).first()
        self.assertEqual(c.name, 'H2O')
        xrs = [dict(namespace=xr.namespace, id=xr.id) for xr in c.cross_references]
        self.assertEqual(sorted(xrs, key=lambda xr: xr['namespace']), [
            dict(namespace='chebi', id='CHEBI:15377'),
            dict(namespace='kegg.compound', id='C00001'),
        ])

        self.assertIsInstance(c.created, datetime.datetime)
        self.assertIsInstance(c.modified, datetime.datetime)
        self.assertLess((datetime.datetime.utcnow() - c.created).total_seconds(), 120)
        self.assertLess((datetime.datetime.utcnow() - c.modified).total_seconds(), 120)
        h20_created = c.created

        c = session.query(Compound).filter_by(id=2562).first()
        self.assertEqual(c.name, 'm n n+m Peptide')
        xrs = [dict(namespace=xr.namespace, id=xr.id) for xr in c.cross_references]
        self.assertEqual(sorted(xrs, key=lambda xr: xr['namespace']), [
            dict(namespace='chebi', id='CHEBI:16670'),
            dict(namespace='kegg.compound', id='C00012'),
        ])

        c = session.query(Compound).filter_by(id=20035).first()
        self.assertEqual(c.name, 'N-(Carbobenzoxy)-Leu-Leu-Phe-trifluoromethylketone')
        self.assertEqual(c.cross_references, [])

        """ enzymes """
        e = session.query(Enzyme).filter_by(id=147631).first()
        self.assertEqual(e.name, 'subtilisin')
        self.assertEqual(e.cross_references, [])

        """ reactions """
        rxn = session.query(Reaction).filter_by(id=6570).first()
        cpd_40 = session.query(Compound).filter_by(id=40).first()
        cpd_2562 = session.query(Compound).filter_by(id=2562).first()
        cpd_20035 = session.query(Compound).filter_by(id=20035).first()
        enz_147631 = session.query(Enzyme).filter_by(id=147631).first()

        self.assertEqual(len(rxn.reactants), 2)
        self.assertEqual(rxn.reactants[0].compound, cpd_2562)
        self.assertEqual(rxn.reactants[0].compartment, None)
        self.assertTrue(rxn.reactants[0].coefficient is None or math.isnan(rxn.reactants[0].coefficient))
        self.assertEqual(rxn.reactants[1].compound, cpd_40)
        self.assertEqual(rxn.reactants[1].coefficient, 1.)
        self.assertEqual(rxn.reactants[1].compartment, None)

        self.assertEqual(len(rxn.products), 2)
        self.assertEqual(rxn.products[0].compound, cpd_2562)
        self.assertEqual(rxn.products[1].compound, cpd_2562)
        self.assertEqual(rxn.reactants[0].compartment, None)
        self.assertEqual(rxn.reactants[1].compartment, None)
        self.assertTrue(rxn.products[0].coefficient is None or math.isnan(rxn.products[0].coefficient))
        self.assertTrue(rxn.products[1].coefficient is None or math.isnan(rxn.products[1].coefficient))

        self.assertEqual([(r.compound, r.coefficient) for r in rxn.kinetic_law.modifiers], [
            (cpd_20035, None),
        ])
        self.assertEqual([dict(namespace=xr.namespace, id=xr.id) for xr in rxn.cross_references], [
            dict(namespace='ec-code', id='3.4.21.62'),
        ])

        """ kinetic laws """
        l = session.query(KineticLaw).filter_by(id=1).first()
        self.assertEqual(l.reaction.id, 6570)
        self.assertEqual(l.enzyme, enz_147631)
        self.assertEqual(l.enzyme_compartment, None)
        self.assertEqual(l.equation, None)

        self.assertEqual(len(l.parameters), 4)

        p = session.query(Parameter).filter_by(kinetic_law=l, type='kcat/Km').first()
        self.assertEqual(p.compound, cpd_2562)
        self.assertEqual(p.compartment, None)
        self.assertEqual(p.value, 120000.)
        self.assertEqual(p.units, 'M^(-1)*s^(-1)')

        p = session.query(Parameter).filter_by(kinetic_law=l, type='kcat').first()
        self.assertEqual(p.compound, None)
        self.assertEqual(p.compartment, None)
        self.assertEqual(p.value, 220.)
        self.assertEqual(p.units, 's^(-1)')

        p = session.query(Parameter).filter_by(kinetic_law=l, type='Ki').first()
        self.assertEqual(p.compound, cpd_20035)
        self.assertEqual(p.compartment, None)
        self.assertEqual(p.value, 9E-9)
        self.assertEqual(p.units, 'M')

        p = session.query(Parameter).filter_by(kinetic_law=l, type='Km').first()
        self.assertEqual(p.compound, cpd_2562)
        self.assertEqual(p.compartment, None)
        self.assertEqual(p.value, 0.0019)
        self.assertEqual(p.units, 'M')

        self.assertEqual(l.taxon, 1467)
        self.assertEqual(l.taxon_wildtype, True)
        self.assertEqual(l.taxon_variant, 'variant DSAI (N76D/N87S/S103A/V104I)')
        self.assertEqual(l.temperature, 25.0)
        self.assertEqual(l.ph, 7.5)
        self.assertEqual(l.media, '50 mM potassium phosphate, 4 % DMSO')
        self.assertEqual([dict(n=ref.namespace, i=ref.id) for ref in l.references], [
            dict(n='pubmed', i='12962477')
        ])

        # download compounds
        downloader.download_compounds()
        c = session.query(Compound).filter_by(id=40).first()
        self.assertEqual(c.name, 'H2O')
        self.assertEqual(c._is_name_ambiguous, False)
        self.assertEqual([s.name for s in c.synonyms], ['Water'])
        self.assertEqual(c.get_inchi_structures(), ['InChI=1S/H2O/h1H2'])
        self.assertEqual(c.get_smiles_structures(), ['[H]O[H]'])
        xrs = [dict(namespace=xr.namespace, id=xr.id) for xr in c.cross_references]
        self.assertEqual(sorted(xrs, key=lambda xr: xr['namespace']), [
            dict(namespace='chebi', id='CHEBI:15377'),
            dict(namespace='kegg.compound', id='C00001'),
            dict(namespace='pubchem.substance', id='3303'),
        ])

        self.assertEqual(c.created, h20_created)
        self.assertIsInstance(c.modified, datetime.datetime)
        self.assertLess((datetime.datetime.utcnow() - c.modified).total_seconds(), 120)

    def test_download_kinetic_laws_multiple(self):
        engine = sabio_rk.get_engine(filename=self.engine_filename)
        session = sabio_rk.get_session(engine=engine, auto_download=False)
        downloader = sabio_rk.Downloader(session, requests_cache_name=self.requests_cache_name)
        downloader.download_kinetic_laws([2, 10026])

        """ compartments """
        self.assertEqual(session.query(Compartment).filter_by(name='Cell').count(), 0)

        cytosol = session.query(Compartment).filter_by(name='cytosol').first()

        """ reactions """
        l = session.query(KineticLaw).filter_by(id=10026).first()
        for part in l.reaction.reactants:
            self.assertEqual(part.compartment, cytosol)
        for part in l.reaction.products:
            self.assertEqual(part.compartment, cytosol)

        """ kinetic laws """
        l = session.query(KineticLaw).filter_by(id=2).first()
        self.assertEqual(l.enzyme_compartment, None)
        self.assertEqual(l.tissue, None)
        self.assertEqual(l.mechanism, None)
        self.assertFalse(l.taxon_wildtype)
        self.assertEqual(l.taxon_variant, 'S156E of subtilisin DSAI (N76D/N87S/S103A/V104I)')
        self.assertEqual(l.equation, None)
        self.assertEqual(l.media, '50 mM potassium phosphate, 4 % DMSO')

        l = session.query(KineticLaw).filter_by(id=10026).first()
        self.assertEqual(l.enzyme_compartment, cytosol)
        self.assertEqual(l.tissue, 'epidermal cell')
        self.assertEqual(l.mechanism, 'Michaelis-Menten')
        self.assertTrue(l.taxon_wildtype)
        self.assertEqual(l.taxon_variant, '')
        self.assertEqual(l.equation, 'Vmax * A / (Km + A)')
        self.assertEqual(l.media, u'1.25 mM CaCl2, 1 mM dithiothretiol, 10 µM tetrahydrobiopterin, 10 units/ml calmodulin')

        for param in l.parameters:
            if param.compound:
                self.assertEqual(param.compartment, cytosol)

    def test_download_kinetic_laws_modifier_type(self):
        engine = sabio_rk.get_engine(filename=self.engine_filename)
        session = sabio_rk.get_session(engine=engine, auto_download=False)
        downloader = sabio_rk.Downloader(session, requests_cache_name=self.requests_cache_name)
        downloader.download_kinetic_laws([10054])

        l = session.query(KineticLaw).filter_by(id=10054).first()
        self.assertEqual(l.enzyme_type, 'Modifier-Catalyst')

    def test_infer_compound_structures_from_names(self):
        engine = sabio_rk.get_engine(filename=self.engine_filename)
        session = sabio_rk.get_session(engine=engine, auto_download=False)
        downloader = sabio_rk.Downloader(session, requests_cache_name=self.requests_cache_name)

        compound_unknown = sabio_rk.Compound(name='Unknown')
        compound_no_pubchem = sabio_rk.Compound(name='a random string: adfkja;uvhayocbvadf')
        compound_one_pubchem = sabio_rk.Compound(name='water')
        compound_multiple_pubchem = sabio_rk.Compound(name='Phenyllactate')
        session.add(compound_unknown)
        session.add(compound_no_pubchem)
        session.add(compound_one_pubchem)
        session.add(compound_multiple_pubchem)

        compounds = [compound_unknown, compound_no_pubchem, compound_one_pubchem, compound_multiple_pubchem]
        downloader.infer_compound_structures_from_names(compounds)

        self.assertEqual(compound_unknown.cross_references, [])
        self.assertEqual(compound_unknown.structures, [])

        self.assertEqual(compound_no_pubchem.cross_references, [])
        self.assertEqual(compound_no_pubchem.structures, [])

        self.assertEqual(set([(xr.namespace, xr.id) for xr in compound_one_pubchem.cross_references]), set([
            ('pubchem.compound', '962'),
        ]))
        self.assertEqual(set([(s.format, s.value) for s in compound_one_pubchem.structures]), set([
            ('inchi', 'InChI=1S/H2O/h1H2'),
        ]))

        self.assertEqual(set([(xr.namespace, xr.id) for xr in compound_multiple_pubchem.cross_references]), set([
            ('pubchem.compound', '3848'),
            ('pubchem.compound', '643327'),
            ('pubchem.compound', '4060207'),
        ]))
        self.assertEqual(set([(s.format, s.value) for s in compound_multiple_pubchem.structures]), set([
            ('inchi', 'InChI=1S/C9H10O3/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8,10H,6H2,(H,11,12)'),
            ('inchi', 'InChI=1S/C9H10O3/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8,10H,6H2,(H,11,12)/t8-/m1/s1'),
            ('inchi', 'InChI=1S/C9H10O3/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8,10H,6H2,(H,11,12)/p-1'),
        ]))

    def test_download(self):
        # get some kinetic laws
        engine = sabio_rk.get_engine(filename=self.engine_filename)
        session = sabio_rk.get_session(engine=engine, auto_download=False)
        downloader = sabio_rk.Downloader(session, requests_cache_name=self.requests_cache_name,
                                         max_laws=9,
                                         index_batch_size=3, webservice_batch_size=3,
                                         excel_batch_size=3, compound_batch_size=3,
                                         verbose=True)
        downloader.download()
        session.close()
        engine.dispose()

        self.assertTrue(os.path.isfile(self.engine_filename))
        engine2 = sabio_rk.get_engine(filename=self.engine_filename)
        session2 = sabio_rk.get_session(engine=engine2, auto_download=False)
        self.assertEqual(session2.query(KineticLaw).count(), downloader.max_laws)

        # get some more
        engine = sabio_rk.get_engine(filename=self.engine_filename)
        session = sabio_rk.get_session(engine=engine, auto_download=False)
        downloader = sabio_rk.Downloader(session, requests_cache_name=self.requests_cache_name,
                                         max_laws=18,
                                         index_batch_size=3, webservice_batch_size=3,
                                         excel_batch_size=3, compound_batch_size=3,
                                         verbose=True)
        downloader.download()
        session.close()
        engine.dispose()

        self.assertTrue(os.path.isfile(self.engine_filename))
        engine2 = sabio_rk.get_engine(filename=self.engine_filename)
        session2 = sabio_rk.get_session(engine=engine2, auto_download=False)
        self.assertEqual(session2.query(KineticLaw).count(), downloader.max_laws)


class TestBackupAndInstall(unittest.TestCase):

    def setUp(self):
        _, self.engine_filename = tempfile.mkstemp()
        os.remove(self.engine_filename)

        _, self.engine_filename_2 = tempfile.mkstemp()
        os.remove(self.engine_filename_2)

        _, self.engine_filename_3 = tempfile.mkstemp()
        os.remove(self.engine_filename_3)

        _, self.engine_filename_4 = tempfile.mkstemp()
        os.remove(self.engine_filename_4)

        _, self.requests_cache_name = tempfile.mkstemp()
        os.remove(self.requests_cache_name)

        # download some kinetic laws
        engine = sabio_rk.get_engine(filename=self.engine_filename)
        session = sabio_rk.get_session(engine=engine, auto_download=False)
        downloader = sabio_rk.Downloader(session, requests_cache_name=self.requests_cache_name,
                                         max_laws=2,
                                         index_batch_size=10, webservice_batch_size=10,
                                         excel_batch_size=10, compound_batch_size=10,
                                         verbose=True)
        downloader.download()
        session.close()
        engine.dispose()

    def tearDown(self):
        if os.path.isfile(self.engine_filename):
            os.remove(self.engine_filename)

        if os.path.isfile(self.engine_filename_2):
            os.remove(self.engine_filename_2)

        if os.path.isfile(self.engine_filename_3):
            os.remove(self.engine_filename_3)

        if os.path.isfile(self.requests_cache_name + '.sqlite'):
            os.remove(self.requests_cache_name + '.sqlite')

    def test(self):
        env = EnvironmentVarGuard()

        if not os.getenv('CODE_SERVER_TOKEN'):
            with open('tests/fixtures/secret/CODE_SERVER_TOKEN', 'r') as file:
                env.set('CODE_SERVER_TOKEN', file.read().rstrip())

        # backup and download
        sabio_rk.backup(filename=self.engine_filename, arcname='test.sabio_rk.sqlite')
        sabio_rk.download(filename=self.engine_filename_2, arcname='test.sabio_rk.sqlite')
        self.assertTrue(os.path.isfile(self.engine_filename_2))

        engine = sabio_rk.get_engine(filename=self.engine_filename_2)
        session = sabio_rk.get_session(engine=engine)
        self.assertEqual(session.query(KineticLaw).count(), 2)

        # setup with download and update
        engine = sabio_rk.get_engine(filename=self.engine_filename_3)
        session = sabio_rk.get_session(engine=engine, auto_download=True, force_update=True, arcname='test.sabio_rk.sqlite', max_laws=4)
        self.assertEqual(session.query(KineticLaw).count(), 4)

        # setup with update
        engine = sabio_rk.get_engine(filename=self.engine_filename_4)
        session = sabio_rk.get_session(engine=engine, auto_download=False, auto_update=True, max_laws=4)
        self.assertEqual(session.query(KineticLaw).count(), 4)


class TestAll(unittest.TestCase):

    @unittest.skip('Skip this test because it is long')
    def test_update_all_and_backup(self):
        session = sabio_rk.get_session(auto_download=False)
        downloader = sabio_rk.Downloader(session, verbose=True)

        downloader.download(update_database=True, update_requests=True)
        self.assertGreaterEqual(session.query(KineticLaw).count(), 55000)

        env = EnvironmentVarGuard()
        if not os.getenv('CODE_SERVER_TOKEN'):
            with open('tests/fixtures/secret/CODE_SERVER_TOKEN', 'r') as file:
                env.set('CODE_SERVER_TOKEN', file.read().rstrip())

        sabio_rk.backup()

    def test_download_full_database(self):
        env = EnvironmentVarGuard()
        if not os.getenv('CODE_SERVER_TOKEN'):
            with open('tests/fixtures/secret/CODE_SERVER_TOKEN', 'r') as file:
                env.set('CODE_SERVER_TOKEN', file.read().rstrip())

        _, filename = tempfile.mkstemp()
        os.remove(filename)

        engine = sabio_rk.get_engine(filename=filename)
        session = sabio_rk.get_session(engine=engine)
        self.assertGreaterEqual(session.query(KineticLaw).count(), 55000)

        session.close()
        engine.dispose()
        os.remove(filename)
