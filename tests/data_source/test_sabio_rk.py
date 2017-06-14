# -*- coding: utf-8 -*-

""" Tests of sabio_rk

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-26
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.data_source import sabio_rk
from kinetic_datanator.data_source.sabio_rk import (Entry, Compartment, Compound, Enzyme,
                                                    ReactionParticipant, KineticLaw, Parameter, Resource)
from kinetic_datanator.util import warning_util
import datetime
import math
import os
import scipy.constants
import shutil
import sqlalchemy
import tempfile
import unittest

warning_util.disable_warnings()


class TestDownloader(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_load_kinetic_law_ids(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backup=False, load_content=False, max_entries=10,
                               webservice_batch_size=1, excel_batch_size=5)
        ids = src.load_kinetic_law_ids()
        self.assertEqual(ids[0:10], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        self.assertGreater(len(ids), 55000)

    def test_load_kinetic_laws_and_compounds(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backup=False, load_content=False)
        src.load_kinetic_laws([1])

        session = src.session

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
        self.assertLess((datetime.datetime.utcnow() - c.created).total_seconds(), 240)
        self.assertLess((datetime.datetime.utcnow() - c.modified).total_seconds(), 240)
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
        law = session.query(KineticLaw).filter_by(id=1).first()
        cpd_40 = session.query(Compound).filter_by(id=40).first()
        cpd_2562 = session.query(Compound).filter_by(id=2562).first()
        cpd_20035 = session.query(Compound).filter_by(id=20035).first()
        enz_147631 = session.query(Enzyme).filter_by(id=147631).first()

        self.assertEqual(next(xr.id for xr in law.cross_references if xr.namespace == 'sabiork.reaction'), '6570')

        self.assertEqual(len(law.reactants), 2)
        self.assertEqual(law.reactants[0].compound, cpd_2562)
        self.assertEqual(law.reactants[0].compartment, None)
        self.assertTrue(law.reactants[0].coefficient is None or math.isnan(law.reactants[0].coefficient))
        self.assertEqual(law.reactants[1].compound, cpd_40)
        self.assertEqual(law.reactants[1].coefficient, 1.)
        self.assertEqual(law.reactants[1].compartment, None)

        self.assertEqual(len(law.products), 2)
        self.assertEqual(law.products[0].compound, cpd_2562)
        self.assertEqual(law.products[1].compound, cpd_2562)
        self.assertEqual(law.products[0].compartment, None)
        self.assertEqual(law.products[1].compartment, None)
        self.assertTrue(law.products[0].coefficient is None or math.isnan(law.products[0].coefficient))
        self.assertTrue(law.products[1].coefficient is None or math.isnan(law.products[1].coefficient))

        self.assertEqual([(r.compound, r.coefficient) for r in law.modifiers], [
            (cpd_20035, None),
        ])

        for mod in law.modifiers:
            self.assertEqual(mod.modifier_kinetic_law, law)

        """ kinetic laws """
        l = session.query(KineticLaw).filter_by(id=1).first()
        self.assertEqual(l.enzyme, enz_147631)
        self.assertEqual(l.enzyme_compartment, None)
        self.assertEqual(l.equation, None)
        self.assertEqual([dict(namespace=xr.namespace, id=xr.id) for xr in l.cross_references], [
            dict(namespace='ec-code', id='3.4.21.62'),
            dict(namespace='sabiork.reaction', id='6570'),
        ])

        self.assertEqual(len(l.parameters), 4)

        p = session.query(Parameter).filter_by(kinetic_law=l, observed_name='kcat/Km').first()
        self.assertEqual(p.compound, cpd_2562)
        self.assertEqual(p.compartment, None)
        self.assertEqual(p.observed_type, 302)
        self.assertEqual(p.observed_value, 120000.)
        self.assertEqual(p.observed_error, None)
        self.assertEqual(p.observed_units, 'M^(-1)*s^(-1)')
        self.assertEqual(p.name, None)
        self.assertEqual(p.type, None)
        self.assertEqual(p.value, None)
        self.assertEqual(p.error, None)
        self.assertEqual(p.units, None)

        p = session.query(Parameter).filter_by(kinetic_law=l, observed_name='kcat').first()
        self.assertEqual(p.compound, None)
        self.assertEqual(p.compartment, None)
        self.assertEqual(p.observed_type, 25)
        self.assertEqual(p.observed_value, 220.)
        self.assertEqual(p.observed_error, None)
        self.assertEqual(p.observed_units, 's^(-1)')

        p = session.query(Parameter).filter_by(kinetic_law=l, observed_name='Ki').first()
        self.assertEqual(p.compound, cpd_20035)
        self.assertEqual(p.compartment, None)
        self.assertEqual(p.observed_type, 261)
        self.assertEqual(p.observed_value, 9E-9)
        self.assertEqual(p.observed_error, None)
        self.assertEqual(p.observed_units, 'M')

        p = session.query(Parameter).filter_by(kinetic_law=l, observed_name='Km').first()
        self.assertEqual(p.compound, cpd_2562)
        self.assertEqual(p.compartment, None)
        self.assertEqual(p.observed_type, 27)
        self.assertEqual(p.observed_value, 0.0019)
        self.assertEqual(p.observed_error, None)
        self.assertEqual(p.observed_units, 'M')

        self.assertEqual(l.taxon, 1467)
        self.assertEqual(l.taxon_wildtype, True)
        self.assertEqual(l.taxon_variant, 'variant DSAI (N76D/N87S/S103A/V104I)')
        self.assertEqual(l.temperature, 25.0)
        self.assertEqual(l.ph, 7.5)
        self.assertEqual(l.media, '50 mM potassium phosphate, 4 % DMSO')
        self.assertEqual([dict(n=ref.namespace, i=ref.id) for ref in l.references], [
            dict(n='pubmed', i='12962477')
        ])
        for ref in l.references:
            self.assertIn(l, ref.kinetic_laws)

        # download compounds
        src.load_compounds()
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
        self.assertLess((datetime.datetime.utcnow() - c.modified).total_seconds(), 240)

    def test_load_and_update_kinetic_laws_and_compounds(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backup=False, load_content=False)
        session = src.session

        src.load_kinetic_laws([1])
        src.load_compounds()
        src.load_missing_kinetic_law_information_from_tsv([1])
        src.normalize_kinetic_laws([1])

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
            dict(namespace='pubchem.substance', id='3303'),
        ])

        self.assertIsInstance(c.created, datetime.datetime)
        self.assertIsInstance(c.modified, datetime.datetime)
        self.assertLess((datetime.datetime.utcnow() - c.created).total_seconds(), 240)
        self.assertLess((datetime.datetime.utcnow() - c.modified).total_seconds(), 240)
        h20_created = c.created

        c = session.query(Compound).filter_by(id=2562).first()
        self.assertEqual(c.name, 'Peptide')
        xrs = [dict(namespace=xr.namespace, id=xr.id) for xr in c.cross_references]
        self.assertEqual(sorted(xrs, key=lambda xr: xr['namespace']), [
            dict(namespace='chebi', id='CHEBI:16670'),
            dict(namespace='kegg.compound', id='C00012'),
            dict(namespace='pubchem.substance', id='3314'),
        ])

        c = session.query(Compound).filter_by(id=20035).first()
        self.assertEqual(c.name, 'N-(Carbobenzoxy)-Leu-Leu-Phe-trifluoromethylketone')
        self.assertEqual(c.cross_references, [])

        """ enzymes """
        e = session.query(Enzyme).filter_by(id=147631).first()
        self.assertEqual(e.name, 'subtilisin')
        self.assertEqual(e.cross_references, [])

        """ reactions """
        law = session.query(KineticLaw).filter_by(id=1).first()
        cpd_40 = session.query(Compound).filter_by(id=40).first()
        cpd_2562 = session.query(Compound).filter_by(id=2562).first()
        cpd_20035 = session.query(Compound).filter_by(id=20035).first()
        enz_147631 = session.query(Enzyme).filter_by(id=147631).first()

        self.assertEqual(next(xr.id for xr in law.cross_references if xr.namespace == 'sabiork.reaction'), '6570')

        self.assertEqual(len(law.reactants), 2)
        self.assertEqual(law.reactants[0].compound, cpd_2562)
        self.assertEqual(law.reactants[0].compartment, None)
        self.assertTrue(law.reactants[0].coefficient is None or math.isnan(law.reactants[0].coefficient))
        self.assertEqual(law.reactants[1].compound, cpd_40)
        self.assertEqual(law.reactants[1].coefficient, 1.)
        self.assertEqual(law.reactants[1].compartment, None)

        self.assertEqual(len(law.products), 2)
        self.assertEqual(law.products[0].compound, cpd_2562)
        self.assertEqual(law.products[1].compound, cpd_2562)
        self.assertEqual(law.products[0].compartment, None)
        self.assertEqual(law.products[1].compartment, None)
        self.assertTrue(law.products[0].coefficient is None or math.isnan(law.products[0].coefficient))
        self.assertTrue(law.products[1].coefficient is None or math.isnan(law.products[1].coefficient))

        self.assertEqual([(r.compound, r.coefficient) for r in law.modifiers], [
            (cpd_20035, None),
        ])

        for mod in law.modifiers:
            self.assertEqual(mod.modifier_kinetic_law, law)

        """ kinetic laws """
        l = session.query(KineticLaw).filter_by(id=1).first()
        self.assertEqual(l.enzyme, enz_147631)
        self.assertEqual(l.enzyme_compartment, None)
        self.assertEqual(l.equation, None)
        self.assertEqual([dict(namespace=xr.namespace, id=xr.id) for xr in l.cross_references], [
            dict(namespace='ec-code', id='3.4.21.62'),
            dict(namespace='sabiork.reaction', id='6570'),
        ])

        self.assertEqual(len(l.parameters), 4)

        p = session.query(Parameter).filter_by(kinetic_law=l, observed_name='kcat/Km').first()
        self.assertEqual(p.compound, cpd_2562)
        self.assertEqual(p.compartment, None)
        self.assertEqual(p.observed_type, 302)
        self.assertEqual(p.observed_value, 120000.)
        self.assertEqual(p.observed_error, None)
        self.assertEqual(p.observed_units, 'M^(-1)*s^(-1)')
        self.assertEqual(p.name, None)
        self.assertEqual(p.type, None)
        self.assertEqual(p.value, None)
        self.assertEqual(p.error, None)
        self.assertEqual(p.units, None)

        p = session.query(Parameter).filter_by(kinetic_law=l, observed_name='kcat').first()
        self.assertEqual(p.compound, None)
        self.assertEqual(p.compartment, None)
        self.assertEqual(p.observed_type, 25)
        self.assertEqual(p.observed_value, 220.)
        self.assertEqual(p.observed_error, None)
        self.assertEqual(p.observed_units, 's^(-1)')
        self.assertEqual(p.name, 'k_cat')
        self.assertEqual(p.type, 25)
        self.assertEqual(p.value, 220.)
        self.assertEqual(p.error, None)
        self.assertEqual(p.units, 's^(-1)')

        p = session.query(Parameter).filter_by(kinetic_law=l, observed_name='Ki').first()
        self.assertEqual(p.compound, cpd_20035)
        self.assertEqual(p.compartment, None)
        self.assertEqual(p.observed_type, 261)
        self.assertEqual(p.observed_value, 1.8E-8)
        self.assertEqual(p.observed_error, 0.0)
        self.assertEqual(p.observed_units, 'M')
        self.assertEqual(p.name, 'k_i')
        self.assertEqual(p.type, 261)
        self.assertEqual(p.value, 1.8E-8)
        self.assertEqual(p.error, 0.0)
        self.assertEqual(p.units, 'M')

        p = session.query(Parameter).filter_by(kinetic_law=l, observed_name='Km').first()
        self.assertEqual(p.compound, cpd_2562)
        self.assertEqual(p.compartment, None)
        self.assertEqual(p.observed_type, 27)
        self.assertEqual(p.observed_value, 0.0019)
        self.assertEqual(p.observed_error, None)
        self.assertEqual(p.observed_units, 'M')
        self.assertEqual(p.name, 'k_m')
        self.assertEqual(p.type, 27)
        self.assertEqual(p.value, 0.0019)
        self.assertEqual(p.error, None)
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
        for ref in l.references:
            self.assertIn(l, ref.kinetic_laws)

        # download compounds
        src.load_compounds()
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
        self.assertLess((datetime.datetime.utcnow() - c.modified).total_seconds(), 240)

    def test_load_kinetic_laws_multiple(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backup=False, load_content=False)
        session = src.session

        ids = [2, 10026]
        src.load_kinetic_laws(ids)
        src.load_compounds()
        src.load_missing_kinetic_law_information_from_tsv(ids)
        src.normalize_kinetic_laws(ids)

        """ compartments """
        self.assertEqual(session.query(Compartment).filter_by(name='Cell').count(), 0)

        cytosol = session.query(Compartment).filter_by(name='cytosol').first()

        """ reactions """
        l = session.query(KineticLaw).filter_by(id=10026).first()
        for part in l.reactants:
            self.assertEqual(part.compartment, cytosol)
        for part in l.products:
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
        self.assertTrue(l.media.startswith('1.25 mM CaCl2'))
        self.assertTrue(l.media.endswith('10 units/ml calmodulin'))

        for param in l.parameters:
            if param.compound:
                self.assertEqual(param.compartment, cytosol)

    def test_load_kinetic_laws_modifier_type(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backup=False, load_content=False)
        src.load_kinetic_laws([10054])

        session = src.session

        l = session.query(KineticLaw).filter_by(id=10054).first()
        self.assertEqual(l.enzyme_type, 'Modifier-Catalyst')

    def test_load_kinetic_laws_with_same_reaction(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backup=False, load_content=False)
        src.load_kinetic_laws([16011, 16013, 16016])

        session = src.session

        l = session.query(KineticLaw).filter_by(id=16011).first()
        self.assertEqual(next(xr.id for xr in l.cross_references if xr.namespace == 'sabiork.reaction'), '9886')

        l = session.query(KineticLaw).filter_by(id=16013).first()
        self.assertEqual(next(xr.id for xr in l.cross_references if xr.namespace == 'sabiork.reaction'), '9886')

        l = session.query(KineticLaw).filter_by(id=16016).first()
        self.assertEqual(next(xr.id for xr in l.cross_references if xr.namespace == 'sabiork.reaction'), '9930')

        r = session.query(Resource).filter_by(namespace='sabiork.reaction', id='9886').first()
        self.assertEqual(set([l.id for l in r.entries]), set([16011, 16013]))

        r = session.query(Resource).filter_by(namespace='sabiork.reaction', id='9930').first()
        self.assertEqual([l.id for l in r.entries], [16016])

    def test_load_kinetic_laws_with_opposite_directions_of_same_reaction(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backup=False, load_content=False)
        src.load_kinetic_laws([46425, 46427])

        session = src.session

        Lactaldehyde = session.query(Compound).filter_by(id=2845).first()

        l = session.query(KineticLaw).filter_by(id=46425).first()
        self.assertEqual(next(xr.id for xr in l.cross_references if xr.namespace == 'sabiork.reaction'), '554')
        self.assertIn(Lactaldehyde, [p.compound for p in l.reactants])
        self.assertNotIn(Lactaldehyde, [p.compound for p in l.products])

        l = session.query(KineticLaw).filter_by(id=46427).first()
        self.assertEqual(next(xr.id for xr in l.cross_references if xr.namespace == 'sabiork.reaction'), '554')
        self.assertNotIn(Lactaldehyde, [p.compound for p in l.reactants])
        self.assertIn(Lactaldehyde, [p.compound for p in l.products])

    def test_load_kinetic_laws_check_units(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backup=False, load_content=False)
        session = src.session

        src.load_kinetic_laws([2023])
        src.load_compounds()
        src.load_missing_kinetic_law_information_from_tsv([2023])
        src.normalize_kinetic_laws([2023])

        law = session.query(KineticLaw).filter_by(id=2023).first()

        self.assertEqual(len(law.parameters), 4)

        param = next(param for param in law.parameters if param.observed_type == 25)
        self.assertEqual(param.observed_name, 'kcat')
        self.assertEqual(param.observed_value, 0.329)
        self.assertEqual(param.observed_error, None)
        self.assertEqual(param.observed_units, 's^(-1)')
        self.assertEqual(param.name, 'k_cat')
        self.assertEqual(param.type, 25)
        self.assertEqual(param.value, 0.329)
        self.assertEqual(param.error, None)
        self.assertEqual(param.units, 's^(-1)')

        param = next(param for param in law.parameters if param.observed_type == 27)
        self.assertEqual(param.observed_name, 'Km')
        self.assertEqual(param.observed_value, 0.232)
        self.assertEqual(param.observed_error, None)
        self.assertEqual(param.observed_units, 'mg/ml')
        self.assertEqual(param.name, None)
        self.assertEqual(param.type, None)
        self.assertEqual(param.value, None)
        self.assertEqual(param.error, None)
        self.assertEqual(param.units, None)

        param = next(param for param in law.parameters if param.observed_type == 261)
        self.assertEqual(param.observed_name, 'Ki')
        self.assertEqual(param.observed_value, 2.4E-9)
        self.assertEqual(param.observed_error, None)
        self.assertEqual(param.observed_units, 'M')
        self.assertEqual(param.name, 'k_i')
        self.assertEqual(param.type, 261)
        self.assertEqual(param.value, 2.4E-9)
        self.assertEqual(param.error, None)
        self.assertEqual(param.units, 'M')

        param = next(param for param in law.parameters if param.observed_type == 302)
        self.assertEqual(param.observed_name, 'kcat/Km')
        self.assertEqual(param.observed_value, 1.42)
        self.assertEqual(param.observed_error, None)
        self.assertEqual(param.observed_units, 'l*g^(-1)*s^(-1)')
        self.assertEqual(param.name, None)
        self.assertEqual(param.type, None)
        self.assertEqual(param.value, None)
        self.assertEqual(param.error, None)
        self.assertEqual(param.units, None)

    def test_load_missing_kinetic_law_information_from_tsv_error_1(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backup=False, load_content=False)
        session = src.session

        src.load_kinetic_laws([3962])
        src.load_compounds()
        src.load_missing_kinetic_law_information_from_tsv([3962])

        p = session.query(Parameter).filter_by(observed_name='kcat').first()
        self.assertEqual(p.observed_value, 0.9)
        self.assertEqual(p.observed_error, 0.6)
        self.assertEqual(p.observed_units, 's^(-1)')

        p = session.query(Parameter).filter_by(observed_name='app_kcat_MEP').first()
        self.assertEqual(p.observed_value, 0.54)
        self.assertEqual(p.observed_error, 0.36)
        self.assertEqual(p.observed_units, 's^(-1)')

    def test_load_missing_kinetic_law_information_from_tsv_error_2(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backup=False, load_content=False)
        session = src.session

        src.load_kinetic_laws([443])
        src.load_compounds()
        src.load_missing_kinetic_law_information_from_tsv([443])

        p = session.query(Parameter).filter_by(observed_name='kcat1').first()
        self.assertEqual(p.observed_value, 28.3333333)
        self.assertEqual(p.observed_error, 1.66667)
        self.assertEqual(p.observed_units, 's^(-1)')

        p = session.query(Parameter).filter_by(observed_name='kcat2').first()
        self.assertEqual(p.observed_value, 27.6666667)
        self.assertEqual(p.observed_error, 0.5)
        self.assertEqual(p.observed_units, 's^(-1)')

    def test_normalize_parameter_value(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backup=False, load_content=False)

        self.assertEqual(src.normalize_parameter_value('k_d', 282, 0.25, 0.15, None),
                         (None, None, None, None, None))
        self.assertEqual(src.normalize_parameter_value('k_cat', 25, None, None, 's^(-1)'),
                         (None, None, None, None, None))

        self.assertEqual(src.normalize_parameter_value('k_cat', 25, 0.25, 0.15, None),
                         (None, None, None, None, None))
        self.assertEqual(src.normalize_parameter_value('k_cat', 25, 0.25, 0.15, 's^(-1)'),
                         ('k_cat', 25, 0.25, 0.15, 's^(-1)'))
        self.assertEqual(src.normalize_parameter_value('k_cat', 25, 0.25, 0.15, 'katal_base'),
                         ('k_cat', 25, 0.25 * scipy.constants.Avogadro, 0.15 * scipy.constants.Avogadro, 's^(-1)'))
        self.assertEqual(src.normalize_parameter_value('k_cat', 25, 0.25, None, 'katal_base'),
                         ('k_cat', 25, 0.25 * scipy.constants.Avogadro, None, 's^(-1)'))
        self.assertEqual(src.normalize_parameter_value('k_cat', 25, 0.25, None, 'mol*s^(-1)*g^(-1)'),
                         (None, None, None, None, None))

        self.assertEqual(src.normalize_parameter_value('v_max', 186, 0.25, None, None),
                         (None, None, None, None, None))
        self.assertEqual(src.normalize_parameter_value('v_max', 186, 0.25, None, 'mol*s^(-1)*g^(-1)'),
                         ('v_max', 186, 0.25, None, 'mol*s^(-1)*g^(-1)'))
        self.assertEqual(src.normalize_parameter_value('v_max', 186, 0.25, None, 's^(-1)'),
                         ('k_cat', 25, 0.25, None, 's^(-1)'))

        self.assertRaises(ValueError, src.normalize_parameter_value, 'k_cat', 25, 0.25, 0.15, 'm')

    def test_infer_compound_structures_from_names(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backup=False, load_content=False)

        session = src.session

        compound_unknown = sabio_rk.Compound(name='Unknown')
        compound_no_pubchem = sabio_rk.Compound(name='a random string: adfkja;uvhayocbvadf')
        compound_one_pubchem = sabio_rk.Compound(name='water')
        compound_multiple_pubchem = sabio_rk.Compound(name='Phenyllactate')
        session.add(compound_unknown)
        session.add(compound_no_pubchem)
        session.add(compound_one_pubchem)
        session.add(compound_multiple_pubchem)

        compounds = [compound_unknown, compound_no_pubchem, compound_one_pubchem, compound_multiple_pubchem]
        src.infer_compound_structures_from_names(compounds)

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

    def test_calc_inchi_formula_connectivity(self):
        s = sabio_rk.CompoundStructure(format='smiles', value='[H]O[H]')
        s.calc_inchi_formula_connectivity()
        self.assertEqual(s._value_inchi, 'InChI=1S/H2O/h1H2')
        self.assertEqual(s._value_inchi_formula_connectivity, 'H2O')

        s = sabio_rk.CompoundStructure(format='inchi', value='InChI=1S/H2O/h1H2')
        s.calc_inchi_formula_connectivity()
        self.assertEqual(s._value_inchi, 'InChI=1S/H2O/h1H2')
        self.assertEqual(s._value_inchi_formula_connectivity, 'H2O')

        s = sabio_rk.CompoundStructure(
            format='inchi', value='InChI=1S/C9H10O3/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8,10H,6H2,(H,11,12)/t8-/m1/s1')
        s.calc_inchi_formula_connectivity()
        self.assertEqual(s._value_inchi, 'InChI=1S/C9H10O3/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8,10H,6H2,(H,11,12)/t8-/m1/s1')
        self.assertEqual(s._value_inchi_formula_connectivity, 'C9H10O3/c10-8(9(11)12)6-7-4-2-1-3-5-7')

    def test_load_content(self):
        # get some kinetic laws
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backup=False,  load_content=True,
                               max_entries=9, webservice_batch_size=1, excel_batch_size=3, verbose=True)
        self.assertTrue(os.path.isfile(src.filename))
        self.assertEqual(src.session.query(KineticLaw).count(), src.max_entries)

        # get some more
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backup=False, load_content=True,
                               max_entries=18, webservice_batch_size=1, excel_batch_size=3, verbose=True)
        self.assertTrue(os.path.isfile(src.filename))
        self.assertEqual(src.session.query(KineticLaw).count(), src.max_entries)

    def test_load_content_commit_intermediate_results(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backup=False,  load_content=True,
                               max_entries=9, commit_intermediate_results=True,
                               webservice_batch_size=1, excel_batch_size=3, verbose=True)
        self.assertEqual(src.session.query(KineticLaw).count(), src.max_entries)


class TestBackupAndInstall(unittest.TestCase):

    def setUp(self):
        self.cache_dirname_1 = tempfile.mkdtemp()
        self.cache_dirname_2 = tempfile.mkdtemp()
        self.cache_dirname_3 = tempfile.mkdtemp()
        self.cache_dirname_4 = tempfile.mkdtemp()
        self.cache_dirname_5 = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.cache_dirname_1)
        shutil.rmtree(self.cache_dirname_2)
        shutil.rmtree(self.cache_dirname_3)
        shutil.rmtree(self.cache_dirname_4)
        shutil.rmtree(self.cache_dirname_5)

    def test(self):
        # backup and download
        src = sabio_rk.SabioRk(name='test.sabio_rk', cache_dirname=self.cache_dirname_1, download_backup=False, load_content=True,
                               max_entries=2, webservice_batch_size=1, excel_batch_size=10, verbose=True, download_request_backup=True)
        src.upload_backup()

        # download
        src2 = sabio_rk.SabioRk(name='test.sabio_rk', cache_dirname=self.cache_dirname_2, download_backup=True, load_content=False,
                                max_entries=2, webservice_batch_size=1, excel_batch_size=10, verbose=True, download_request_backup=True)
        self.assertTrue(os.path.isfile(src2.filename))
        self.assertGreater(os.stat(src2.requests_cache_filename).st_size, 1e6)
        self.assertEqual(src2.session.query(KineticLaw).count(), 2)

        # setup with download and update
        src3 = sabio_rk.SabioRk(name='test.sabio_rk', cache_dirname=self.cache_dirname_3, download_backup=True, load_content=True,
                                max_entries=4, webservice_batch_size=1, excel_batch_size=10, verbose=True, download_request_backup=True)
        self.assertEqual(src3.session.query(KineticLaw).count(), 4)

        # setup with update
        src4 = sabio_rk.SabioRk(name='test.sabio_rk', cache_dirname=self.cache_dirname_4, download_backup=False, load_content=True,
                                max_entries=4, webservice_batch_size=1, excel_batch_size=10, verbose=True, download_request_backup=True)
        self.assertEqual(src4.session.query(KineticLaw).count(), 4)

        # no download or update
        src5 = sabio_rk.SabioRk(name='test.sabio_rk', cache_dirname=self.cache_dirname_5, download_backup=False, load_content=False,
                                max_entries=2, webservice_batch_size=1, excel_batch_size=10, verbose=True, download_request_backup=True)
        self.assertEqual(src5.session.query(KineticLaw).count(), 0)


class TestStats(unittest.TestCase):

    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_export_stats(self):
        src = sabio_rk.SabioRk()
        stats = src.calc_stats()
        filename = os.path.join(self.dirname, 'stats.xlsx')
        src.export_stats(stats, filename=filename)
        self.assertTrue(os.path.isfile(filename))


class TestAll(unittest.TestCase):

    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    @unittest.skip('Skip this test because it is long')
    def test_update_all_and_backup(self):
        src = sabio_rk.SabioRk(download_backup=False, load_content=True, clear_content=True,
                               verbose=True, clear_requests_cache=False, max_entries=float('inf'))
        self.assertGreaterEqual(src.session.query(KineticLaw).count(), 55000)

        src.upload_backup()

    def test_download_full_database(self):
        src = sabio_rk.SabioRk(cache_dirname=self.dirname, download_backup=True, load_content=False,
                               verbose=True)
        self.assertGreaterEqual(src.session.query(KineticLaw).count(), 55000)
