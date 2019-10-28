# -*- coding: utf-8 -*-

""" Tests of sabio_rk

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-26
:Copyright: 2017, Karr Lab
:License: MIT
"""

from datanator.data_source import sabio_rk
from datanator.data_source.sabio_rk import (Entry, Compartment, Compound, Enzyme,
                                            ReactionParticipant, KineticLaw, Parameter, Resource)
from datanator.util import warning_util
import capturer
import datetime
import ftputil
import datanator.config
import math
import mock
import numpy
import os
import quilt
import scipy.constants
import shutil
import sqlalchemy
import sys
import tempfile
import unittest
import wc_utils.config
import wc_utils.quilt

warning_util.disable_warnings()


class TestDownloader(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_load_kinetic_law_ids(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backups=False, load_content=False, max_entries=10,
                               webservice_batch_size=1, excel_batch_size=5, verbose=True)
        ids = src.load_kinetic_law_ids()
        self.assertEqual(ids[0:10], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        self.assertGreater(len(ids), 55000)

    def test_load_kinetic_laws_and_compounds(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backups=False, load_content=False, verbose=True)
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
        self.assertLess((datetime.datetime.utcnow() - c.created).total_seconds(), 3600)
        self.assertLess((datetime.datetime.utcnow() - c.modified).total_seconds(), 3600)
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
        self.assertEqual(e.subunits, [])
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
        self.assertEqual(p.observed_value, 1.8E-8)
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
        self.assertNotIn(dict(namespace='None', id='View all entries for compound H2O'), xrs)
        self.assertIn(dict(namespace='chebi', id='CHEBI:15377'), xrs)
        self.assertIn(dict(namespace='kegg.compound', id='C00001'), xrs)
        self.assertIn(dict(namespace='pubchem.substance', id='3303'), xrs)

        self.assertEqual(c.created, h20_created)
        self.assertIsInstance(c.modified, datetime.datetime)
        self.assertLess((datetime.datetime.utcnow() - c.modified).total_seconds(), 3600)

    def test_load_and_update_kinetic_laws_and_compounds(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backups=False, load_content=False, verbose=True)
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
        self.assertNotIn(dict(namespace='None', id='View all entries for compound H2O'), xrs)
        self.assertIn(dict(namespace='chebi', id='CHEBI:15377'), xrs)
        self.assertIn(dict(namespace='kegg.compound', id='C00001'), xrs)
        self.assertIn(dict(namespace='pubchem.substance', id='3303'), xrs)

        self.assertIsInstance(c.created, datetime.datetime)
        self.assertIsInstance(c.modified, datetime.datetime)
        self.assertLess((datetime.datetime.utcnow() - c.created).total_seconds(), 3600)
        self.assertLess((datetime.datetime.utcnow() - c.modified).total_seconds(), 3600)
        h20_created = c.created

        c = session.query(Compound).filter_by(id=2562).first()
        self.assertTrue('Peptide' in c.name)
        xrs = [dict(namespace=xr.namespace, id=xr.id) for xr in c.cross_references]
        self.assertNotIn(dict(namespace='None', id='View all entries for compound H2O'), xrs)
        self.assertIn(dict(namespace='chebi', id='CHEBI:16670'), xrs)
        self.assertIn(dict(namespace='kegg.compound', id='C00012'), xrs)
        self.assertIn(dict(namespace='pubchem.substance', id='3314'), xrs)

        c = session.query(Compound).filter_by(id=20035).first()
        self.assertEqual(c.name, 'N-(Carbobenzoxy)-Leu-Leu-Phe-trifluoromethylketone')
        xrs = {xr.namespace: xr.id for xr in c.cross_references}
        self.assertNotIn(None, xrs)
        self.assertNotIn('None', xrs)

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
        self.assertEqual(p.observed_error, None)
        self.assertEqual(p.observed_units, 'M')
        self.assertEqual(p.name, 'k_i')
        self.assertEqual(p.type, 261)
        self.assertEqual(p.value, 1.8E-8)
        self.assertEqual(p.error, None)
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
        self.assertNotIn(dict(namespace='None', id='View all entries for compound H2O'), xrs)
        self.assertIn(dict(namespace='chebi', id='CHEBI:15377'), xrs)
        self.assertIn(dict(namespace='kegg.compound', id='C00001'), xrs)
        self.assertIn(dict(namespace='pubchem.substance', id='3303'), xrs)

        self.assertEqual(c.created, h20_created)
        self.assertIsInstance(c.modified, datetime.datetime)
        self.assertLess((datetime.datetime.utcnow() - c.modified).total_seconds(), 3600)

    def test_load_kinetic_laws_multiple(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backups=False, load_content=False, verbose=True)
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
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backups=False, load_content=False, verbose=True)
        src.load_kinetic_laws([10054])

        session = src.session

        l = session.query(KineticLaw).filter_by(id=10054).first()
        self.assertEqual(l.enzyme_type, 'Modifier-Catalyst')

    def test_load_kinetic_laws_with_same_reaction(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backups=False, load_content=False, verbose=True)
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
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backups=False, load_content=False)
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

    def test_parse_complex_subunit_structure(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backups=False, load_content=False, verbose=True)

        self.assertEqual(src.parse_complex_subunit_structure((
            '(<a href="http://www.uniprot.org/uniprot/Q59669" target="_blank">Q59669</a>)'
        )), {'Q59669': 1})

        self.assertEqual(src.parse_complex_subunit_structure((
            '(<a href="http://www.uniprot.org/uniprot/Q59669" target="_blank">Q59669</a>)*2'
        )), {'Q59669': 2})

        self.assertEqual(src.parse_complex_subunit_structure((
            '('
            '(<a href="http://www.uniprot.org/uniprot/Q59669" target="_blank">Q59669</a>)'
            '(<a href="http://www.uniprot.org/uniprot/Q59670" target="_blank">Q59670</a>)'
            ')'
        )), {'Q59669': 1, 'Q59670': 1})

        self.assertEqual(src.parse_complex_subunit_structure((
            '('
            '(<a href="http://www.uniprot.org/uniprot/Q59669" target="_blank">Q59669</a>)*2'
            '(<a href="http://www.uniprot.org/uniprot/Q59670" target="_blank">Q59670</a>)*3'
            ')*4'
        )), {'Q59669': 8, 'Q59670': 12})

        self.assertEqual(src.parse_complex_subunit_structure((
            '('
            '(<a href="http://www.uniprot.org/uniprot/Q59669" target="_blank">Q59669</a>)'
            '(<a href="http://www.uniprot.org/uniprot/Q59670" target="_blank">Q59670</a>)*2'
            ')*3'
        )), {'Q59669': 3, 'Q59670': 6})

        self.assertEqual(src.parse_complex_subunit_structure((
            '<a href="http://www.uniprot.org/uniprot/P09219" target=_blank>P09219</a>; '
            '<a href="http://www.uniprot.org/uniprot/P07677" target=_blank>P07677</a>; '
        )), {'P09219': 1, 'P07677': 1})

        self.assertEqual(src.parse_complex_subunit_structure((
            '<a href="http://www.uniprot.org/uniprot/P09219" target=_blank>P09219</a>; '
            '<a href="http://www.uniprot.org/uniprot/P07677" target=_blank>P07677</a>; '
        )), {'P09219': 1, 'P07677': 1})

        self.assertEqual(src.parse_complex_subunit_structure((
            '(<a href="http://www.uniprot.org/uniprot/P19112" target=_blank>P19112</a>)*4; '
            '<a href="http://www.uniprot.org/uniprot/Q9Z1N1" target=_blank>Q9Z1N1</a>; '
        )), {'P19112': 4, 'Q9Z1N1': 1})

        self.assertEqual(src.parse_complex_subunit_structure((
            '((<a href="http://www.uniprot.org/uniprot/P16924" target="_blank">P16924</a>)*2'
            '(<a href="http://www.uniprot.org/uniprot/P09102" target="_blank">P09102</a>)*2); '
            '((<a href="http://www.uniprot.org/uniprot/Q5ZLK5" target="_blank">Q5ZLK5</a>)*2'
            '(<a href="http://www.uniprot.org/uniprot/P09102" target="_blank">P09102</a>)*2);'
        )), {'P16924': 2, 'P09102': 4, 'Q5ZLK5': 2})

        self.assertEqual(src.parse_complex_subunit_structure((
            '((<a href="http://www.uniprot.org/uniprot/Q03393" target=_blank>Q03393</a>)*3)*2); '
        )), {'Q03393': 6})

    def test_loading_enzymes(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backups=False, load_content=False, verbose=True)
        ids = [11021, 2139, 1645]
        src.load_kinetic_laws(ids)
        src.load_compounds()
        src.load_missing_kinetic_law_information_from_tsv(ids)
        src.load_missing_enzyme_information_from_html(ids)
        src.calc_enzyme_molecular_weights(src.session.query(Enzyme).all())
        src.normalize_kinetic_laws(ids)

        # 11021, one subunit
        law = src.session.query(KineticLaw).filter_by(id=11021).first()
        self.assertEqual(len(law.enzyme.subunits), 0)
        self.assertEqual(len(law.enzyme.subunits[0].cross_references), 1)
        self.assertEqual(law.enzyme.subunits[0].cross_references[0].namespace, 'uniprot')
        self.assertEqual(law.enzyme.subunits[0].cross_references[0].id, 'P19631')
        self.assertEqual(law.enzyme.subunits[0].coefficient, 2)
        self.assertEqual(law.enzyme.subunits[0].sequence, (
            'STAGKVIKCKAAVLWEANKPFSLEEVEVAPPKAHEVRIKIVATGICRSDDHVVTGALAMP'
            'FPVILGHEAAGVVESVGEKVTLLKPGDAVIPLFVPQCGECRSCLSTKGNLCIKNDLSSSP'
            'TGLMADGTTRFTCKGKAIHHFIGTSTFTEYTVVHETAAAKIDSAAPLEKVCLIGCGFSTG'
            'YGAVLQTAKVEPGSTCAVFGLGGVGLSVVMGCKAAGASRIIAIDINKDKFAKAKELGATE'
            'CVNPKDFKKPIHEVLTEMTGKGVDYSFEVIGRIETMTEALASCHYNYGVSVIVGVPPAAQ'
            'KISFDPMLIFSGRTWKGSVFGGWKSKDAVPKLVADYMKKKFVLDPLITHTLPFTKINEGF'
            'DLLRTGKSIRTVLVL'
        ))
        numpy.testing.assert_approx_equal(law.enzyme.subunits[0].molecular_weight, 39810, significant=3)
        self.assertEqual(law.enzyme.cross_references, [])

        # 2139, two subunits
        law = src.session.query(KineticLaw).filter_by(id=2139).first()
        self.assertEqual(len(law.enzyme.subunits), 2)

        subunit = next(subunit for subunit in law.enzyme.subunits if subunit.cross_references[0].id == 'P07677')
        self.assertEqual(len(subunit.cross_references), 1)
        self.assertEqual(subunit.cross_references[0].namespace, 'uniprot')
        self.assertEqual(subunit.cross_references[0].id, 'P07677')
        self.assertEqual(subunit.coefficient, 1)
        self.assertEqual(subunit.sequence, (
            'MTRGRVIQVMGPVVDVKFENGHLPAIYNALKIQHKARNENEVDIDLTLEVALHLGDDTVR'
            'TIAMASTDGLIRGMEVIDTGAPISVPVGQVTLGRVFNVLGEPIDLEGDIPADARRDPIHR'
            'PAPKFEELATEVEILETGIKVVDLLAPYIKGGKIGLFGGAGVGKTVLIQELIHNIAQEHG'
            'GISVFAGVGERTREGNDLYHEMKDSGVISKTAMVFGQMNEPPGARMRVALTGLTMAEYFR'
            'DEQGQDGLLFIDNIFRFTQAGSEVSALLGRMPSAIGYQPTLATEMGQLQERITSTAKGSI'
            'TSIQAIYVPADDYTDPAPATTFSHLDATTNLERKLAEMGIYPAVDPLVSTSRALAPEIVG'
            'EEHYQVARKVQQTLERYKELQDIIAILGMDELSDEDKLVVHRARRIQFFLSQNFHVAEQF'
            'TGQPGSYVPVKETVRGFKEILEGKYDHLPEDRFRLVGRIEEVVEKAKAMGVEV'
        ))
        numpy.testing.assert_approx_equal(subunit.molecular_weight, 51950, significant=3)

        subunit = next(subunit for subunit in law.enzyme.subunits if subunit.cross_references[0].id == 'P09219')
        self.assertEqual(len(subunit.cross_references), 1)
        self.assertEqual(subunit.cross_references[0].namespace, 'uniprot')
        self.assertEqual(subunit.cross_references[0].id, 'P09219')
        self.assertEqual(subunit.coefficient, 1)
        self.assertEqual(subunit.sequence, (
            'MSIRAEEISALIKQQIENYESQIQVSDVGTVIQVGDGIARAHGLDNVMSGEAVEFANAVM'
            'GMALNLEENNVGIVILGPYTGIKEGDEVRRTGRIMEVPVGETLIGRVVNPLGQPVDGLGP'
            'VETTETRPIESRAPGVMDRRSVHEPLQTGIKAIDALVPIGRGQRELIIGDRQTGKTSVAI'
            'DTIINQKDQNMICIYVAIGQKESTVATVVETLAKHGAPDYTIVVTASASQPAPLLFLAPY'
            'AGVAMGEYFMIMGKHVLVVIDDLSKQAAAYRQLSLLLRRPPGREAYPGDIFYLHSRLLER'
            'AAKLSDAKGGGSLTALPFVETQAGDISAYIPTNVISITDGQIFLQSDLFFSGVRPAINAG'
            'LSVSRVGGAAQIKAMKKVAGTLRLDLAAYRELEAFAQFGSDLDKATQANVARGARTVEVL'
            'KQDLHQPIPVEKQVLIIYALTRGFLDDIPVEDVRRFEKEFYLWLDQNGQHLLEHIRTTKD'
            'LPNEDDLNQAIEAFKKTFVVSQ'
        ))
        numpy.testing.assert_approx_equal(subunit.molecular_weight, 54600, significant=3)

        self.assertEqual(law.enzyme.cross_references, [])

        # 1645, one subunit not in SBML
        law = src.session.query(KineticLaw).filter_by(id=1645).first()
        self.assertEqual(len(law.enzyme.subunits), 1)
        self.assertEqual(len(law.enzyme.subunits[0].cross_references), 1)
        self.assertEqual(law.enzyme.subunits[0].cross_references[0].namespace, 'uniprot')
        self.assertEqual(law.enzyme.subunits[0].cross_references[0].id, 'Q70GK9')
        self.assertEqual(law.enzyme.subunits[0].coefficient, 1)
        self.assertEqual(law.enzyme.cross_references, [])

    def test_load_kinetic_laws_check_units(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backups=False, load_content=False, verbose=True)
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
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backups=False, load_content=False, verbose=True)
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
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backups=False, load_content=False, verbose=True)
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
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backups=False, load_content=False, verbose=True)

        self.assertEqual(src.normalize_parameter_value('k_d', 282, 0.25, 0.15, None, None),
                         (None, None, None, None, None))
        self.assertEqual(src.normalize_parameter_value('k_cat', 25, None, None, 's^(-1)', None),
                         (None, None, None, None, None))

        self.assertEqual(src.normalize_parameter_value('k_cat', 25, 0.25, 0.15, None, None),
                         (None, None, None, None, None))
        self.assertEqual(src.normalize_parameter_value('k_cat', 25, 0.25, 0.15, 's^(-1)', None),
                         ('k_cat', 25, 0.25, 0.15, 's^(-1)'))
        self.assertEqual(src.normalize_parameter_value('k_cat', 25, 0.25, 0.15, 'katal_base', None),
                         (None, None, None, None, None))
        self.assertEqual(src.normalize_parameter_value('k_cat', 25, 0.25, None, 'katal_base', None),
                         (None, None, None, None, None))
        self.assertEqual(src.normalize_parameter_value('k_cat', 25, 0.25, None, 'mol*s^(-1)*g^(-1)', None),
                         (None, None, None, None, None))

        self.assertEqual(src.normalize_parameter_value('v_max', 186, 0.25, None, None, 500),
                         (None, None, None, None, None))
        self.assertEqual(src.normalize_parameter_value('v_max', 186, 0.25, None, 'mol*s^(-1)*g^(-1)', 500),
                         ('k_cat', 25, 0.25 * 500, None, 's^(-1)'))
        self.assertEqual(src.normalize_parameter_value('v_max', 186, 0.25, None, 'mol*s^(-1)*g^(-1)', None),
                         ('v_max', 186, 0.25, None, 'mol*s^(-1)*g^(-1)'))
        self.assertEqual(src.normalize_parameter_value('v_max', 186, 0.25, 0.15, 'mol*s^(-1)*g^(-1)', 500),
                         ('k_cat', 25, 0.25 * 500, 0.15 * 500, 's^(-1)'))
        self.assertEqual(src.normalize_parameter_value('v_max', 186, 0.25, 0.15, 'mol*s^(-1)*g^(-1)', None),
                         ('v_max', 186, 0.25, 0.15, 'mol*s^(-1)*g^(-1)'))
        self.assertEqual(src.normalize_parameter_value('v_max', 186, 0.25, None, 's^(-1)', 500),
                         ('k_cat', 25, 0.25, None, 's^(-1)'))
        self.assertEqual(src.normalize_parameter_value('v_max', 186, 0.25, None, 's^(-1)', None),
                         ('k_cat', 25, 0.25, None, 's^(-1)'))

        self.assertRaises(ValueError, src.normalize_parameter_value, 'k_cat', 25, 0.25, 0.15, 'm', None)

    def test_enzyme_parameter(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname,
                               download_backups=False,
                               load_content=False, verbose=True)
        session = src.session

        src.load_kinetic_laws([213])
        src.load_compounds()
        src.load_missing_kinetic_law_information_from_tsv([213])

        law = session.query(KineticLaw).filter_by(id=213).first()

        params = list(filter(lambda param: isinstance(param.enzyme, sabio_rk.Enzyme), law.parameters))
        self.assertEqual(len(params), 1)
        self.assertEqual(params[0].enzyme.id, 1000)
        self.assertEqual(params[0].enzyme.name, 'inorganic diphosphatase')

    def test_infer_compound_structures_from_names(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backups=False, load_content=False)

        session = src.session

        compound_unknown = sabio_rk.Compound(name='Unknown')
        compound_no_pubchem = sabio_rk.Compound(name='a random string: adfkja;uvhayocbvadf')
        compound_one_pubchem = sabio_rk.Compound(name='water')
        compound_multiple_pubchem = sabio_rk.Compound(name='Aspartate')
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
            ('pubchem.compound', '5960'),
            ('pubchem.compound', '5460541'),
            ('pubchem.compound', '5460294'),
        ]))
        self.assertEqual(set([(s.format, s.value) for s in compound_multiple_pubchem.structures]), set([
            ('inchi', 'InChI=1S/C4H7NO4/c5-2(4(8)9)1-3(6)7/h2H,1,5H2,(H,6,7)(H,8,9)/t2-/m0/s1'),
            ('inchi', 'InChI=1S/C4H7NO4/c5-2(4(8)9)1-3(6)7/h2H,1,5H2,(H,6,7)(H,8,9)/p-2/t2-/m0/s1'),
            ('inchi', 'InChI=1S/C4H7NO4/c5-2(4(8)9)1-3(6)7/h2H,1,5H2,(H,6,7)(H,8,9)/p-1/t2-/m0/s1'),
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

    def test_full_kinetic_laws(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname,
                               download_backups=False,
                               load_content=False, verbose=True)
        session = src.session

        src.load_kinetic_laws([23637])
        src.load_compounds()
        src.load_missing_kinetic_law_information_from_tsv([23637])

        law = session.query(KineticLaw).filter_by(id=23637).first()
        self.assertEqual(law.equation, 'V * S / (Km + S + S^2 / Ki)')

    def test_load_content(self):
        # get some kinetic laws
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backups=False,  load_content=True,
                               max_entries=9, webservice_batch_size=1, excel_batch_size=3, verbose=True)
        self.assertTrue(os.path.isfile(src.filename))
        self.assertEqual(src.session.query(KineticLaw).count(), src.max_entries)

        # get some more
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backups=False, load_content=True,
                               max_entries=18, webservice_batch_size=1, excel_batch_size=3, verbose=True)
        self.assertTrue(os.path.isfile(src.filename))
        self.assertEqual(src.session.query(KineticLaw).count(), src.max_entries)

    def test_load_content_commit_intermediate_results(self):
        src = sabio_rk.SabioRk(cache_dirname=self.cache_dirname, download_backups=False,  load_content=True,
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

        self.owner = datanator.config.get_config()['datanator']['quilt']['owner']
        self.package = 'datanator_test__'
        self.owner_package = '{}/{}'.format(self.owner, self.package)
        self.token = wc_utils.config.get_config()['wc_utils']['quilt']['token']

        self.tmp_dirname = tempfile.mkdtemp()

        self.delete_test_package_locally()
        self.delete_test_package_remotely()

    def tearDown(self):
        shutil.rmtree(self.cache_dirname_1)
        shutil.rmtree(self.cache_dirname_2)
        shutil.rmtree(self.cache_dirname_3)
        shutil.rmtree(self.cache_dirname_4)
        shutil.rmtree(self.cache_dirname_5)

        shutil.rmtree(self.tmp_dirname)

        self.delete_test_package_locally()
        self.delete_test_package_remotely()

    def delete_test_package_locally(self):
        # remove local package
        with capturer.CaptureOutput(relay=False):
            quilt.rm('{}/{}'.format(self.owner, self.package), force=True)

    def delete_test_package_remotely(self):
        # delete package from Quilt server, if it exists
        try:
            with capturer.CaptureOutput(relay=False):
                quilt.access_list(self.owner_package)
                with mock.patch('quilt.tools.command.input', return_value=self.owner_package):
                    quilt.delete(self.owner_package)
        except quilt.tools.command.HTTPResponseException as err:
            if str(err) != 'Package {} does not exist'.format(self.owner_package):
                raise(err)

    def test(self):
        # create test Quilt package
        os.mkdir(os.path.join(self.tmp_dirname, 'up'))
        os.mkdir(os.path.join(self.tmp_dirname, 'up', 'subdir'))
        with open(os.path.join(self.tmp_dirname, 'up', 'subdir', 'README.md'), 'w') as file:
            file.write('# datanator_test__\n')
            file.write('Test package for datanator\n')

        manager = wc_utils.quilt.QuiltManager(os.path.join(self.tmp_dirname, 'up'), self.package, owner=self.owner)
        manager.upload()

        # backup and download
        src = sabio_rk.SabioRk(name='test.sabio_rk', cache_dirname=self.cache_dirname_1, download_backups=False, load_content=True,
                               max_entries=2, webservice_batch_size=1, excel_batch_size=10, verbose=True, download_request_backup=True,
                               quilt_package=self.package)
        src.upload_backups()

        # download
        src2 = sabio_rk.SabioRk(name='test.sabio_rk', cache_dirname=self.cache_dirname_2, download_backups=True, load_content=False,
                                max_entries=2, webservice_batch_size=1, excel_batch_size=10, verbose=True, download_request_backup=True,
                                quilt_package=self.package)
        self.assertTrue(os.path.isfile(src2.filename))
        self.assertGreater(os.stat(src2.requests_cache_filename).st_size, 1e6)
        self.assertEqual(src2.session.query(KineticLaw).count(), 2)

        # setup with download and update
        src3 = sabio_rk.SabioRk(name='test.sabio_rk', cache_dirname=self.cache_dirname_3, download_backups=True, load_content=True,
                                max_entries=4, webservice_batch_size=1, excel_batch_size=10, verbose=True, download_request_backup=True,
                                quilt_package=self.package)
        self.assertEqual(src3.session.query(KineticLaw).count(), 4)

        # setup with update
        src4 = sabio_rk.SabioRk(name='test.sabio_rk', cache_dirname=self.cache_dirname_4, download_backups=False, load_content=True,
                                max_entries=4, webservice_batch_size=1, excel_batch_size=10, verbose=True, download_request_backup=True,
                                quilt_package=self.package)
        self.assertEqual(src4.session.query(KineticLaw).count(), 4)

        # no download or update
        src5 = sabio_rk.SabioRk(name='test.sabio_rk', cache_dirname=self.cache_dirname_5, download_backups=False, load_content=False,
                                max_entries=2, webservice_batch_size=1, excel_batch_size=10, verbose=True, download_request_backup=True,
                                quilt_package=self.package)
        self.assertEqual(src5.session.query(KineticLaw).count(), 0)

        # check that README still in package
        os.mkdir(os.path.join(self.tmp_dirname, 'down'))
        manager = wc_utils.quilt.QuiltManager(os.path.join(self.tmp_dirname, 'down'), self.package, owner=self.owner)
        manager.download(system_path='subdir/README.md', sym_links=True)

        with open(os.path.join(self.tmp_dirname, 'down', 'subdir', 'README.md'), 'r') as file:
            self.assertEqual(file.readline(), '# datanator_test__\n')
            self.assertEqual(file.readline(), 'Test package for datanator\n')


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
        src = sabio_rk.SabioRk(download_backups=False, load_content=True, clear_content=True,
                               verbose=True, clear_requests_cache=False, max_entries=float('inf'))
        self.assertGreaterEqual(src.session.query(KineticLaw).count(), 55000)

        src.upload_backups()

    def test_download_full_database(self):
        src = sabio_rk.SabioRk(cache_dirname=self.dirname, download_backups=True, load_content=False,
                               verbose=True)
        self.assertGreaterEqual(src.session.query(KineticLaw).count(), 55000)
