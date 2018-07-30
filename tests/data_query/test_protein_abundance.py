""" Tests of common_schemy protein abundance queries

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-09-19
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model
from kinetic_datanator.data_query import protein_abundance
from kinetic_datanator.core import models, common_schema
import tempfile
import shutil
import unittest

class TestProteinAbundanceQuery(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        flk = common_schema.FlaskCommonSchema(cache_dirname=self.cache_dirname)

        self.protein_P00323 = flk.session.query(models.ProteinSubunit).filter_by(uniprot_id = 'P00323').first()
        self.protein_Q42025 = flk.session.query(models.ProteinSubunit).filter_by(uniprot_id = 'Q42025').first()
        self.q = protein_abundance.ProteinAbundanceQuery(cache_dirname=self.cache_dirname)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)

    def test_filter_observed_results(self):

        vals = self.q.get_observed_result(self.protein_P00323)
        self.assertEqual(set(items.value for items in vals), set([1003.0, 1100.0, 2323.0, 1336.0, 639.0]))
        self.assertEqual(set(items.units for items in vals), set(['PPM']))

        for items in vals:
            self.assertEqual(items.observable.specie.name, self.protein_P00323.subunit_name)
            self.assertEqual(items.observable.specie.uniprot_id, self.protein_P00323.uniprot_id)
            self.assertEqual(items.observable.specie.entrez_id, self.protein_P00323.entrez_id)
            self.assertEqual(items.observable.specie.gene_name, self.protein_P00323.gene_name)
            self.assertEqual(items.observable.specie.length, self.protein_P00323.length)
            self.assertEqual(items.observable.specie.mass, self.protein_P00323.mass)
            self.assertEqual(items.observable.specie.sequence, self.protein_P00323.canonical_sequence)
            if items.observable.specie.cross_references[0].id == '882/882-Desulfo_Form_Exp_SC_zhang_2006.txt':
                self.assertEqual(items.value, 1003.0)
            self.assertIn(items.metadata.genetics.taxon, ['Desulfovibrio vulgaris str. Hildenborough','Synechocystis sp. PCC 6803'])
            self.assertEqual(set([ref.namespace for ref in items.metadata.cross_references]), set(['url']))


        vals = self.q.get_observed_result(self.protein_Q42025)
        self.assertEqual(set(items.value for items in vals), set([240.0, 140.0, 358.0]))
        self.assertEqual(set(items.units for items in vals), set(['PPM']))

        for items in vals:
            self.assertEqual(items.observable.specie.name, self.protein_Q42025.subunit_name)
            self.assertEqual(items.observable.specie.uniprot_id, self.protein_Q42025.uniprot_id)
            self.assertEqual(items.observable.specie.entrez_id, self.protein_Q42025.entrez_id)
            self.assertEqual(items.observable.specie.gene_name, self.protein_Q42025.gene_name)
            self.assertEqual(items.observable.specie.length, self.protein_Q42025.length)
            self.assertEqual(items.observable.specie.mass, self.protein_Q42025.mass)
            self.assertEqual(items.observable.specie.sequence, self.protein_Q42025.canonical_sequence)
            if items.observable.specie.cross_references[0].id == '3702/3702-2.1ArabidopsisLeaf_waterDefict_SOW.txt':
                self.assertEqual(items.value, 140.0)



    def test_get_abundance_by_uniprot(self):
        uniprot = self.protein_P00323.uniprot_id
        abundances = self.q.get_abundance_by_uniprot(uniprot).filter(
            models.AbundanceData.pax_load.in_([1, 2])).all()
        self.assertEqual(set(c.abundance for c in abundances),
                         set([1003.0, 1336.0]))

    def test_get_abundance_by_gene_name(self):
        gene_name = 'rplO'
        abundances = self.q.get_abundance_by_gene_name(gene_name).filter(
            models.AbundanceData.pax_load.in_([1, 2])).all()
        self.assertEqual(set(c.abundance for c in abundances),
                         set([1027.0, 1861.0]))

    def test_get_abundance_by_sequence(self):

        sequence = self.protein_P00323.canonical_sequence
        abundances = self.q.get_abundance_by_sequence(sequence).filter(
            models.AbundanceData.pax_load.in_([1, 2])).all()
        self.assertEqual(set(c.abundance for c in abundances),
                         set([1003.0, 1336.0]))

    def test_get_abundance_by_entrez(self):
        entrez_id = self.protein_P00323.entrez_id
        abundances = self.q.get_abundance_by_entrez(entrez_id).filter(
            models.AbundanceData.pax_load.in_([1, 2])).all()
        self.assertEqual(set(c.abundance for c in abundances),
                         set([1003.0, 1336.0]))

    def test_get_abundance_by_mass(self):
        mass = self.protein_P00323.mass
        abundances = self.q.get_abundance_by_mass(mass).filter(
            models.AbundanceData.pax_load.in_([1, 2])).all()
        self.assertEqual(set(c.abundance for c in abundances),
                         set([1003.0, 1336.0]))

    def test_get_abundance_by_length(self):
        length = self.protein_P00323.length
        abundances = self.q.get_abundance_by_length(length).filter(
            models.AbundanceData.pax_load.in_([1, 2])).all()
        self.assertEqual(set(c.abundance for c in abundances),
                         set([1003.0, 1336.0, 1027.0, 1861.0]))
