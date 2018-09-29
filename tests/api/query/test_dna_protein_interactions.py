""" Tests of DNA and Protein interaction queries

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-09-27
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model
from kinetic_datanator.api.query import dna_protein_interactions as dpi
from kinetic_datanator.core import models, common_schema
import unittest
from Bio import motifs
import tempfile
import shutil
import csv


class TestProteintoDNAInteractionQuery(unittest.TestCase):
    """
    Tests for Protein to DNA interactions

    Case: You have a protein and want to find the DNA binding site

    """

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        flk = common_schema.CommonSchema(cache_dirname=cls.cache_dirname)

        cls.arnt  = flk.session.query(models.ProteinSubunit).filter_by(uniprot_id = 'P53762').first()
        cls.q = dpi.ProteintoDNAInteractionQuery(cache_dirname=cls.cache_dirname)

    def test_get_observed_result(self):

        result = self.q.get_observed_result(self.arnt)

        for version in result:
            self.assertEqual(version.specie.sequence, 'CACGTG')
            self.assertEqual(version.specie.cross_references[0].id, '7592839')
            self.assertEqual(version.metadata.genetics.taxon, 'Mus musculus')
            self.assertEqual(version.metadata.method.name,'SELEX')
            break


    def test_get_DNA_by_protein(self):

        position = self.q.get_DNA_by_protein(self.arnt)
        for version in position:
            self.assertEqual(set(c.frequency_a for c in version), set([0,19,4]))
            self.assertEqual(set(c.position for c in version), set([1, 2, 3, 4, 5, 6]))

class TestDNAtoProteinInteractionQuery(unittest.TestCase):
    """
    Tests for DNA to Protein interactions

    Case: You have a DNA segment and want to find binding Protein

    """
    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.dna_segment1 = data_model.DnaSpecie(sequence = 'CCTTTGTT')
        cls.dna_segment2 = data_model.DnaSpecie(sequence = 'AAGGTCAA')

        cls.q = dpi.DNAtoProteinInteractionQuery(cache_dirname=cls.cache_dirname)

    def test_get_observed_result(self):
        
        observe = self.q.get_observed_result(self.dna_segment2)
        self.assertEqual(set(c.specie.gene_name for c in observe), set(['NR4A2', 'TRP(MYB) class']))


    def test_get_protein_by_binding_matrix(self):

        query = self.q.get_protein_by_DNA_sequence(self.dna_segment1.sequence)

        self.assertEqual(set(c[0].subunit_name for c in query), set(['pan', 'Sox2', 'DOF5.6']))
        self.assertEqual(set(c[1] for c in query), set([0,0,-8,-8]))
