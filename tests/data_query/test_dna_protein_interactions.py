""" Tests of DNA and Protein interaction queries

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-09-27
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model
from kinetic_datanator.data_query import dna_protein_interactions as dpi
from kinetic_datanator.core import models, flask_common_schema
import unittest
from Bio import motifs
import tempfile
import shutil
import csv


class TestProteintoDNAInteractionQueryGenerator(unittest.TestCase):
    """
    Tests for Protein to DNA interactions

    Case: You have a protein and want to find the DNA binding site

    """

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        flk = flask_common_schema.FlaskCommonSchema(cache_dirname=self.cache_dirname)

        self.arnt  = flk.session.query(models.ProteinSubunit).filter_by(uniprot_id = 'P53762').first()

    def test_filter_observed_values(self):
        q = dpi.ProteintoDNAInteractionQueryGenerator(cache_dirname=self.cache_dirname)

        observable = q.get_observed_result(self.arnt)

        self.assertEqual(observable[0].specie.sequence, 'CACGTG')
        self.assertEqual(observable[0].specie.cross_references[0].id, '7592839')


    def test_get_DNA_by_protein(self):
        q = dpi.ProteintoDNAInteractionQueryGenerator(cache_dirname=self.cache_dirname)

        position = q.get_DNA_by_protein(self.arnt)

        self.assertEqual(set(c.frequency_a for c in position[0]), set([0,19,4]))
        self.assertEqual(set(c.position for c in position[0]), set([1, 2, 3, 4, 5, 6]))

class TestDNAtoProteinInteractionQueryGenerator(unittest.TestCase):
    """
    Tests for DNA to Protein interactions

    Case: You have a DNA segment and want to find binding Protein

    """
    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.dna_segment1 = data_model.DnaSpecie(sequence = 'CCTTTGTT')
        self.dna_segment2 = data_model.DnaSpecie(sequence = 'AAGGTCAA')

    def test_filter_observed_values(self):
        q = dpi.DNAtoProteinInteractionQueryGenerator(cache_dirname=self.cache_dirname)

        observe = q.get_observed_result(self.dna_segment2)
        self.assertEqual(set(c.specie.gene_name for c in observe), set(['NR4A2', 'TRP(MYB) class']))


    def test_get_protein_by_binding_matrix(self):
        q = dpi.DNAtoProteinInteractionQueryGenerator(cache_dirname=self.cache_dirname)

        query = q.get_protein_by_DNA_sequence(self.dna_segment1.sequence)

        self.assertEqual(set(c[0].subunit_name for c in query), set(['pan', 'Sox2', 'DOF5.6']))
        self.assertEqual(set(c[1] for c in query), set([0,0,-8,-8]))
