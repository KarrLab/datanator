""" Tests of DNA and Protein interaction queries

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-09-27
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model, common_schema
from kinetic_datanator.data_query import dna_protein_interactions as dpi
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
    def setUp(self):
        self.arnt  = data_model.ProteinSpecie(uniprot_id = 'P53762', gene_name = 'ARNT')

    def test_filter_observed_values(self):
        q = dpi.ProteintoDNAInteractionQueryGenerator()

        observable = q.get_observed_values(self.arnt)

        self.assertEqual(observable[0].specie.sequence, 'CACGTG')
        self.assertEqual(observable[0].specie.cross_references[0].id, '7592839')


    def test_get_DNA_by_protein(self):
        q = dpi.ProteintoDNAInteractionQueryGenerator()

        position = q.get_DNA_by_protein(self.arnt)

        self.assertEqual(set(c.frequency_a for c in position[0]), set([0,19,4]))
        self.assertEqual(set(c.position for c in position[0]), set([1, 2, 3, 4, 5, 6]))


# class TestDNAtoProteinInteractionQueryGenerator(unittest.TestCase):
#     """
#     Tests for DNA to Protein interactions
#
#     Case: You have a DNA segment and want to find binding Protein
#
#     """
#     @classmethod
#     def setUp(self):
#         q = dpi.ProteintoDNAInteractionQueryGenerator()
#         self.srf = data_model.ProteinSpecie(uniprot_id = 'P11831', sequence = 'SRF')
#         observable = q.get_observed_values(self.srf)
#         self.dna_segment = data_model.DnaSpecie(binding_matrix = observable[0].specie.binding_matrix,
#                                                 sequence = observable[0].specie.sequence)
#
#     def test_get_protein_by_binding_matrix(self):
#         """
#
#         """
#
#         q = dpi.DNAtoProteinInteractionQueryGenerator()
#
#         query = q.get_protein_by_binding_matrix(self.dna_segment.binding_matrix)
#
#         self.assert
