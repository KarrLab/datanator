""" Test of text metabolite manager

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2018-08-08
:Copyright: 2018, Karr Lab
:License: MIT
"""

from kinetic_datanator.api.lib.subunit.manager import subunit_manager
from kinetic_datanator.core import common_schema, models
import unittest
import tempfile
import shutil

class TestProteinSubunitManager(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.protein_P00323 = subunit_manager.data_source.session.query(models.ProteinSubunit).filter_by(uniprot_id = 'P00323').first()

    def test_get_abundance_by_uniprot(self):
        uniprot = self.protein_P00323.uniprot_id
        abundances = subunit_manager.get_abundance_by_uniprot(uniprot).filter(
            models.AbundanceData.pax_load.in_([1, 2])).all()
        self.assertEqual(set(c.abundance for c in abundances),
                         set([1003.0, 1336.0]))

    def test_get_abundance_by_gene_name(self):
        gene_name = 'rplO'
        abundances = subunit_manager.get_abundance_by_gene_name(gene_name).filter(
            models.AbundanceData.pax_load.in_([1, 2])).all()
        self.assertEqual(set(c.abundance for c in abundances),
                         set([1027.0, 1861.0]))

    def test_get_abundance_by_sequence(self):

        sequence = self.protein_P00323.canonical_sequence
        abundances = subunit_manager.get_abundance_by_sequence(sequence).filter(
            models.AbundanceData.pax_load.in_([1, 2])).all()
        self.assertEqual(set(c.abundance for c in abundances),
                         set([1003.0, 1336.0]))

    @unittest.skip('entrez needs to be fixed')
    def test_get_abundance_by_entrez(self):
        entrez_id = self.protein_P00323.entrez_id
        abundances = subunit_manager.get_abundance_by_entrez(entrez_id).filter(
            models.AbundanceData.pax_load.in_([1, 2])).all()
        self.assertEqual(set(c.abundance for c in abundances),
                         set([1003.0, 1336.0]))

    def test_get_abundance_by_mass(self):
        mass = self.protein_P00323.mass
        abundances = subunit_manager.get_abundance_by_mass(mass).filter(
            models.AbundanceData.pax_load.in_([1, 2])).all()
        self.assertEqual(set(c.abundance for c in abundances),
                         set([1003.0, 1336.0]))

    def test_get_abundance_by_length(self):
        length = self.protein_P00323.length
        abundances = subunit_manager.get_abundance_by_length(length).filter(
            models.AbundanceData.pax_load.in_([1, 2])).all()
        self.assertEqual(set(c.abundance for c in abundances),
                         set([1003.0, 1336.0, 1027.0, 1861.0]))
