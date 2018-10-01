""" Test of taxonomy utilities

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-10
:Copyright: 2017, Karr Lab
:License: MIT
"""

from datanator.util import taxonomy_util
import unittest


class TestTaxonomyUtil(unittest.TestCase):

    def test_setup_database(self):
        taxonomy_util.setup_database()

    def test_Taxon_init_from_name(self):
        taxon = taxonomy_util.Taxon(name='mycoplasma')
        self.assertEqual(taxon.name, 'Mycoplasma')
        self.assertEqual(taxon.id_of_nearest_ncbi_taxon, 2093)
        self.assertEqual(taxon.distance_from_nearest_ncbi_taxon, 0)
        self.assertEqual(taxon.additional_name_beyond_nearest_ncbi_taxon, '')

        taxon = taxonomy_util.Taxon(name='mycoplasma genitalium')
        self.assertEqual(taxon.name, 'Mycoplasma genitalium')
        self.assertEqual(taxon.id_of_nearest_ncbi_taxon, 2097)
        self.assertEqual(taxon.distance_from_nearest_ncbi_taxon, 0)
        self.assertEqual(taxon.additional_name_beyond_nearest_ncbi_taxon, '')

        taxon = taxonomy_util.Taxon(name='mycoplasma genitalium G37')
        self.assertEqual(taxon.name, 'Mycoplasma genitalium G37')
        self.assertEqual(taxon.id_of_nearest_ncbi_taxon, 243273)
        self.assertEqual(taxon.distance_from_nearest_ncbi_taxon, 0)
        self.assertEqual(taxon.additional_name_beyond_nearest_ncbi_taxon, '')

        taxon = taxonomy_util.Taxon(name='mycoplasma genitalium XXX')
        self.assertEqual(taxon.name, 'Mycoplasma genitalium XXX')
        self.assertEqual(taxon.id_of_nearest_ncbi_taxon, 2097)
        self.assertEqual(taxon.distance_from_nearest_ncbi_taxon, 1)
        self.assertEqual(taxon.additional_name_beyond_nearest_ncbi_taxon, ' XXX')

        taxon = taxonomy_util.Taxon(name='mycoplasma XXX')
        self.assertEqual(taxon.name, 'Mycoplasma XXX')
        self.assertEqual(taxon.id_of_nearest_ncbi_taxon, 2093)
        self.assertEqual(taxon.distance_from_nearest_ncbi_taxon, 1)
        self.assertEqual(taxon.additional_name_beyond_nearest_ncbi_taxon, ' XXX')

        taxon = taxonomy_util.Taxon(name='mycoplasma XXX YYY')
        self.assertEqual(taxon.name, 'Mycoplasma XXX YYY')
        self.assertEqual(taxon.id_of_nearest_ncbi_taxon, 2093)
        self.assertEqual(taxon.distance_from_nearest_ncbi_taxon, 2)
        self.assertEqual(taxon.additional_name_beyond_nearest_ncbi_taxon, ' XXX YYY')

    def test_Taxon_init_from_id(self):
        taxon = taxonomy_util.Taxon(ncbi_id=2093)
        self.assertEqual(taxon.name, 'Mycoplasma')
        self.assertEqual(taxon.id_of_nearest_ncbi_taxon, 2093)
        self.assertEqual(taxon.distance_from_nearest_ncbi_taxon, 0)
        self.assertEqual(taxon.additional_name_beyond_nearest_ncbi_taxon, '')

        taxon = taxonomy_util.Taxon(ncbi_id=2097)
        self.assertEqual(taxon.name, 'Mycoplasma genitalium')
        self.assertEqual(taxon.id_of_nearest_ncbi_taxon, 2097)
        self.assertEqual(taxon.distance_from_nearest_ncbi_taxon, 0)
        self.assertEqual(taxon.additional_name_beyond_nearest_ncbi_taxon, '')

        taxon = taxonomy_util.Taxon(ncbi_id=243273)
        self.assertEqual(taxon.name, 'Mycoplasma genitalium G37')
        self.assertEqual(taxon.id_of_nearest_ncbi_taxon, 243273)
        self.assertEqual(taxon.distance_from_nearest_ncbi_taxon, 0)
        self.assertEqual(taxon.additional_name_beyond_nearest_ncbi_taxon, '')

        self.assertRaises(ValueError, taxonomy_util.Taxon, ncbi_id=-1)

    def test_Taxon_get_ncbi_id(self):
        taxon = taxonomy_util.Taxon(ncbi_id=2093)
        self.assertEqual(taxon.get_ncbi_id(), 2093)

        taxon = taxonomy_util.Taxon(name='mycoplasma XXX')
        self.assertEqual(taxon.get_ncbi_id(), None)

    def test_Taxon_get_parent_taxa(self):
        parents = taxonomy_util.Taxon(name='mycoplasma genitalium').get_parent_taxa()
        self.assertEqual(parents[-1].name, 'Mycoplasma')
        self.assertEqual(parents[-2].name, 'Mycoplasmataceae')

        parents = taxonomy_util.Taxon(name='mycoplasma genitalium G37').get_parent_taxa()
        self.assertEqual(parents[-1].name, 'Mycoplasma genitalium')
        self.assertEqual(parents[-2].name, 'Mycoplasma')
        self.assertEqual(parents[-3].name, 'Mycoplasmataceae')

        parents = taxonomy_util.Taxon(name='mycoplasma genitalium XXX').get_parent_taxa()
        self.assertEqual(parents[-1].name, 'Mycoplasma genitalium')
        self.assertEqual(parents[-2].name, 'Mycoplasma')
        self.assertEqual(parents[-3].name, 'Mycoplasmataceae')

    def test_Taxon_get_rank(self):
        self.assertEqual(taxonomy_util.Taxon(ncbi_id=2093).get_rank(), 'genus')
        self.assertEqual(taxonomy_util.Taxon(ncbi_id=2097).get_rank(), 'species')
        self.assertEqual(taxonomy_util.Taxon(ncbi_id=243273).get_rank(), None)

    def test_Taxon_get_common_ancestor(self):
        taxon_a = taxonomy_util.Taxon(name='mycoplasma genitalium')
        taxon_b = taxonomy_util.Taxon(name='mycoplasma pneumoniae')
        self.assertEqual(taxon_a.get_common_ancestor(taxon_b).name, 'Mycoplasma')

        taxon_a = taxonomy_util.Taxon(name='mycoplasma genitalium G37')
        taxon_b = taxonomy_util.Taxon(name='mycoplasma genitalium XXX')
        self.assertEqual(taxon_a.get_common_ancestor(taxon_b).name, 'Mycoplasma genitalium')

        taxon_a = taxonomy_util.Taxon(name='mycoplasma genitalium G37')
        taxon_b = taxonomy_util.Taxon(name='mycoplasma pneumoniae')
        self.assertEqual(taxon_a.get_common_ancestor(taxon_b).name, 'Mycoplasma')

        taxon_a = taxonomy_util.Taxon(name='mycoplasma genitalium G37')
        taxon_b = taxonomy_util.Taxon(name='mycoplasma')
        self.assertEqual(taxon_a.get_common_ancestor(taxon_b).name, 'Mycoplasma')

        taxon_a = taxonomy_util.Taxon(name='mycoplasma genitalium XXX')
        taxon_b = taxonomy_util.Taxon(name='mycoplasma')
        self.assertEqual(taxon_a.get_common_ancestor(taxon_b).name, 'Mycoplasma')

        taxon_a = taxonomy_util.Taxon(name='mycoplasma genitalium G37')
        taxon_b = taxonomy_util.Taxon(name='escherichia coli')
        self.assertEqual(taxon_a.get_common_ancestor(taxon_b).name, 'Bacteria')

    def test_get_distance_to_common_ancestor(self):
        taxon_a = taxonomy_util.Taxon(name='mycoplasma genitalium')
        taxon_b = taxonomy_util.Taxon(name='mycoplasma pneumoniae')
        self.assertEqual(taxon_a.get_distance_to_common_ancestor(taxon_b), 1)

        taxon_a = taxonomy_util.Taxon(name='mycoplasma genitalium G37')
        taxon_b = taxonomy_util.Taxon(name='mycoplasma genitalium XXX')
        self.assertEqual(taxon_a.get_distance_to_common_ancestor(taxon_b), 1)

        taxon_a = taxonomy_util.Taxon(name='mycoplasma genitalium G37')
        taxon_b = taxonomy_util.Taxon(name='mycoplasma pneumoniae')
        self.assertEqual(taxon_a.get_distance_to_common_ancestor(taxon_b), 2)

        taxon_a = taxonomy_util.Taxon(name='mycoplasma genitalium G37')
        taxon_b = taxonomy_util.Taxon(name='mycoplasma')
        self.assertEqual(taxon_a.get_distance_to_common_ancestor(taxon_b), 2)

        taxon_a = taxonomy_util.Taxon(name='mycoplasma genitalium XXX')
        taxon_b = taxonomy_util.Taxon(name='mycoplasma')
        self.assertEqual(taxon_a.get_distance_to_common_ancestor(taxon_b), 2)

        taxon_a = taxonomy_util.Taxon(name='mycoplasma genitalium G37')
        taxon_b = taxonomy_util.Taxon(name='escherichia coli')
        self.assertEqual(taxon_a.get_distance_to_common_ancestor(taxon_b), 8)

        taxon_a = taxonomy_util.Taxon(name='mycoplasma genitalium XXX')
        taxon_b = taxonomy_util.Taxon(name='escherichia coli')
        self.assertEqual(taxon_a.get_distance_to_common_ancestor(taxon_b), 8)

        taxon_a = taxonomy_util.Taxon(name='mycoplasma genitalium')
        taxon_b = taxonomy_util.Taxon(name='escherichia coli')
        self.assertEqual(taxon_a.get_distance_to_common_ancestor(taxon_b), 7)

    def test_get_distance_to_root(self):
        self.assertEqual(taxonomy_util.Taxon(name='mycoplasma').get_distance_to_root(), 8)
        self.assertEqual(taxonomy_util.Taxon(name='mycoplasma genitalium').get_distance_to_root(), 9)
        self.assertEqual(taxonomy_util.Taxon(name='mycoplasma genitalium G37').get_distance_to_root(), 10)
        self.assertEqual(taxonomy_util.Taxon(name='mycoplasma genitalium XXX').get_distance_to_root(), 10)
        self.assertEqual(taxonomy_util.Taxon(name='mycoplasma genitalium XXX YYY').get_distance_to_root(), 11)

    def test_get_max_distance_to_common_ancestor(self):
        self.assertEqual(taxonomy_util.Taxon(name='mycoplasma').get_max_distance_to_common_ancestor(), 8)
        self.assertEqual(taxonomy_util.Taxon(name='mycoplasma genitalium').get_max_distance_to_common_ancestor(), 9)
        self.assertEqual(taxonomy_util.Taxon(name='mycoplasma genitalium G37').get_max_distance_to_common_ancestor(), 10)
        self.assertEqual(taxonomy_util.Taxon(name='mycoplasma genitalium XXX').get_max_distance_to_common_ancestor(), 10)
        self.assertEqual(taxonomy_util.Taxon(name='mycoplasma genitalium XXX YYY').get_max_distance_to_common_ancestor(), 11)
