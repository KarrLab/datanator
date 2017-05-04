""" Tests of kinetic_datanator

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-06
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator import datanator
from kinetic_datanator.core import data_structs
from os import path
import os
import sys
import unittest


class TestDatanator(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.fixtures_dir = path.join(path.dirname(__file__), "fixtures")

        out_dir = cls.output_dir = path.join(path.dirname(__file__), "output")
        if not path.isdir(out_dir):
            os.makedirs(out_dir)

    def test_init(self):
        d = datanator.Datanator(max_taxon_dist=10, include_mutants=True,
                                min_temp=34, max_temp=40, min_ph=6, max_ph=8)
        self.assertEqual(d.max_taxon_dist, 10)
        self.assertEqual(d.include_mutants, True)
        self.assertEqual(d.min_temp, 34)
        self.assertEqual(d.max_temp, 40)
        self.assertEqual(d.min_ph, 6)
        self.assertEqual(d.max_ph, 8)

    def test_read_model(self):
        i_file = path.join(self.fixtures_dir, "ump_kinase.xlsx")
        d = datanator.Datanator()
        model = d.read_model(i_file)
        taxon, compartments, molecules, reactions = model
        self.assertEqual(taxon.name, 'Mycoplasma pneumoniae M129')

    def test_annotate_model(self):
        i_file = path.join(self.fixtures_dir, "ump_kinase.xlsx")
        d = datanator.Datanator()
        model = d.read_model(i_file)
        d.annotate_model(model)

        taxon, compartments, molecules, reactions = model

        mol = next(mol for mol in molecules if mol.id == 'UMP')
        ids = set(xr.id for xr in mol.cross_references if xr.source == 'sabio-id')
        names = set(xr.id for xr in mol.cross_references if xr.source == 'sabio-name')
        self.assertEqual(ids, set([7458, 1300, 25061]))
        self.assertEqual(names, set(["2,4-Dioxotetrahydropyrimidine D-ribonucleotide", "UMP", "Uridine 5'-phosphate"]))

        rxn = next(rxn for rxn in reactions if rxn.id == 'ump_kinase')
        self.assertEqual(rxn.get_ec_number(), '2.7.4')  # true EC is 2.7.4.22

    def test_annotate_molecules(self):
        i_file = path.join(self.fixtures_dir, "ump_kinase.xlsx")
        d = datanator.Datanator()
        model = d.read_model(i_file)
        d.annotate_molecules(model)

        taxon, compartments, molecules, reactions = model

        mol = next(mol for mol in molecules if mol.id == 'UMP')
        ids = set(xr.id for xr in mol.cross_references if xr.source == 'sabio-id')
        names = set(xr.id for xr in mol.cross_references if xr.source == 'sabio-name')
        self.assertEqual(ids, set([7458, 1300, 25061]))
        self.assertEqual(names, set(["2,4-Dioxotetrahydropyrimidine D-ribonucleotide", "UMP", "Uridine 5'-phosphate"]))

    def test_annotate_reactions(self):
        i_file = path.join(self.fixtures_dir, "ump_kinase.xlsx")
        d = datanator.Datanator()
        model = d.read_model(i_file)
        d.annotate_reactions(model)

        taxon, compartments, molecules, reactions = model
        rxn = next(rxn for rxn in reactions if rxn.id == 'ump_kinase')
        self.assertEqual(rxn.get_ec_number(), '2.7.4')  # true EC is 2.7.4.22

    @unittest.skip('me')
    def test_datanator(self):
        # Test getting kinetic data about a reaction that is relevant for a species

        i_file = path.join(self.fixtures_dir, "ump_kinase.xlsx")
        o_file = path.join(self.output_dir, "ump_kinase.xlsx")

        # get kinetic data for reations in Excel sheet
        rxns = datanator.get_kinetic_data(i_file, o_file)

        # test correct number of output reactions
        self.assertEqual(len(rxns), 1)

        # test outputted reaction is UMP kinases
        rxn = rxns[0]
        self.assertEqual(rxn.id, 'ump_kinase')

        # test the formatted data fields
        self.assertEqual(rxn.reaction_ids, ['201'])

        # median_km_entry is an object of the Entry class found in the sabio_rk module
        # the following are a series of tests that make sure that the entry fields are
        #filled in properly
        median_entry = rxn.km_data.median_entry
        self.assertNotEqual(median_entry, None)
        self.assertEqual(median_entry.entry_id, '17927')
        self.assertEqual(median_entry.vmax, '0.00665')
        self.assertEqual(median_entry.proximity, 6)
