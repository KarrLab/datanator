# -*- coding: utf-8 -*-

""" Test of common schema

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-07-31
:Copyright: 2017, Karr Lab
:License: MIT
"""
import unittest
from sqlalchemy.orm import sessionmaker
from kinetic_datanator.core import common_schema
import tempfile
import shutil


class ShortTestCommonSchema(unittest.TestCase):
    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_working(self):
        cs = common_schema.CommonSchema(cache_dirname = self.cache_dirname , clear_content = True,
            load_content=False, download_backup=False, max_entries = 5)
        # cs = common_schema.CommonSchema(name = 'aggregate', clear_content = True,
        #                                 load_content=False, download_backup=False,
        #                                 max_entries = 5)
        cs.load_content()
        session = cs.session

        #TODO: Build More Tests

        subunit = session.query(common_schema.ProteinSubunit).filter_by(uniprot_id = 'P41182').first()
        self.assertEqual(subunit.entrez_id, 604)

        taxon = session.query(common_schema.Taxon).get(882)
        self.assertEqual(taxon.name, 'D.vulgaris')


# class LongTestCommonSchema(unittest.TestCase):
#     def setUp(self):
#         self.cache_dirname = tempfile.mkdtemp()
#
#     def tearDown(self):
#         shutil.rmtree(self.cache_dirname)
#
#     def test_working(self):
#         # cs = common_schema.CommonSchema(cache_dirname = self.cache_dirname , clear_content = True, load_content=False, download_backup=False)
#         cs = common_schema.CommonSchema(name = 'aggregate', clear_content = True,
#                                         load_content=False, download_backup=False
#                                         max_entries = 20)
#         cs.load_content()
#         session = cs.session
#
#
#         subunit = session.query(common_schema.Subunit).filter_by(uniprot_id = 'P41182').first()
#         self.assertEqual(subunit.entrez_id, 604)
#
#         taxon = session.query(common_schema.Taxon).get(882)
#         self.assertEqual(taxon.species_name, 'D.vulgaris')
