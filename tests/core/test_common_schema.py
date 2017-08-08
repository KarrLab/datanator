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


class TestCommonSchema(unittest.TestCase):
    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_working(self):
        cs = common_schema.CommonSchema(cache_dirname = self.cache_dirname , clear_content = True, load_content=False, download_backup=False)
        # cs = CommonSchema(name = 'aggregate', clear_content = True, load_content=False, download_backup=False)
        cs.load_content()
        session = cs.session

        taxon = session.query(common_schema.Taxon).get(882)
        self.assertEqual(taxon.species_name, 'D.vulgaris')