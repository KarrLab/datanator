# -*- coding: utf-8 -*-

""" Test of pax database

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-08-1
:Copyright: 2017, Karr Lab
:License: MIT
"""
import unittest
from kinetic_datanator.data_source import pax
from sqlalchemy.orm import sessionmaker
import tempfile
import shutil




class TestQuery(unittest.TestCase):
        def setUp(self):
            self.cache_dirname = tempfile.mkdtemp()

        def tearDown(self):
            shutil.rmtree(self.cache_dirname)

        def test_query(self):
            src = pax.Pax(cache_dirname = self.cache_dirname, clear_content = False, load_content=False, download_backup=False, verbose = True)
            ## Fraction of DB to load
            src.fraction = .005
            src.load_content()
            session = src.session


            z = session.query(pax.Dataset).get(2)
            self.assertIsInstance(z.score, float)
            self.assertIsInstance(z.taxon_ncbi_id, int)
            self.assertEqual(z.score, 0.03)

            y = session.query(pax.Observation).get(5)
            self.assertIsInstance(y.abundance, float)
            self.assertEqual(y.dataset_id, 1)

            x = session.query(pax.Taxon).get(882)
            self.assertEqual(x.species_name, 'D.vulgaris')
