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
            src = pax.Pax(cache_dirname = self.cache_dirname , clear_content = False, load_content=False, download_backup=False, verbose = True)
            ## Fraction of DB to load
            src.fraction = .005
            src.load_content()
            session = src.session


            z = session.query(pax.Dataset).filter(pax.Dataset.score == 9.76).first()
            self.assertEqual(z.file_name, '10090/10090-Adrenal_gland_geiger_2013.txt')

            y = session.query(pax.Observation).get(658)
            self.assertEqual(round(y.abundance,2), 0.34)
            self.assertEqual(y.dataset_id, 1)
            self.assertEqual(y.protein_id, 2100299)

            x = session.query(pax.Taxon).get(10090)
            self.assertEqual(x.species_name, 'M.musculus')
