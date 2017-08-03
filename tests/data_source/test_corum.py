# -*- coding: utf-8 -*-

""" Test of corum database

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-08-1
:Copyright: 2017, Karr Lab
:License: MIT
"""

import unittest
from kinetic_datanator.data_source import corum
from sqlalchemy.orm import sessionmaker
import tempfile
import shutil

class TestCorumDBCreation(unittest.TestCase):
        def setUp(self):
            self.cache_dirname = tempfile.mkdtemp()

        def tearDown(self):
            shutil.rmtree(self.cache_dirname)

        def test_load_some_content(self):
            src = corum.Corum(cache_dirname = self.cache_dirname, clear_content = True, verbose = False, max_entries = 10)
            src.load_content()
            session = src.session

            subunit = session.query(corum.Subunit).get(3)
            self.assertEqual(subunit.su_uniprot, 'P41182')
            self.assertEqual(str(subunit.protein_name), 'B-cell lymphoma 6 protein ')

            complx = subunit.complex

            self.assertEqual(complx.complex_id, 2)
            self.assertEqual(complx.complex_name, 'BCL6-HDAC5 complex')

            obs = complx.observation
            self.assertEqual(obs.id,2)
            self.assertEqual(obs.pubmed_id, 11929873)

            tax = obs.taxon

            self.assertEqual(tax.swissprot_id, 'Homo sapiens (Human)')


        def test_load_all_content(self):
            src = corum.Corum(cache_dirname = self.cache_dirname, clear_content = True, verbose = False)
            src.load_content()
            session = src.session

            c = session.query(corum.Complex).get(80)
            self.assertEqual(c.complex_name, 'Ubiquitin E3 ligase (SKP1A, SKP2, CUL1, RBX1)')

            s = session.query(corum.Subunit).filter(corum.Subunit.su_uniprot == 'Q9UQL6').first()
            self.assertEqual(s.protein_name, 'Histone deacetylase 5')

            t = session.query(corum.Taxon).filter(corum.Taxon.ncbi_id == 9606).first()
            self.assertEqual(t.swissprot_id, 'Homo sapiens (Human)')
