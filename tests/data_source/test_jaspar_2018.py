"""
This module tests all aspects of jaspar.py

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Author: Jonathan Karr <jonrkar@gmail.com>
:Date: 2017-07-24
:Copyright: 2017, Karr Lab
:License: MIT

"""

from kinetic_datanator.data_source import jaspar2018
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import random
import requests
import shutil
import tempfile
import unittest


class TestLoad(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.jaspar = jaspar2018.JASPAR2018(cache_dirname = self.cache_dirname)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)

    def test_load(self):

        self.assertTrue(self.cache_dirname+'/JASPAR2018.sqlite')
        q = self.jaspar.session.query(jaspar2018.Matrix).all()
        self.assertTrue(q)

    def test_entries(self):

        q = self.jaspar.session.query(jaspar2018.Matrix).get(9237)

        self.assertEqual(q.COLLECTION,'CORE')
        self.assertEqual(q.BASE_ID, 'MA0009')
        self.assertEqual(q.NAME,'T')

        q = self.jaspar.session.query(jaspar2018.Annotation).filter_by(ID = 9249).all()

        self.assertEqual(set(c.VAL for c in q), \
        set(['C2H2 zinc finger factors', '-', 'Dof-type', '10074718', 'plants', 'SELEX', '28']))

        q = self.jaspar.session.query(jaspar2018.Species).filter_by(ID = 9238).first()

        self.assertEqual(q.TAX_ID, '7227')
