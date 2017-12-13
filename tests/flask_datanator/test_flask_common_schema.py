# -*- coding: utf-8 -*-

""" Test of common schema

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-07-31
:Copyright: 2017, Karr Lab
:License: MIT
"""
import unittest
from kinetic_datanator.flask_datanator import flask_common_schema, models
import flask_whooshalchemy
import tempfile
import shutil
import random



class LoadingTestCommonSchema(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.cs = flask_common_schema.FlaskCommonSchema(cache_dirname= self.cache_dirname,
                                clear_content = True, load_entire_small_DBs = False,
                                download_backup= False, load_content = True, max_entries = 10,
                                verbose = True)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)


    def test_whoosh(self):
        self.assertEqual(set([c.name for c in models.Compound.query.whoosh_search('adenine').all()]),
            set(['Adenosine', 'Adenosine monophosphate', 'Cyclic AMP', "Adenosine 3',5'-diphosphate", 'Adenine']))
