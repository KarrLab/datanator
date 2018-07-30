""" Test of text search

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-12-18
:Copyright: 2017, Karr Lab
:License: MIT
"""


from kinetic_datanator.core import common_schema, models
from kinetic_datanator.data_query import text_search
import tempfile
import shutil
import flask_whooshalchemy
import unittest

class TestTextSearchSession(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.sesh = text_search.TextSearchSession(db_cache_dirname=self.cache_dirname)

        flaskdb = common_schema.CommonSchema(cache_dirname = self.cache_dirname)

        for item in flaskdb.text_indicies:
            flask_whooshalchemy.whoosh_index(models.app, item)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)



    def test_list_return_objects(self):
        """
        Tests ability to collect objects from full text search of database
        """
        listed, dict_search = self.sesh.return_search('2-Oxopentanoate')
        for c in models.Compound.query.whoosh_search('2-Oxopentanoate').all():
            self.assertIn(c, listed)

        listed, dict_search = self.sesh.return_search('MCM complex')
        for c in set(models.ProteinComplex.query.whoosh_search('MCM complex').all()):
            self.assertIn(c, listed)

    def test_dict_return_objects(self):
        """
        Tests ability to collect objects from full text search of database
        """
        list_db_models, dict_search = self.sesh.return_search('2-Oxopentanoate')
        self.assertGreater(len(dict_search['Compound']), 0)
        self.assertEqual(len(dict_search['ProteinComplex']), 0)
        self.assertEqual(len(dict_search['ProteinSubunit']), 0)
        self.assertGreater(len(dict_search['Reaction']), 0)

        rank, dict_search =  self.sesh.return_search('MCM complex')
        self.assertEqual(len(dict_search['Compound']), 0)
        self.assertGreater(len(dict_search['ProteinComplex']), 0)
        self.assertEqual(len(dict_search['ProteinSubunit']), 0)
        self.assertEqual(len(dict_search['Reaction']), 0)
