""" Test of text search

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-12-18
:Copyright: 2017, Karr Lab
:License: MIT
"""


from kinetic_datanator.app import text_search, flask_common_schema, models
import tempfile
import shutil
import flask_whooshalchemy
import unittest

class TestTextSearchSession(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.sesh = text_search.TextSearchSession(db_cache_dirname=self.cache_dirname)

        flaskdb = flask_common_schema.FlaskCommonSchema()

        for item in flaskdb.text_indicies:
            flask_whooshalchemy.whoosh_index(models.app, item)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)



    def test_return_objects(self):
        """
        Tests ability to collect objects from full text search of database
        """
        search_dict= self.sesh.return_search('2-Oxopentanoate')
        for c in search_dict['Compound']:
            self.assertIn(c, models.Compound.query.whoosh_search('2-Oxopentanoate').all())

        search_dict= self.sesh.return_search('MCM complex')
        self.assertEqual(set([c for c in search_dict['ProteinComplex']]),
                         set(models.ProteinComplex.query.whoosh_search('MCM complex').all()))

        search_dict = self.sesh.return_search('P49418')
        for c in search_dict['ProteinInteractions']:
            self.assertIn(c, set(models.ProteinInteractions.query.whoosh_search('P49418').all()))
