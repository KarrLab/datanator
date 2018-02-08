""" Test of text search

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-12-18
:Copyright: 2017, Karr Lab
:License: MIT
"""


from kinetic_datanator.core import flask_common_schema, models
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

        flaskdb = flask_common_schema.FlaskCommonSchema(cache_dirname = self.cache_dirname)

        for item in flaskdb.text_indicies:
            flask_whooshalchemy.whoosh_index(models.app, item)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)



    def test_return_objects(self):
        """
        Tests ability to collect objects from full text search of database
        """
        search= self.sesh.return_search('2-Oxopentanoate')
        for c in search:
            self.assertIn(c, models.Compound.query.whoosh_search('2-Oxopentanoate').all())

        search= self.sesh.return_search('MCM complex')
        for c in search:
            self.assertIn(c, set(models.ProteinComplex.query.whoosh_search('MCM complex').all()))
