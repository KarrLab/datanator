""" Test of text search

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-12-18
:Copyright: 2017, Karr Lab
:License: MIT
"""


from datanator.core import common_schema, models
from datanator.api.query import text_search
import tempfile
import shutil
import unittest

@unittest.skip('Major refactoring in `backend` branch')
class TestTextSearchSession(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.sesh = text_search.TextSearchSession(db_cache_dirname=cls.cache_dirname)

        flaskdb = common_schema.CommonSchema(cache_dirname = cls.cache_dirname)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    def test_list_return_objects(self):
        """
        Tests ability to collect objects from full text search of database
        """
        listed, dict_search = self.sesh.return_search('2-Oxopentanoate')
        for c in models.Metabolite.query.search('2-Oxopentanoate').all():
            self.assertIn(c, listed)

        listed, dict_search = self.sesh.return_search('MCM complex')
        for c in set(models.ProteinComplex.query.search('MCM complex').all()):
            self.assertIn(c, listed)

    def test_dict_return_objects(self):
        """
        Tests ability to collect objects from full text search of database
        """
        list_db_models, dict_search = self.sesh.return_search('2-Oxopentanoate')

        self.assertGreater(len(dict_search['Metabolite']), 0)
        self.assertEqual(len(dict_search['ProteinComplex']), 0)
        self.assertEqual(len(dict_search['ProteinSubunit']), 0)
        self.assertGreater(len(dict_search['Reaction']), 0)

        rank, dict_search =  self.sesh.return_search('MCM complex')
        print(rank)
        print(dict_search)
        self.assertEqual(len(dict_search['Metabolite']), 0)
        self.assertGreater(len(dict_search['ProteinComplex']), 0)
        self.assertGreater(len(dict_search['ProteinSubunit']), 0)
        self.assertEqual(len(dict_search['Reaction']), 0)
