import unittest
from datanator.flask import front_end_query
import tempfile
import shutil

class TestQueryFrontEnd(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.src = front_end_query.QueryFrontEnd()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    def string_query(self):
        string = '2-Ketobutyric acid (alpha ketobutyric acid)'
        result = self.src.string_query(string)
        self.assertEqual(result[0]['m2m_id'], 'M2MDB000001')