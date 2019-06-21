import unittest
from datanator.flask.query import front_end_query
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

    def test_string_query(self):
        string = '2-Ketobutyric acid (alpha ketobutyric acid)'
        results = self.src.string_query(string)
        self.assertEqual(results[0]['m2m_id'], 'M2MDB001633')