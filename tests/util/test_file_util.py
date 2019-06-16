import unittest
from datanator.util import file_util
import tempfile
import shutil


class TestFileUtil(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.src = file_util.FileUtil()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)


    def test_flatten_json(self):
        dic = {
            "a": 1,
            "b": 2,
            "c": [{"d": [2, 3, 4], "e": [{"f": 1, "g": 2}]}],
            'h': [1, 2, 3]
        }
        flat = self.src.flatten_json(dic)
        self.assertEqual(flat['c_0_d_0'], 2)
        self.assertEqual(flat['c_0_e_0_f'], 1)
        self.assertEqual(flat['h_0'], 1)

    def test_extract_values(self):
        dic = {
            "a": 1,
            "b": 2,
            "c": [{"d": [2, 3, 4], "e": [{"f": 1, "g": 2}]}],
            'h': [1, 2, 3]
        }
        values = self.src.extract_values(dic, 'g')
        self.assertEqual([2], values)