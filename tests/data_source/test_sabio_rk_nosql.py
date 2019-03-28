'''Tests of sqlite_to_json
'''

import unittest
import shutil
import tempfile
from datanator.data_source import sabio_rk_nosql
from datanator.util import warning_util
warning_util.disable_warnings()


class TestSabioRkNoSQL(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.file_directory = '../../datanator/data_source/cache/SabioRk/'
        cls.client = 'mongodb://localhost:27017/'
        cls.db = 'test'
        cls.src = sabio_rk_nosql.SabioRkNoSQL(
            cls.file_directory, cls.db, cls.client, verbose = True, output_directory=cls.file_directory )

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    def setUp(self):
        (self.file_names, self.file_dict) = self.src.load_json()
        self.collection = self.src.con_db()

    @unittest.skip("test_con_db")
    def test_con_db(self):
        self.assertNotEqual(self.collection, 'Server not available')

    @unittest.skip("test_load_json")
    def test_load_json(self):
        null = None
        self.assertTrue('compartment' in self.file_names)
        self.assertTrue('synonym' in self.file_names)
        self.assertEqual(self.file_dict['entry'][0], {
            "_id": 1,
            "_type": "compound",
            "id": 2562,
            "name": "Peptide",
            "created": "2018-11-13 15:34:45.639403",
            "modified": "2018-11-13 18:06:01.266605"
        })
        self.assertEqual(self.file_dict['kinetic_law'][1], {
            "_id": 11,
            "enzyme_id": 10,
            "enzyme_compartment_id": null,
            "enzyme_type": "",
            "tissue": null,
            "mechanism": null,
            "equation": null,
            "taxon": 1467,
            "taxon_wildtype": 0,
            "taxon_variant": "S156E of subtilisin DSAI (N76D/N87S/S103A/V104I)",
            "temperature": 25.0,
            "ph": 7.5,
            "media": "50 mM potassium phosphate, 4 % DMSO"
        })

    def test_make_doc(self):
        session = self.src.make_doc(self.file_names, self.file_dict)
