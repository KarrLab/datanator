'''Tests of sqlite_to_json

:Author: Zhouyang Lian <zhouyang.lian@familian.life>
:Author: Jonathan <jonrkarr@gmail.com>
:Date: 2019-04-02
:Copyright: 2019, Karr Lab
:License: MIT
'''

import unittest
import shutil
import tempfile
from datanator.data_source import sabio_rk_nosql

class TestSabioRkNoSQL(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'test'
        cls.MongoDB = 'mongodb://mongo:27017/'
        cls.quilt_package = 'sabiork_nosql'
        cls.src = sabio_rk_nosql.SabioRkNoSQL(
            db = cls.db, MongoDB = cls.MongoDB, cache_directory = cls.cache_dirname,
            quilt_package = cls.quilt_package, verbose = True, max_entries = 20)
        (cls.file_names, cls.file_dict) = cls.src.load_json()
        cls.collection = cls.src.con_db('sabio_rk')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        # cls.collection.drop()

    @unittest.skip("passed")
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

    # @unittest.skip("test_make_doc")
    # def test_make_doc(self):
    #     session = self.src.make_doc(self.file_names, self.file_dict)

    def test_add_deprot_inchi(self):
        self.src.add_deprot_inchi()
        cursor = self.collection.find({'kinlaw_id': 2})
        doc = cursor[0]
        self.assertEqual(doc['reaction_participant'][0]['substrate'][0]['inchi_deprot'], "InChI=1S/H2O")
