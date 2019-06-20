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
import datanator.config.core

class TestSabioRkNoSQL(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'test'
        username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
        port = datanator.config.core.get_config()['datanator']['mongodb']['port']
        replSet = datanator.config.core.get_config()['datanator']['mongodb']['replSet']
        cls.quilt_package = 'sabiork_nosql'
        cls.src = sabio_rk_nosql.SabioRkNoSQL(
            db = cls.db, MongoDB =MongoDB, cache_directory = cls.cache_dirname,
            quilt_package = cls.quilt_package, verbose = True, username = username,
            password = password, max_entries = 10)
        cls.client, cls.db_obj, cls.collection = cls.src.con_db('sabio_rk')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.client.close()

    @unittest.skip("passed")
    def test_load_json(self):
        null = None
        file_names, file_dict = self.src.load_json()
        self.assertTrue('compartment' in file_names)
        self.assertTrue('synonym' in file_names)
        self.assertEqual(file_dict['entry'][0], {
            "_id": 1,
            "_type": "compound",
            "id": 2562,
            "name": "Peptide",
            "created": "2018-11-13 15:34:45.639403",
            "modified": "2018-11-13 18:06:01.266605"
        })
        self.assertEqual(file_dict['kinetic_law'][1], {
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
    #     self.src.make_doc(self.file_names, self.file_dict)

    @unittest.skip("test_make_doc")
    def test_add_deprot_inchi(self):
        self.src.add_deprot_inchi()
        cursor = self.collection.find_one({'kinlaw_id': 2})
        self.assertEqual(cursor['reaction_participant'][0]['substrate'][0]['inchi_deprot'], "InChI=1S/H2O")
        self.assertEqual(cursor['reaction_participant'][0]['substrate'][0]['hashed_inchi'], 
            'ae41b33d263ff93c034384992c84ed94a4d73cb0d200e1ebfed3871c')