'''Tests of metabolite_nosql

:Author: Zhouyang Lian <zhouyang.lian@familian.life>
:Author: Jonathan <jonrkarr@gmail.com>
:Date: 2019-04-02
:Copyright: 2019, Karr Lab
:License: MIT
'''

import unittest
import shutil
import tempfile
from datanator.data_source import metabolite_nosql
import datanator.config.core
import os
import json


class TestMetaboliteNoSQL(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.source = 'ecmdb' # 'ymdb' or 'ecmdb'
        cls.db = 'datanator'
        username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
        port = datanator.config.core.get_config()['datanator']['mongodb']['port']
        replSet = datanator.config.core.get_config()['datanator']['mongodb']['replSet']
        cls.output_directory = cls.cache_dirname # directory to store JSON files
        cls.src = metabolite_nosql.MetaboliteNoSQL(cls.output_directory,
            cls.source, MongoDB, cls.db, verbose = True, max_entries=20,
            username = username, password = password)
        cls.client, cls.db_obj, cls.collection = cls.src.con_db(cls.source)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.client.close()

    # @unittest.skip('no json files in circle')
    def test_write_to_json(self):
        session = self.src.write_to_json()
        null = None
        if self.source == 'ymdb':
            ymdb_6 = self.collection.find({"ymdb_id": "YMDB00006"})[0]
            self.assertEqual(ymdb_6['ymdb_id'], "YMDB00006")
            self.assertEqual(ymdb_6['species'], "Saccharomyces cerevisiae")
            self.assertEqual(ymdb_6['name'], "1D-Myo-inositol 1,4,5,6-tetrakisphosphate")

            ymdb_10 = self.collection.find({"ymdb_id": "YMDB00010"})[0]
            self.assertEqual(ymdb_10['ymdb_id'], "YMDB00010")
            self.assertEqual(ymdb_10['species'], "Saccharomyces cerevisiae")
            self.assertEqual(ymdb_10['wikipedia'], None)

            file_name = self.output_directory + '/' + 'YMDB00003.json'
            with open (file_name, 'r') as f:
                data = json.load(f)
            self.assertEqual(data['ymdb_id'], "YMDB00003")
            self.assertEqual(data['name'], "Urea")
            self.assertEqual(data['state'], "Solid")


        elif self.source == 'ecmdb':
            ecmdb_5 = self.collection.find({"m2m_id": "M2MDB000005"})[0]
            self.assertEqual(ecmdb_5['accession'], "ECMDB00023")
            self.assertEqual(ecmdb_5['name'], "3-Hydroxyisobutyric acid")
            self.assertEqual(ecmdb_5['chemical_formula'], "C4H8O3")

            ecmdb_10 = self.collection.find({"m2m_id": "M2MDB000010"})[0]
            self.assertEqual(ecmdb_10['accession'], "ECMDB00034")
            self.assertEqual(ecmdb_10['name'], "Adenine")
            self.assertEqual(ecmdb_10['chemical_formula'], "C5H5N5")

            file_name = self.output_directory + '/' + 'M2MDB000003.json'
            with open (file_name, 'r') as f:
                data = json.load(f)
            self.assertEqual(data['accession'], "ECMDB00014")
            self.assertEqual(data['name'], "Deoxycytidine")
            self.assertEqual(data['wikipedia'], "Deoxycytidine")

        else:
            print("Database source has to be 'ecmdb' or 'ymdb'")


