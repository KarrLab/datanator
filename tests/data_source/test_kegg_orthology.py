import unittest
from datanator.data_source import kegg_orthology
import tempfile
import shutil
import pymongo
import json

class TestKeggOrthology(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'test'
        cls.MongoDB = 'mongodb://mongo:27017/'
        cls.collection_str = 'kegg_orthology'
        cls.src = kegg_orthology.KeggOrthology(
            cls.cache_dirname, cls.MongoDB, cls.db, replicaSet=None, verbose=True, max_entries=20)
        cls.client, cls.db, cls.collection = cls.src.con_db(cls.collection_str)
        cls.data = {
            "name":"ko00001",
            "children":[
            {
                "name":"09100 Metabolism",
                "children":[
                {
                    "name":"09101 Carbohydrate metabolism",
                    "children":[
                    {
                        "name":"00010 Glycolysis \/ Gluconeogenesis [PATH:ko00010]",
                        "children":[
                        {
                            "name":"K00844  HK; hexokinase [EC:2.7.1.1]"
                        },
                        {
                            "name":"K12407  GCK; glucokinase [EC:2.7.1.2]"
                        },
                        {
                            "name":"K00845  glk; glucokinase [EC:2.7.1.2]"
                        },
                        {
                            "name":"K01810  GPI, pgi; glucose-6-phosphate isomerase [EC:5.3.1.9]"
                        }
                        ]
                    }]
                }],
                "name":"00020 Citrate cycle (TCA cycle) [PATH:ko00020]",
                "children":[
                {
                    "name":"K01647  CS, gltA; citrate synthase [EC:2.3.3.1]"
                },
                {
                    "name":"K01648  ACLY; ATP citrate (pro-S)-lyase [EC:2.3.3.8]"
                },
                {
                    "name":"K15230  aclA; ATP-citrate lyase alpha-subunit [EC:2.3.3.8]"
                },
                {
                    "name":"K15231  aclB; ATP-citrate lyase beta-subunit [EC:2.3.3.8]"
                }
                ]
            }]}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.client.close()

    def test_con_db(self):
        self.assertNotEqual(self.collection, 'Server not available')

    def test_extract_values(self):
        data = json.dumps(self.data)
        loaded_data = json.loads(data)
        name_list = self.src.extract_values(loaded_data, 'name')
        self.assertEqual(name_list[0], 'ko00001')
        self.assertEqual(name_list[-1][:6], 'K15231')
    
    def test_parse_ko_txt(self):
        pass

