import unittest
from datanator.data_source import kegg_reaction
from datanator.core import query_nosql
import tempfile
import shutil
import pymongo
import json
import os

class TestKeggOrthology(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        # cls.cache_dirname = './datanator/data_source/cache/'
        cls.db = 'test'
        cls.MongoDB = 'mongodb://mongo:27017/'
        cls.collection_str = 'kegg_reaction'
        cls.src = kegg_reaction.KeggReaction(
            cls.cache_dirname, cls.MongoDB, cls.db, replicaSet=None, verbose=True, max_entries=20)
        cls.client, cls.db, cls.collection = cls.src.con_db(cls.collection_str)
        path = os.path.join(cls.cache_dirname, cls.collection_str)
        os.makedirs(path, exist_ok=True)
        cls.data = {
            "name":"br08204",
            "children":[
            {
                "name":"1. Oxidoreductase reactions",
                "children":[
                {
                    "name":"1.1  Acting on the CH-OH group of donors",
                    "children":[
                    {
                        "name":"1.1.1  With NAD+ or NADP+ as acceptor",
                        "children":[
                        {
                            "name":"1.1.1.1",
                            "children":[
                            {
                                "name":"RC00050"
                            },
                            {
                                "name":"RC00087"
                            },
                            {
                                "name":"RC00088"
                            },
                            {
                                "name":"RC00099"
                            },
                            {
                                "name":"RC00116"
                            },
                            {
                                "name":"RC00649"
                            },
                            {
                                "name":"RC00955"
                            },
                            {
                                "name":"RC01734"
                            },
                            {
                                "name":"RC02273"
                            }
                            ]
                        } ] } ] } ] } ] }

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.client.close()

    # @unittest.skip('passed')
    def test_extract_values(self):
        data = json.dumps(self.data)
        loaded_data = json.loads(data)
        name_list = self.src.extract_values(loaded_data, 'name')
        self.assertEqual(name_list[0], 'br08204')
        self.assertEqual(name_list[-1][:7], 'RC02273')
    
    # @unittest.skip('passed')
    def test_download_rxn_cls(self):
        file_name = 'RC02273'
        self.src.download_rxn_cls(file_name + '.txt')
        path_to_file = os.path.join(self.cache_dirname, self.collection_str)
        self.assertTrue(os.path.exists(path_to_file+'/'+file_name+'.txt'))

    # @unittest.skip('passed')
    def test_pars_root_json(self):
        names = self.src.parse_root_json()
        self.assertTrue(len(names) > 100)
        self.assertTrue(names[:6], ['RC00001','RC00002','RC00003','RC00004','RC00005','RC00006'])

    @unittest.skip('passed')
    def test_parse_ko_txt(self):
        file_name = 'K03014'
        self.src.download_ko(file_name + '.txt')
        doc = self.src.parse_ko_txt(file_name + '.txt')
        self.assertEqual(doc['kegg_orthology_id'], 'K03014')
        self.assertEqual(doc['gene_ortholog'][0], {'organism': 'HSA', 'gene_id': '5435(POLR2F)'})
        self.assertEqual(doc['gene_ortholog'][-1], {'organism': 'VG', 'gene_id': ['22220299(C147L)', '9887894(crov491)']})
        self.assertEqual(doc['reference'][0], {'namespace': 'PMID', 'id': '19896367'})

        file_name3 = 'K01810'
        self.src.download_ko(file_name3 + '.txt')
        doc3 = self.src.parse_ko_txt(file_name3+'.txt')
        self.assertEqual(doc3['kegg_orthology_id'], 'K01810')
        self.assertEqual(doc3['gene_ortholog'][0], {'organism': 'HSA', 'gene_id': '2821(GPI)'})
        self.assertEqual(doc3['gene_ortholog'][-1], {'organism': 'LOKI', 'gene_id': ['Lokiarch_08040(pgi_1)', 'Lokiarch_21890(pgi_2)']})
        self.assertEqual(doc3['reference'][0], {'namespace': 'PMID', 'id': '2387591'})

    @unittest.skip('hold up a min')
    def test_load_content(self):
        self.src.load_content()
        col = self.collection
        self.assertEqual(col.count_documents({}), 20)
        self.assertEqual(col.count_documents({'kegg_orthology_id': 'K00001'}), 1)
        doc = col.find_one({'kegg_orthology_id': 'K00001'})
        self.assertEqual(doc['gene_name'],  ["E1.1.1.1", "adh"])