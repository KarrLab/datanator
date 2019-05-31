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
        cls.collection_str = 'kegg_reaction_class'
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
    def test_parse_root_json(self):
        names = self.src.parse_root_json()
        self.assertTrue(len(names) > 100)
        self.assertTrue(names[:6], ['RC00001','RC00002','RC00003','RC00004','RC00005','RC00006'])

    # @unittest.skip('passed')
    def test_parse_rc_multiline(self):
        lines = ['REACTION    R00021 R00093 R00114 R00243 R00248 R00250 R00258 R00279, R00667']
        rxn = self.src.parse_rc_multiline(lines)
        self.assertEqual(rxn[0], 'R00021')
        self.assertEqual(rxn[-1], 'R00667')
        line = ['DEFINITION  C1c-C5a:N1a+*-*+O5a:C1b+C6a-C1b+C6a', 'S1a-S3x:*-S3x:C1b-C1x']
        definition = self.src.parse_rc_multiline(line)
        self.assertEqual(definition[1], 'S1a-S3x:*-S3x:C1b-C1x')

    def test_parse_rc_orthology(self):
        lines = ['ORTHOLOGY   K00260  glutamate dehydrogenase [EC:1.4.1.2]',
        'K00261  glutamate dehydrogenase (NAD(P)+) [EC:1.4.1.3]',
        'K00816  kynurenine---oxoglutarate transaminase / cysteine-S-conjugate beta-lyase / glutamine---phenylpyruvate transaminase [EC:2.6.1.7 4.4.1.13 2.6.1.64]']
        ko_id, names = self.src.parse_rc_orthology(lines)
        self.assertEqual(ko_id, ['K00260', 'K00261', 'K00816'])
        self.assertEqual(names[0], 'glutamate dehydrogenase')
        self.assertEqual(names[2], ['kynurenine---oxoglutarate transaminase',
                                    'cysteine-S-conjugate beta-lyase',
                                    'glutamine---phenylpyruvate transaminase'])


    @unittest.skip('hold up a min')
    def test_load_content(self):
        self.src.load_content()
        col = self.collection
        self.assertEqual(col.count_documents({}), 14)
