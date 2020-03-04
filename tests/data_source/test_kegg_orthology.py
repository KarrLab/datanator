import unittest
from unittest.mock import patch, mock_open
from datanator.data_source import kegg_orthology
import tempfile
import shutil
import json
import os
import datanator.config.core

class TestKeggOrthology(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'test'
        username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        cls.MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
        port = datanator.config.core.get_config()['datanator']['mongodb']['port']
        replSet = datanator.config.core.get_config()['datanator']['mongodb']['replSet']
        cls.collection_str = 'kegg_orthology_new'
        cls.src = kegg_orthology.KeggOrthology(
                                    cls.cache_dirname, cls.MongoDB, cls.db, 
                                    replicaSet=replSet, verbose=True, max_entries=20,
                                    username = username, password = password)
        cls.client, cls.db, cls.collection = cls.src.con_db(cls.collection_str)
        path = os.path.join(cls.cache_dirname, cls.collection_str)
        os.makedirs(path, exist_ok=True)
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
    
    @unittest.skip('passed')
    def test_download_ko(self):
        file_name = 'K03014'
        self.src.download_ko(file_name + '.txt')
        path_to_file = os.path.join(self.cache_dirname, self.collection_str)
        self.assertTrue(os.path.exists(path_to_file+'/'+file_name+'.txt'))

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

    @unittest.skip('hold up a min')
    def test_parse_definition(self):
        line = 'Definition fructose-bisphosphate aldolase, class II [EC:4.1.2.13]'
        name_list, ec_list = self.src.parse_definition(line)
        self.assertEqual(['fructose-bisphosphate aldolase, class II'], name_list)
        self.assertEqual(['4.1.2.13'], ec_list)

    # @unittest.skip('passed')
    def test_parse_pathway_disease(self):
        input_string = '''PATHWAY     ko00260  Glycine, serine and threonine metabolism
            ko00270  Cysteine and methionine metabolism
            ko00300  Lysine biosynthesis
            ko01100  Metabolic pathways
            ko01110  Biosynthesis of secondary metabolites
            ko01120  Microbial metabolism in diverse environments
            ko01130  Biosynthesis of antibiotics
            ko01230  Biosynthesis of amino acids'''
        with patch('builtins.open', mock_open(read_data=input_string)) as m:
            with open('foo') as h:
                lines = h.readlines()
        result = self.src.parse_pathway_disease(lines)
        exp = [{'kegg_pathway_code': 'ko00260', 'pathway_description': 'Glycine, serine and threonine metabolism'}, 
        {'kegg_pathway_code': 'ko00270', 'pathway_description': 'Cysteine and methionine metabolism'}, 
        {'kegg_pathway_code': 'ko00300', 'pathway_description': 'Lysine biosynthesis'}, 
        {'kegg_pathway_code': 'ko01100', 'pathway_description': 'Metabolic pathways'}, 
        {'kegg_pathway_code': 'ko01110', 'pathway_description': 'Biosynthesis of secondary metabolites'}, 
        {'kegg_pathway_code': 'ko01120', 'pathway_description': 'Microbial metabolism in diverse environments'}, 
        {'kegg_pathway_code': 'ko01130', 'pathway_description': 'Biosynthesis of antibiotics'}, 
        {'kegg_pathway_code': 'ko01230', 'pathway_description': 'Biosynthesis of amino acids'}]
        self.assertEqual(result, exp)

    def test_parse_gene(self):
        input_string = '''GENES       ENR: H650_03765
            CRO: ROD_12481
            KSA: C813_04375
            KOR: AWR26_00195
            KRD: A3780_24485
            LEH: C3F35_03305
            LAZ: A8A57_10860
            PSHI: SAMEA2665130_2199(fabG_4)
            XFA: XF_0319
            XFT: PD_0262(fabG)
            XFM: Xfasm12_0283
            XFN: XfasM23_0254
            XFF: XFLM_06665(fabG)
            XFL: P303_10495
            XFS: D934_02370
            XFH: XFHB_01355
            XTW: AB672_01530
            XCC: XCC2294(phbB) XCC3355(fabG)'''
        with patch('builtins.open', mock_open(read_data=input_string)) as m:
            with open('foo') as h:
                lines = h.readlines()
        result = self.src.parse_gene(lines)
        self.assertEqual(result[-1], {'organism': 'XCC', 'genetic_info': [{'locus_id': 'XCC2294', 'gene_id': 'phbB'}, {'locus_id': 'XCC3355', 'gene_id': 'fabG'}]})        