import unittest
from datanator.core import query_taxon_tree
import tempfile
import shutil
import configparser
import os
import json

class TestQueryTaxonTree(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'datanator'
        parser = configparser.ConfigParser(allow_no_value=True)
        parser.read(os.path.expanduser('~/.wc/datanator.ini'))
        cls.username = parser.get('mongodb', 'user')
        cls.password = parser.get('mongodb', 'password')
        cls.MongoDB = parser.get('mongodb', 'server')
        port = int(parser.get('mongodb', 'port'))
        replSet = parser.get('mongodb', 'replSet')

        cls.src = query_taxon_tree.QueryTaxonTree(
            cache_dirname=cls.cache_dirname, MongoDB=cls.MongoDB, db=cls.db,
                 verbose=True, max_entries=20, username = cls.username, password = cls.password)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    # @unittest.skip('passed')
    def test_get_all_species(self):
        result = []
        generator = self.src.get_all_species()
        for i, name in enumerate(generator):
            if i > 100:
                break
            result.append(json.loads(name))
        self.assertEqual(len(result), 101)
        self.assertEqual(result[38]['tax_name'], "'Amaranthus retroflexus' phytoplasma")

    def test_get_name_by_id(self):
        ids = [743725, 2107591]
        names = self.src.get_name_by_id(ids)
        self.assertEqual(names[0], 'Candidatus Diapherotrites')


    # @unittest.skip('passed')
    def test_get_anc_by_name(self):
        names = ['Candidatus Diapherotrites', 'Candidatus Forterrea multitransposorum CG_2015-17_Forterrea_25_41']
        result_ids, result_names = self.src.get_anc_by_name(names)
        self.assertEqual(result_ids[0], [131567, 2157, 1783276])

    # @unittest.skip('passed')
    def test_get_anc_by_id(self):
        ids = [743725, 2107591]
        result_ids, result_names = self.src.get_anc_by_id(ids)
        self.assertEqual(result_ids[0], [131567, 2157, 1783276])
        self.assertEqual(result_ids[1], [131567, 2157, 1783276, 743725, 2107589, 2107590])

    # @unittest.skip('passed')
    def test_get_common_ancestor(self):
        names = ['Candidatus Diapherotrites', 'Candidatus Forterrea multitransposorum CG_2015-17_Forterrea_25_41']
        ancestor, distances = self.src.get_common_ancestor(names[0], names[1])
        self.assertEqual(ancestor, 1783276)
        self.assertEqual(distances[0], 1)
        self.assertEqual(distances[1], 4)

    def test_get_rank(self):
        ids = [131567, 2759, 33154, 33208, 6072, 33213, 33511, 7711, 9526, 314295, 9604, 207598, 9605, 9606]
        ranks = self.src.get_rank(ids)
        self.assertEqual(ranks[3], 'kingdom')
        self.assertEqual(ranks[1], '+')