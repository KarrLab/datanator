import unittest
from datanator.data_source import taxon_tree
import datanator.config.core
import tempfile
import shutil
import os
import json


class TestTaxonTree(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'test'
        username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
        port = datanator.config.core.get_config()['datanator']['mongodb']['port']
        replSet = datanator.config.core.get_config()['datanator']['mongodb']['replSet']
        cls.collection_str = 'taxon_tree'
        cls.src = taxon_tree.TaxonTree(
            cls.cache_dirname, MongoDB, cls.db, replicaSet=None, 
            verbose=True, max_entries=10, username = username, password = password)
        cls.path = os.path.join(cls.cache_dirname, cls.collection_str)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.src.client.close()

    # @unittest.skip('passed')
    def test_download_dump(self):
        noi = 'division.dmp'
        my_file = os.path.join(self.path, noi)
        self.src.download_dump()
        self.assertTrue(os.path.isfile(my_file))

    # @unittest.skip('passed')
    def test_parse_fullname_line(self):
        line1 = '1936272    |   Candidatus Heimdallarchaeota    |   cellular organisms; Archaea; Asgard group;  |'
        line2 = '2012493    |   Candidatus Heimdallarchaeota archaeon B3_Heim   |   cellular organisms; Archaea; Asgard group; Candidatus Heimdallarchaeota;    |'
        line3 = '1935183    |   Asgard group    |   cellular organisms; Archaea;    |'

        self.assertEqual(self.src.parse_fullname_line(line1), [
                         '1936272', 'Candidatus Heimdallarchaeota', ['cellular organisms', 'Archaea', 'Asgard group']])
        self.assertEqual(self.src.parse_fullname_line(line2)[:2], [
                         '2012493', 'Candidatus Heimdallarchaeota archaeon B3_Heim'])
        self.assertEqual(self.src.parse_fullname_line(line3)
                         [1], 'Asgard group')
        self.assertEqual(self.src.parse_fullname_line(line3)[
                         2], ['cellular organisms', 'Archaea'])

    # @unittest.skip('passed')
    def test_parse_taxid_line(self):
        line1 = '1841596\t|\t131567 2157 1935183 1936272 \t|\n'
        self.assertEqual(self.src.parse_taxid_line(line1), [
                         '131567', '2157', '1935183', '1936272'])

    # @unittest.skip('passed')
    def test_parse_fullname_taxid(self):
        self.src.parse_fullname_taxid()
        doc = self.src.collection.find_one({'tax_id': 1935183})
        self.assertEqual(doc['anc_id'], [131567, 2157])

    # @unittest.skip('passed')
    def test_parse_nodes(self):
        self.src.parse_nodes()
        doc = self.src.collection.find_one({'tax_id': 1})
        self.assertEqual(doc['tax_name'], 'root')
        self.assertEqual(doc['division_id'], 8)

    # @unittest.skip('passed')
    def test_parse_division(self):
        self.src.parse_division()

    # @unittest.skip('passed')
    def test_parse_names(self):
        self.src.parse_names()

    # @unittest.skip('passed')
    def test_parse_gencode(self):
        self.src.parse_gencode()

    def test_load_content(self):
        self.src.load_content()
