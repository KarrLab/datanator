import unittest
from datanator.data_source import taxon_tree
import tempfile
import shutil
import os
import json

class TestTaxonTree(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        # cls.cache_dirname = './datanator/data_source/cache/'
        cls.db = 'test'
        cls.MongoDB = 'mongodb://mongo:27017/'
        cls.collection_str = 'taxon_tree'
        cls.src = taxon_tree.TaxonTree(
            cls.cache_dirname, cls.MongoDB, cls.db, replicaSet=None, verbose=True, max_entries=20)
        cls.path = os.path.join(cls.cache_dirname, cls.collection_str)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.src.client.close()

    @unittest.skip('passed')
    def test_download_dump(self):
    	noi = 'division.dmp'
    	my_file = os.path.join(self.path, noi)
    	self.src.download_dump()
    	self.assertTrue(os.path.isfile(my_file))

    @unittest.skip('passed')
    def test_parse_fullname_line(self):
    	line1 = '1431752	|	Anaxipha thomasi	|	allotype	|	SWRC Ae_151	|'
    	line2 = '1969370	|	131567 2157 1935183 1706441 	|'
    	line3 = '2026747	|	Candidatus Heimdallarchaeota archaeon	|		|		|		|		|		|	Candidatus Heimdallarchaeota	|		|	Archaea	|'

    	self.assertEqual(self.src.parse_fullname_line(line1), ['1431752', 'Anaxipha thomasi', 'allotype', 'SWRC Ae_151'])
    	self.assertEqual(self.src.parse_fullname_line(line2), ['1969370', '131567 2157 1935183 1706441 '])
    	self.assertEqual(self.src.parse_fullname_line(line3)[1], 'Candidatus Heimdallarchaeota archaeon' )
    	self.assertEqual(self.src.parse_fullname_line(line3)[2], '' )
    	self.assertEqual(self.src.parse_fullname_line(line3)[-1], 'Archaea' )

    @unittest.skip('passed')
    def test_parse_taxid_line(self):
    	line1 = '1841598	|	131567 2157 1935183 1936272 	|'
    	self.assertEqual(self.src.parse_taxid_line(line1), ['131567', '2157', '1935183', '1936272'])

    @unittest.skip('passed')
    def test_parse_fullname_taxid(self):
    	self.src.parse_fullname_taxid()
    	doc = self.src.collection.find_one( {'tax_id': '1935183'})
    	self.assertEqual(doc['anc_id'], ['131567', '2157'])

    @unittest.skip('passed')
    def test_parse_nodes(self):
    	self.src.parse_nodes()

    @unittest.skip('passed')
    def test_parse_division(self):
    	self.src.parse_division()

    @unittest.skip('passed')
    def test_parse_names(self):
    	self.src.parse_names()

    @unittest.skip('passed')
    def test_parse_gencode(self):
    	self.src.parse_gencode()



