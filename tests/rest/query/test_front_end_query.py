import unittest
from datanator.rest.query import front_end_query
import tempfile
import shutil

class TestQueryFrontEnd(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.src = front_end_query.QueryFrontEnd()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    @unittest.skip('passed')
    def test_string_query(self):
        string = '2-Ketobutyric acid (alpha ketobutyric acid)'
        results = self.src.string_query(string)
        self.assertEqual(results[0]['m2m_id'], 'M2MDB001633')

    @unittest.skip('passed')
    def test_inchi_query_metabolite(self):
    	inchi = 'InChI=1S/C3H8O10P2/c4-3(5)2(13-15(9,10)11)1-12-14(6,7)8'
    	results = self.src.inchi_query_metabolite(inchi)
    	self.assertEqual(results[0]['m2m_id'], 'M2MDB000326')
    	self.assertEqual(results[1]['ymdb_id'], 'YMDB00671')

    def test_inchi_query_organism(self):
    	inchi = 'InChI=1S/C3H8O10P2/c4-3(5)2(13-15(9,10)11)1-12-14(6,7)8'
    	organism = 'Candidatus Pacearchaeota archaeon'
    	results = self.src.inchi_query_organism(inchi, organism)
    	self.assertEqual(results[0]['m2m_id'], 'M2MDB000326')
    	self.assertEqual(results[0]['taxon_distance'], -1)