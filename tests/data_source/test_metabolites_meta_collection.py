import unittest
from datanator.data_source import metabolites_meta_collection
import tempfile
import shutil


class TestMongoUtil(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'datanator'
        cls.MongoDB = 'mongodb://mongo:27017/'
        cls.collection_str = 'ecmdb'
        cls.src = metabolites_meta_collection.MetabolitesMeta(cache_dirname=cls.cache_dirname,
                                                              MongoDB=cls.MongoDB, replicaSet='rs0', db=cls.db,
                                                              collection_str=cls.collection_str, verbose=True, max_entries=float('inf'))
        # cls.client, cls.db, cls.collection_obj = cls.src.con_db(
        #     cls.collection_str)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        # cls.client.drop_database(cls.db)
        cls.src.client_obj.close()

    @unittest.skip('passed')
    def test_parse_inchi(self):
        s = """InChI=1S/C9H13N3O4/c10-7-1-2-12(9(15)11-7)8-3-5(14)6(4-13)16-8/h1-2,
        5-6,8,13-14H,3-4H2,(H2,10,11,15)/t5-, 6+,8+/m0/s1"""
        self.assertEqual(self.src.parse_inchi(
            s), 'InChI=1S/C9H13N3O4/c10-7-1-2-12(9(15)11-7)8-3-5(14)6(4-13)16-8')
        s1 = "InChI=1S/H2O/h1H2"
        self.assertEqual(self.src.parse_inchi(s1), "InChI=1S/H2O")

    @unittest.skip('passed')
    def test_find_metabollite_inchi(self):
        doc = self.src.collection_obj.find_one({'name': 'Ureidopropionic acid'})
        inchi = self.src.find_metabolite_inchi(doc)
        self.assertEqual(inchi, 'InChI=1S/C4H8N2O3/c5-4(9)6-2-1-3(7)8')
        doc1 = {}
        self.assertEqual(self.src.find_metabolite_inchi(doc1), 
        		'No key named "inchi" in given document')

    def test_find_rxn_id(self):
    	inchi = 'InChI=1S/C9H13N3O4/c10-7-1-2-12(9(15)11-7)8-3-5(14)6(4-13)16-8'
    	self.src.find_rxn_id(inchi)