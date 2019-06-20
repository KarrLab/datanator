import unittest
from datanator.data_source import metabolites_meta_collection
import datanator.config.core
import pymongo
import tempfile
import shutil


class TestMetabolitesMeta(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'datanator'
        cls.meta_loc = 'test'
        cls.username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        cls.password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        cls.MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
        port = datanator.config.core.get_config()['datanator']['mongodb']['port']
        replSet = datanator.config.core.get_config()['datanator']['mongodb']['replSet']
        cls.src = metabolites_meta_collection.MetabolitesMeta(cache_dirname=cls.cache_dirname,
                                                              MongoDB=cls.MongoDB,  db=cls.db,
                                                              verbose=True, max_entries=20, username = cls.username,
                                                              password = cls.password, meta_loc = cls.meta_loc)
        cls.client, cls.db_obj, cls.collection_obj = cls.src.con_db(cls.db)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.client.close()


    def test_get_metabolite_fields(self):
        dict_list = self.src.get_metabolite_fields(
        	fields = ['m2m_id', 'inchi'], collection_str = 'ecmdb')
        self.assertEqual(
            dict_list[0]['inchi'], 'InChI=1S/C4H6O3/c1-2-3(5)4(6)7')

    def test_load_content(self):
        self.src.load_content()
        client = pymongo.MongoClient(
            self.MongoDB, username = self.username, 
            password = self.password)
        meta_db = client[self.meta_loc]
        collection = meta_db['metabolites_meta']
        cursor = collection.find_one({'inchi': 'InChI=1S/C5H10N2O3S/c6-3(2-11)5(10)7-1-4(8)9'})
        self.assertEqual(cursor['inchi_hashed'], '95b3e8d30d2dcf0041a5add611bcc3760ca70f62055e9644bb023bde')
        
        client.close()