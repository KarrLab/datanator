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
        cls.meta_loc = 'datanator'
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


    # def test_fill_metabolite_fields(self):
    #     dict_list = self.src.fill_metabolite_fields(
    #     	fields = ['m2m_id', 'inchi'], collection_src = 'ecmdb', collection_des = 'metabolites_meta')
    #     self.assertEqual(
    #         dict_list[0]['inchi'], 'InChI=1S/C4H6O3/c1-2-3(5)4(6)7')

    def test_load_content(self):
        self.src.load_content()
        client = pymongo.MongoClient(
            self.MongoDB, username = self.username, 
            password = self.password)
        meta_db = client[self.meta_loc]
        collection = meta_db['metabolites_meta']
        cursor = collection.find_one({'inchi': 'InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2'})
        self.assertEqual(cursor['inchi_hashed'], '2b24980b46fd11b024dce40fef6ff4453f0e30682adb1cb2b3fb60e1')
        
        client.close()