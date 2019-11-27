import unittest
from karr_lab_aws_manager.elasticsearch_kl import batch_load
from datanator_query_python.config import config
import tempfile
import shutil
import requests

class TestMongoToES(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dir = tempfile.mkdtemp()
        cls.src = batch_load.MongoToES(profile_name='es-poweruser', credential_path='~/.wc/third_party/aws_credentials',
                config_path='~/.wc/third_party/aws_config', elastic_path='~/.wc/third_party/elasticsearch.ini',
                cache_dir=cls.cache_dir, service_name='es', index='test', max_entries=float('inf'), verbose=True)
        cls.url = cls.src.es_endpoint + '/' + cls.src.index
        requests.delete(cls.url, auth=cls.src.awsauth)
        conf = config.Config()
        cls.username = conf.USERNAME
        cls.password = conf.PASSWORD
        cls.server = conf.SERVER
        cls.authDB = conf.AUTHDB
        cls.db = 'datanator'

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dir)
        requests.delete(cls.url, auth=cls.src.awsauth)

    def test_connection(self):
        result = self.src.client.list_domain_names()
        self.assertEqual(result['ResponseMetadata']['HTTPStatusCode'], 200)
        self.assertTrue('datanator-elasticsearch' in self.src.es_endpoint)

    def test_data_from_mongo(self):
        count, _ = self.src.data_from_mongo_protein(self.server, self.db, self.username, 
                                                    self.password, authSource=self.authDB)
        self.assertTrue(count >= 1000)

    def test_data_from_metabolite(self):
        _, count_0, _, count_1 = self.src.data_from_mongo_metabolite(self.server, self.db, self.username, 
                                            self.password, authSource=self.authDB)
        self.assertTrue(count_0 >= 1000)
        self.assertTrue(count_1 >= 1000)

    def test_data_from_metabolites_meta(self):
        doc = self.src.data_from_mongo_metabolites_meta(self.server, self.db, self.username,
                                                            self.password, authSource=self.authDB)
        result = []
        for i in range(5):
            result.append(doc)
        self.assertEqual(len(result), 5)