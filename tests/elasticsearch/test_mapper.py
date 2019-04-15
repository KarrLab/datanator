'''Tests of mapper

:Author: Zhouyang Lian <zhouyang.lian@familian.life>
:Author: Jonathan <jonrkarr@gmail.com>
:Date: 2019-04-02
:Copyright: 2019, Karr Lab
:License: MIT
'''

import unittest
import shutil
import tempfile
from datanator.elasticsearch import mapper
import os
import json


class TestMapMongo(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.MongoDB = 'mongodb://mongo:27017/'
        cls.source = 'uniprot'
        cls.elastic = 'http://elasticsearch:9200'
        cls.src = mapper.Mapper(
            cls.elastic, cls.MongoDB, cls.source, cls.cache_dirname, verbose=True)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    @unittest.skip("test_elasticsearch")
    def test_elasticsearch(self):
        self.assertEqual(self.src.elastic_con().status_code, 200)

    def test_mapto(self):
        self.src.map_to()