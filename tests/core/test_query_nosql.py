import unittest
from datanator.core import query_nosql
import tempfile
import shutil
import pprint

class TestQueryNoSQL(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'test'
        cls.MongoDB = 'mongodb://mongo:27017/'
        cls.src = query_nosql.DataQuery(
            cls.cache_dirname, cls.MongoDB, None, cls.db, verbose=True, max_entries=20)
        cls.collection_str = 'ecmdb'
        cls.client, cls.db, cls.collection_obj = cls.src.con_db(
            cls.collection_str)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        # cls.client.drop_database(cls.db)
        cls.client.close()

    @unittest.skip('skip to testing for h1_hesc')
    def test_doc_feeder(self):
        query = {'m2m_id': {
            '$in': ["M2MDB000004", "M2MDB000005", "M2MDB000006"]}}
        col = self.src.doc_feeder(self.collection_str, query=query)
        for doc in col:
            if doc['m2m_id'] == "M2MDB000005":
                self.assertEqual(doc['accession'], "ECMDB00023")
            elif doc['accession'] == "ECMDB00019":
                self.assertEqual(doc['m2m_id'], "M2MDB000004")

        taxon_range = list(range(100, 20000))
        null = None
        query = {'taxon': {'$in': taxon_range}, 
                 'taxon_wildtype': {'$eq': 1}, 
                 'temperature': {'$lte':45, '$gte': 20},
                 'ph': {'$lte': 8, '$gte': 6},
                 'parameter.value': {'$ne': null },
                 'parameter.sbo_type': {'$in': [25,27] } }

        # only return these fileds in projection
        projection = {'kinlaw_id':1, 'resource': 1, 'enzymes.enzyme': 1, 
                    'enzymes.subunit': 1, 'parameter': 1}

        collection = self.src.doc_feeder('sabio_rk', query = query, projection=projection)
        with self.assertRaises(KeyError):
            next(collection)['reaction_participant']


    '''Testing queries in h1_hesc
    '''
    @unittest.skip('skip to testing for h1_hesc')
    def test_doc_feeder(self):
        target = 'C8H8O3/c9-7(8(10)11)6-4-2-1-3-5-6'
        query = {'$or': [ {'reaction_participant.substrate.structure.inchi_connectivity': target }, 
                          {'reaction_participant.product.structure.inchi_connectivity': target }  ] }
        projection = {'reaction_participant.substrate.sabio_compound_id': 1,
                     'reaction_participant.product.sabio_compound_id': 1 }
        collection = self.src.doc_feeder('sabio_rk', query = query, projection=projection)