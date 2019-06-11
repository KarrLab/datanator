import unittest
from datanator.core import query_nosql
import tempfile
import shutil
from datanator.util import server_util

class TestQueryNoSQL(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'datanator'
        cls.MongoDB = 'mongodb://mongo:27017/'
        cls.collection_str = 'ecmdb'
        cls.src = query_nosql.DataQuery(
            cache_dirname=cls.cache_dirname, MongoDB=cls.MongoDB, replicaSet=None, db=cls.db,
                 verbose=True, max_entries=20)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

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


class TestQueryMetabolitesMeta(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'datanator'
        cls.MongoDB = 'mongodb://mongo:27017/'
        cls.src = query_nosql.QueryMetabolitesMeta(
            cache_dirname=cls.cache_dirname, MongoDB=cls.MongoDB, replicaSet=None, db=cls.db,
                 verbose=True, max_entries=20)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    @unittest.skip('passed')
    def test_find_synonyms(self):
        c = "ATP"
        compounds = ["ATP", "Oxygen", '']
        ran = 'something'
        rxn, syn = self.src.find_synonyms(c)
        rxn_ran, syn_ran = self.src.find_synonyms(ran)
        rxns, syns = self.src.find_synonyms(compounds)
        self.assertTrue("5'-ATP" in syn[c])
        self.assertTrue("Oxygen molecule" in syns['Oxygen'])
        self.assertEqual(None, syns['synonyms'])
        self.assertTrue('does not exist' in syn_ran[ran])

        empty = ''
        rxn, syn = self.src.find_synonyms(empty)
        self.assertEqual(syn, {'synonyms': None})

    @unittest.skip('passed')
    def test_find_rxn_by_participant(self):
        substrates = ["Undecanal", 'H2O', 'NADP+']
        products = ['Undecanoate', "NADPH", 'H+']
        ids = self.src.find_rxn_by_participant(substrates, products)
        self.assertTrue(28926 in ids)


class TestQuerySabio(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'datanator'
        config_file = '/root/host/karr_lab/datanator/.config/config.ini'
        username, password, server, port = server_util.ServerUtil(config_file = config_file).get_admin_config()
        cls.MongoDB = server
        cls.username = username
        cls.password = password
        cls.src = query_nosql.QuerySabio(
            cache_dirname=cls.cache_dirname, MongoDB=cls.MongoDB, db=cls.db,
                 verbose=True, max_entries=20, username = cls.username, password = cls.password)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    @unittest.skip('passed')
    def test_find_reaction_participants(self):
        _id = [31, 32, 33, 34]
        rxns = self.src.find_reaction_participants(_id)
        self.assertEqual(rxns[0], {'substrates': ["2-Hydroxybutyrate", "Riboflavin-5-phosphate"], 
                        'products': ['2-Oxobutyrate', 'Reduced FMN']})

        self.assertEqual(rxns[3], {'substrates': ['Riboflavin-5-phosphate', '4-Chloromandelate'],
                                    'products': ['Reduced FMN', '4-Chloro-2-Oxobenzeneacetic acid'] } )

    def test_get_kinlawid_by_inchi(self):
        inchi = ['InChI=1S/C8H16O3/c1-2-3-4-5-6-7(9)8(10)11/h7,9H,2-6H2,1H3,(H,10,11)',
        'InChI=1S/C17H21N4O9P/c1-7-3-9-10(4-8(7)2)21']
        kinlaw_id = self.src.get_kinlawid_by_inchi(inchi)
        self.assertTrue(28 in kinlaw_id)
