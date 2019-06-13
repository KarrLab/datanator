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
    def test_get_metabolite_synonyms(self):
        c = "ATP"
        compounds = ["ATP", "Oxygen", '']
        ran = 'something'
        rxn, syn = self.src.get_metabolite_synonyms(c)
        rxn_ran, syn_ran = self.src.get_metabolite_synonyms(ran)
        rxns, syns = self.src.get_metabolite_synonyms(compounds)
        self.assertTrue("5'-ATP" in syn[c])
        self.assertTrue("Oxygen molecule" in syns['Oxygen'])
        self.assertEqual(None, syns['synonyms'])
        self.assertTrue('does not exist' in syn_ran[ran])

        empty = ''
        rxn, syn = self.src.get_metabolite_synonyms(empty)
        self.assertEqual(syn, {'synonyms': None})

    # @unittest.skip('passed')
    def test_get_metabolite_inchi(self):
        compounds = ['atp']
        inchi = self.src.get_metabolite_inchi(compounds)
        self.assertEqual(inchi[0], 'InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20')


class TestQuerySabio(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'datanator'
        config_file = '/root/host/karr_lab/datanator/.config/config.ini'
        username, password, server, port = server_util.ServerUtil(config_file = config_file).get_user_config()
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

    @unittest.skip('passed')
    def test_get_kinlawid_by_inchi_slow(self):
        inchi = ['InChI=1S/C8H8O3/c9-7(8(10)11)6-4-2-1-3-5-6/h1-5,7,9H,(H,10,11)/t7-/m0/s1',
                'InChI=1S/C17H21N4O9P/c1-7-3-9-10(4-8(7)2)21(15-13(18-9)16(25)20-17(26)19-15)5-11(22)14(24)12(23)6-30-31(27,28)29/h3-4,11-12,14,22-24H,5-6H2,1-2H3,(H,20,25,26)(H2,27,28,29)/t11-,12+,14-/m0/s1',
                'InChI=1S/C8H6O3/c9-7(8(10)11)6-4-2-1-3-5-6/h1-5H,(H,10,11)/p-1']
        kinlaw_id = self.src.get_kinlawid_by_inchi_slow(inchi)
        print(kinlaw_id)
        self.assertTrue(9 in kinlaw_id)

    @unittest.skip('passed')
    def test_get_kinlawid_by_inchi(self):
        inchi = ['InChI=1S/C8H8O3/c9-7(8(10)11)6-4-2-1-3-5-6/h1-5,7,9H,(H,10,11)/t7-/m0/s1',
                'InChI=1S/C17H21N4O9P/c1-7-3-9-10(4-8(7)2)21(15-13(18-9)16(25)20-17(26)19-15)5-11(22)14(24)12(23)6-30-31(27,28)29/h3-4,11-12,14,22-24H,5-6H2,1-2H3,(H,20,25,26)(H2,27,28,29)/t11-,12+,14-/m0/s1',
                'InChI=1S/C8H6O3/c9-7(8(10)11)6-4-2-1-3-5-6/h1-5H,(H,10,11)/p-1']
        rxn = self.src.get_kinlawid_by_inchi(inchi)
        self.assertTrue(9 in rxn)
        self.assertTrue(21016 in rxn)

    @unittest.skip('passed')
    def test_get_kinlawid_by_rxn(self):
        substrates = ['InChI=1S/C8H8O3/c9-7(8(10)11)6-4-2-1-3-5-6/h1-5,7,9H,(H,10,11)/t7-/m0/s1',
                    'InChI=1S/C17H21N4O9P/c1-7-3-9-10(4-8(7)2)21(15-13(18-9)16(25)20-17(26)19-15)5-11(22)14(24)12(23)6-30-31(27,28)29/h3-4,11-12,14,22-24H,5-6H2,1-2H3,(H,20,25,26)(H2,27,28,29)/t11-,12+,14-/m0/s1']
        products = ['InChI=1S/C8H6O3/c9-7(8(10)11)6-4-2-1-3-5-6/h1-5H,(H,10,11)/p-1']
        _id = self.src.get_kinlawid_by_rxn(substrates,products)
        self.assertTrue(21016 in _id)