import unittest
from datanator.core import query_nosql
import tempfile
import shutil
import datanator.config.core


class TestQueryNoSQL(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'datanator'
        username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
        port = datanator.config.core.get_config()['datanator']['mongodb']['port']
        replSet = datanator.config.core.get_config()['datanator']['mongodb']['replSet']
        cls.collection_str = 'ecmdb'
        cls.src = query_nosql.DataQuery(
            cache_dirname=cls.cache_dirname, MongoDB=MongoDB, replicaSet=None, db=cls.db,
                 verbose=True, max_entries=20, username = username, password = password)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    # @unittest.skip('skip to testing for h1_hesc')
    def test_doc_feeder(self):
        query = {'m2m_id': {
            '$in': ["M2MDB000004", "M2MDB000005", "M2MDB000006"]}}
        col = self.src.doc_feeder(self.collection_str, query=query)
        for doc in col:
            if doc['m2m_id'] == "M2MDB000005":
                self.assertEqual(doc['accession'], "ECMDB00023")
            elif doc['accession'] == "ECMDB00019":
                self.assertEqual(doc['m2m_id'], "M2MDB000004")


class TestQueryMetabolitesMeta(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'datanator'
        username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
        port = datanator.config.core.get_config()['datanator']['mongodb']['port']
        replSet = datanator.config.core.get_config()['datanator']['mongodb']['replSet']
        cls.src = query_nosql.QueryMetabolitesMeta(
            cache_dirname=cls.cache_dirname, MongoDB=MongoDB, replicaSet=replSet, db=cls.db,
                 verbose=True, max_entries=20, username = username, password = password)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    # @unittest.skip('passed')
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
        inchis = self.src.get_metabolite_inchi(compounds)
        self.assertEqual(inchis[0]['inchi'], 'InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)'+
            '15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)' +
            '27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)' +
            '(H2,18,19,20)/t4-,6-,7-,10-/m1/s1')
        compound = ['β-D-Ribopyranose']
        inchi = self.src.get_metabolite_inchi(compound)
        self.assertEqual(inchi[0]['m2m_id'], 'M2MDB001173')

    def test_get_ids_from_hash(self):
        hashed_inchi_1 = '09fab91d3708097484215d419c9326290150f37e7c1bcc48a1bb4c7b'
        hasehd_inchi_2 = 'cc46a6b3360bd6a51da1ec1b3da746456c7c05ed8ff63f930452bf9f'
        result_1 = self.src.get_ids_from_hash(hashed_inchi_1)
        result_2 = self.src.get_ids_from_hash(hasehd_inchi_2)
        self.assertEqual(result_1, {'m2m_id': 'M2MDB000016', 'ymdb_id': 'YMDB00058'})
        self.assertEqual(result_2, {'m2m_id': 'M2MDB000006', 'ymdb_id': None})

    # @unittest.skip('passed')
    def test_get_metabolite_name_by_hash(self):
        compounds = ['f56e0c2c16f3a2549c65be52179ed860b7cb375e4037061d238e433d',
                    'afaca8d9351843f37d9d010c8eb15601eb385121c988ca671ac0db31']
        result = self.src.get_metabolite_name_by_hash(compounds)
        self.assertEqual(result[0], 'Unispheres Q 10')
        self.assertEqual(result[1], 'Guanosine diphosphoric acid fucose')

    # @unittest.skip('passed')
    def test_get_metabolite_hashed_inchi(self):
        compounds = ['atp', 'Pyruvic acid', 'pyruvic Acid']
        hashed_inchi = self.src.get_metabolite_hashed_inchi(compounds)
        self.assertEqual(hashed_inchi[1], 'adb60a29296db43210c3a5456fc49c9a6e53605ce271773c79598595')
        self.assertEqual(hashed_inchi[1], hashed_inchi[2])
        compound = ['methionine']
        hashed_inchi = self.src.get_metabolite_hashed_inchi(compound)
        self.assertEqual(hashed_inchi, ['394c8c92dbe2f6c8c24fba086d06658fe8408ceb352fbf3a44732282'])

    # @unittest.skip('passed')
    def test_get_metabolite_similar_compounds(self):

        compound = ['methionine']
        raw1, result1 = self.src.get_metabolite_similar_compounds(compound, num = 3, threshold = 0.6)
        self.assertTrue(list(raw1[0][0].keys())[0], '5c40a5a611421d5a2fdb8d29e9d334009d59955909ff755db7491cf4')
        raw2, result2 = self.src.get_metabolite_similar_compounds(compound, num = 10, threshold = 0.75)
        self.assertEqual(len(raw2[0][0].keys()), 1)
        raw3, result3 = self.src.get_metabolite_similar_compounds(compound, num = 10, threshold = 0.9)
        # self.assertEqual(raw3[0]['raw'], -1)
        compound = ['β-D-Ribopyranose']
        raw4, result4 = self.src.get_metabolite_similar_compounds(compound, num = 30, threshold = 0.6)
        self.assertTrue('d' not in list(result4[0].keys()))
        # self.assertTrue('Auto inducer 2' in list(result4[0].keys()))

class TestQuerySabio(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'datanator'
        username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        server = datanator.config.core.get_config()['datanator']['mongodb']['server']
        port = datanator.config.core.get_config()['datanator']['mongodb']['port']
        replSet = datanator.config.core.get_config()['datanator']['mongodb']['replSet']
        cls.MongoDB = server
        cls.username = username
        cls.password = password
        cls.src = query_nosql.QuerySabio(
            cache_dirname=cls.cache_dirname, MongoDB=cls.MongoDB, db=cls.db,
                 verbose=True, max_entries=20, username = cls.username, password = cls.password)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    # @unittest.skip('passed')
    def test_find_reaction_participants(self):
        _id = [31, 32, 33, 34]
        rxns = self.src.find_reaction_participants(_id)
        self.assertEqual(rxns[0], {'substrates': ["2-Hydroxybutyrate", "Riboflavin-5-phosphate"], 
                        'products': ['2-Oxobutyrate', 'Reduced FMN']})

        self.assertEqual(rxns[3], {'substrates': ['Riboflavin-5-phosphate', '4-Chloromandelate'],
                                    'products': ['Reduced FMN', '4-Chloro-2-Oxobenzeneacetic acid'] } )

    # @unittest.skip('passed')
    def test_get_kinlawid_by_inchi(self):
        inchi = ['InChI=1S/C8H8O3/c9-7(8(10)11)6-4-2-1-3-5-6/h1-5,7,9H,(H,10,11)/t7-/m0/s1',
                'InChI=1S/C17H21N4O9P/c1-7-3-9-10(4-8(7)2)21(15-13(18-9)16(25)20-17(26)19-15)5-11(22)14(24)12(23)6-30-31(27,28)29/h3-4,11-12,14,22-24H,5-6H2,1-2H3,(H,20,25,26)(H2,27,28,29)/t11-,12+,14-/m0/s1',
                'InChI=1S/C8H6O3/c9-7(8(10)11)6-4-2-1-3-5-6/h1-5H,(H,10,11)/p-1']
        rxn = self.src.get_kinlawid_by_inchi(inchi)
        self.assertTrue(9 in rxn)
        self.assertTrue(21016 in rxn)

    # @unittest.skip('passed')
    def test_get_kinlawid_by_rxn(self):
        substrates = ['InChI=1S/C8H8O3/c9-7(8(10)11)6-4-2-1-3-5-6/h1-5,7,9H,(H,10,11)/t7-/m0/s1',
                    'InChI=1S/C17H21N4O9P/c1-7-3-9-10(4-8(7)2)21(15-13(18-9)16(25)20-17(26)19-15)5-11(22)14(24)12(23)6-30-31(27,28)29/h3-4,11-12,14,22-24H,5-6H2,1-2H3,(H,20,25,26)(H2,27,28,29)/t11-,12+,14-/m0/s1']
        products = ['InChI=1S/C8H6O3/c9-7(8(10)11)6-4-2-1-3-5-6/h1-5H,(H,10,11)/p-1']
        _id = self.src.get_kinlawid_by_rxn(substrates,products)
        self.assertTrue(21016 in _id)


class TestQueryTaxonTree(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'datanator'
        username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        server = datanator.config.core.get_config()['datanator']['mongodb']['server']
        port = datanator.config.core.get_config()['datanator']['mongodb']['port']
        replSet = datanator.config.core.get_config()['datanator']['mongodb']['replSet']
        cls.MongoDB = server
        cls.username = username
        cls.password = password
        cls.src = query_nosql.QueryTaxonTree(
            cache_dirname=cls.cache_dirname, MongoDB=cls.MongoDB, db=cls.db,
                 verbose=True, max_entries=20, username = cls.username, password = cls.password)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    # @unittest.skip('passed')
    def test_get_all_species(self):
        result = self.src.get_all_species()
        self.assertTrue('Candidatus Heimdallarchaeota archaeon AB_125' in result)

    def test_get_name_by_id(self):
        ids = [743725, 2107591]
        names = self.src.get_name_by_id(ids)
        self.assertEqual(names[0], 'Candidatus Diapherotrites')


    # @unittest.skip('passed')
    def test_get_anc_by_name(self):
        names = ['Candidatus Diapherotrites', 'Candidatus Forterrea multitransposorum CG_2015-17_Forterrea_25_41']
        result_ids, result_names = self.src.get_anc_by_name(names)
        self.assertEqual(result_ids[0], [131567, 2157, 1783276])

    # @unittest.skip('passed')
    def test_get_anc_by_id(self):
        ids = [743725, 2107591]
        result_ids, result_names = self.src.get_anc_by_id(ids)
        self.assertEqual(result_ids[0], [131567, 2157, 1783276])
        self.assertEqual(result_ids[1], [131567, 2157, 1783276, 743725, 2107589, 2107590])

    # @unittest.skip('passed')
    def test_get_common_ancestor(self):
        names = ['Candidatus Diapherotrites', 'Candidatus Forterrea multitransposorum CG_2015-17_Forterrea_25_41']
        ancestor, distances = self.src.get_common_ancestor(names[0], names[1])
        self.assertEqual(ancestor, 1783276)
        self.assertEqual(distances[0], 1)
        self.assertEqual(distances[1], 4)

    def test_get_rank(self):
        ids = [131567, 2759, 33154, 33208, 6072, 33213, 33511, 7711, 9526, 314295, 9604, 207598, 9605, 9606]
        ranks = self.src.get_rank(ids)
        self.assertEqual(ranks[3], 'kingdom')
        self.assertEqual(ranks[1], '+')