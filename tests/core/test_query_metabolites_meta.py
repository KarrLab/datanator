import unittest
from datanator_query_python.query import query_metabolites_meta
import tempfile
import shutil
import configparser
import os


class TestQueryMetabolitesMeta(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.db = 'datanator'
        parser = configparser.ConfigParser(allow_no_value=True)
        parser.read(os.path.expanduser('~/.wc/datanator.ini'))
        username = parser.get('mongodb', 'user')
        password = parser.get('mongodb', 'password')
        MongoDB = parser.get('mongodb', 'server')
        port = int(parser.get('mongodb', 'port'))
        replSet = parser.get('mongodb', 'replSet')
        cls.src = query_metabolites_meta.QueryMetabolitesMeta(
            cache_dirname=cls.cache_dirname, MongoDB=MongoDB, replicaSet=replSet, db=cls.db,
                 verbose=True, max_entries=20, username = username, password = password)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    # @unittest.skip('passed')
    def test_get_metabolite_synonyms(self):
        c = "Isopropylglyoxylic acid"
        compounds = ["Isopropylglyoxylic acid", "3-OH-iso-but", '']
        ran = 'something'
        rxn, syn = self.src.get_metabolite_synonyms(c)
        rxn_ran, syn_ran = self.src.get_metabolite_synonyms(ran)
        rxns, syns = self.src.get_metabolite_synonyms(compounds)
        self.assertTrue("Ketovaline" in syn[c])
        self.assertTrue("3-OH-isobutyrate" in syns['3-OH-iso-but'])
        self.assertEqual(None, syns['synonyms'])
        self.assertTrue('does not exist' in syn_ran[ran])

        empty = ''
        rxn, syn = self.src.get_metabolite_synonyms(empty)
        self.assertEqual(syn, {'synonyms': None})

    # @unittest.skip('passed')
    def test_get_metabolite_inchi(self):
        compounds = ['Ketovaline']
        inchis = self.src.get_metabolite_inchi(compounds)
        self.assertEqual(inchis[0]['inchi'], 'InChI=1S/C5H8O3/c1-3(2)4(6)5(7)8/h3H,1-2H3,(H,7,8)')
        compound = ['3-Hydroxy-2-methylpropanoic acid']
        inchi = self.src.get_metabolite_inchi(compound)
        self.assertEqual(inchi[0]['m2m_id'], 'M2MDB006130')

    def test_get_ids_from_hash(self):
        hashed_inchi_1 = 'QHKABHOOEWYVLI-UHFFFAOYSA-N'
        hasehd_inchi_2 = 'YBJHBAHKTGYVGT-ZKWXMUAHSA-N'
        result_1 = self.src.get_ids_from_hash(hashed_inchi_1)
        result_2 = self.src.get_ids_from_hash(hasehd_inchi_2)
        self.assertEqual(result_1, {'m2m_id': 'M2MDB000606', 'ymdb_id': 'YMDB00365'})
        self.assertEqual(result_2, {'m2m_id': 'M2MDB000008', 'ymdb_id': 'YMDB00282'})

    # @unittest.skip('passed')
    def test_get_metabolite_name_by_hash(self):
        compounds = ['QHKABHOOEWYVLI-UHFFFAOYSA-N',
                    'YBJHBAHKTGYVGT-ZKWXMUAHSA-N']
        result = self.src.get_metabolite_name_by_hash(compounds)
        self.assertEqual(result[0], 'Ketovaline')
        self.assertEqual(result[1], 'Vitamin-h')
        compound = ['TYEYBOSBBBHJIV-UHFFFAOYSA-N']
        result = self.src.get_metabolite_name_by_hash(compound)

    # @unittest.skip('passed')
    def test_get_metabolite_hashed_inchi(self):
        compounds = ['alpha-Ketoisovaleric acid', 'delta-Biotin factor S', 'Rovimix H 2']
        hashed_inchi = self.src.get_metabolite_hashed_inchi(compounds)
        self.assertEqual(hashed_inchi[1], 'YBJHBAHKTGYVGT-ZKWXMUAHSA-N')
        self.assertEqual(hashed_inchi[1], hashed_inchi[2])
        compound = ['3-Hydroxy-2-methylpropanoic acid']
        hashed_inchi = self.src.get_metabolite_hashed_inchi(compound)
        self.assertEqual(hashed_inchi, ['DBXBTMSZEOQQDU-VKHMYHEASA-N'])

    # @unittest.skip('passed')
    def test_get_metabolite_similar_compounds(self):

        compound = ['Methyl-Pyruvic acid']
        raw1, result1 = self.src.get_metabolite_similar_compounds(compound, num = 3, threshold = 0.6)
        self.assertTrue(list(raw1[0].keys())[0], 'FERIUCNNQQJTOY-UHFFFAOYSA-N')
        raw2, result2 = self.src.get_metabolite_similar_compounds(compound, num = 10, threshold = 0.7)
        self.assertEqual(len(raw2), 10)
        raw3, result3 = self.src.get_metabolite_similar_compounds(compound, num = 10, threshold = 0.9)
        # self.assertEqual(raw3[0]['raw'], -1)
        # compound = ['β-D-Ribopyranose']
        # raw4, result4 = self.src.get_metabolite_similar_compounds(compound, num = 30, threshold = 0.6)
        # self.assertTrue('d' not in list(result4[0].keys()))
        # self.assertTrue('Auto inducer 2' in list(result4[0].keys()))

        