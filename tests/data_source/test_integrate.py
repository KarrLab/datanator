import unittest
from datanator.data_source import integrate
import tempfile
import shutil


class TestIntegrate(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        # '../../datanator/data_source/cache'
        cls.db = 'datanator'
        cls.MongoDB = 'mongodb://mongo:27017/'
        cls.src = integrate.Integrate(
            cls.cache_dirname, cls.MongoDB, cls.db, verbose=True, max_entries=20)
        cls.client, cls.db = cls.src.con_db()
        cls.e_coli = 562  # E. Coli taxon id
        cls.s_cere = 4932  # Taxon ID for Saccharomyces cerevisiae

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.client.close()

    @unittest.skip('passed')
    def test_con_db(self):
        self.assertNotEqual(self.db, 'Server not available')

    @unittest.skip('passed')
    def test_retrieve_sabio(self):
        collection = self.db['sabio_rk']
        cursor = self.src.retrieve_sabio(collection, self.e_coli)
        self.assertEqual(cursor.count(with_limit_and_skip=True), 20)
        self.assertEqual(cursor[0]['kinlaw_id'], 47)
        self.assertEqual(cursor[0]['parameter'][0]['observed_value'], 0.074)

    #@unittest.skip('passed')
    def test_find_chebi(self):
        collection = self.db['sabio_rk']
        cursor = self.src.retrieve_sabio(collection, self.e_coli)
        substrate_kegg, product_kegg = self.src.find_sabio_kegg(cursor)
        self.assertTrue('C00019' in substrate_kegg)

    #@unittest.skip('passed')
    def test_find_metabolite_kegg(self):
        collection = self.db['ecmdb']
        substrate_kegg = ['C01092', 'C00019', 'C00019', 'C01092', 'C01092', 'C00019', 'C00019', 'C01092', 'C01137', 'C01092', 'C00019', 'C01092', 'C00019', 'C01092', 'C00001', 'C00012', 'C01092', 'C01137', 'C01092',
                          'C00019', 'C00012', 'C00001', 'C00012', 'C00001', 'C00012', 'C00001', 'C00001', 'C00012', 'C00001', 'C00012', 'C00001', 'C00012', 'C00001', 'C00012', 'C00001', 'C00012', 'C00001', 'C00012', 'C00012', 'C00001']
        product_kegg = ['C04425', 'C01037', 'C04425', 'C01037', 'C04425', 'C01037', 'C01037', 'C04425', None, 'C01037', 'C04425', 'C01037', 'C04425', 'C01037', 'C00012', 'C00012', None, 'C01037', 'C04425', 'C01037',
                        'C00012', 'C00012', 'C00012', 'C00012', 'C00012', 'C00012', 'C00012', 'C00012', 'C00012', 'C00012', 'C00012', 'C00012', 'C00012', 'C00012', 'C00012', 'C00012', 'C00012', 'C00012', 'C00012', 'C00012']
        self.assertEqual(len(substrate_kegg), len(product_kegg))
        fileds = {'accession': 1, 'm2m_id': 1, 'name': 1, 'synonyms': 1, '_id': 0}
        cursor = self.src.find_metabolite_kegg(collection, substrate_kegg, fileds)
        self.assertEqual(cursor[0]['m2m_id'], 'M2MDB000142')
        self.assertEqual(cursor[0].get('description'), None)
        self.assertEqual(cursor[1]['m2m_id'], 'M2MDB000289')
        self.assertEqual(cursor[1]['synonyms']['synonym'][4], "5'-Deoxyadenosine-5'-L-methionine disulfate ditosylate")


