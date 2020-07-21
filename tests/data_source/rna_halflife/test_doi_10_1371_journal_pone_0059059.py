import unittest
from datanator.data_source.rna_halflife import doi_10_1371_journal_pone_0059059
import tempfile
import shutil
from datanator_query_python.config import config
import tabula


class TestProteinAggregate(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        des_db = 'test'
        src_db = 'datanator'
        cls.protein_col = 'uniprot'
        cls.rna_col = 'rna_halflife'
        cls.cache_dir = tempfile.mkdtemp()
        conf = config.TestConfig()
        username = conf.USERNAME
        password = conf.PASSWORD
        MongoDB = conf.SERVER    
        cls.src = doi_10_1371_journal_pone_0059059.Halflife(server=MongoDB, src_db=src_db,
        protein_col=cls.protein_col, authDB='admin', readPreference='nearest',
        username=username, password=password, verbose=True, max_entries=20,
        des_db=des_db, rna_col=cls.rna_col, cache_dir=cls.cache_dir)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dir)
        cls.src.uniprot_collection_manager.db_obj.drop_collection(cls.protein_col)
        cls.src.db_obj.drop_collection(cls.rna_col)
        cls.src.uniprot_collection_manager.client.close()
        cls.src.client.close()
        cls.src.uniprot_query_manager.client.close()

    @unittest.skip('Needs acadamic IPs.')
    def test_fill_rna_halflife(self):
        url = 'https://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0059059.s002'
        df = tabula.read_pdf(url, pandas_options={'header': None, 'na_values': 'ND'}, pages='all')
        df.columns = ['gene_name', 'a', 'b', 'c', 'd', 'e']
        self.src.fill_rna_half_life(df, ['Lactococcus lactis subsp. lactis Il1403', 272623])