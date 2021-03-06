import unittest
import shutil
import tempfile
from datanator.data_source import gene_ortholog
import datanator.config.core


class TestKeggOrgCode(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        db = 'test'
        username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
        cls.src = gene_ortholog.KeggGeneOrtholog(MongoDB, des_db=db, src_db='datanator', max_entries=10, username=username, password=password,
                                                readPreference='nearest', authSource='admin', verbose=True)
        cls.query = 'aly:ARALYDRAFT_486312'

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.src.des_db.drop_collection(cls.src.collection_str)

    @unittest.skip('passed')
    def test_parse_html(self):
        soup = self.src.get_html(self.query)
        results = self.src.parse_html(soup)
        for i, result in enumerate(results):
            if i == self.src.max_entries:
                break
            print(result)

    @unittest.skip('passed')
    def test_uniprot_to_org_gene(self):
        uniprot_id = 'Q05758'
        result = self.src.uniprot_to_org_gene(uniprot_id)
        self.assertEqual('ath:AT3G58610', result)
        uniprot_id = 'Q8N7E2'
        result = self.src.uniprot_to_org_gene(uniprot_id)
        print(result)

    def test_parse_gene_info(self):
        result = self.src.parse_gene_info('100008727')
        self.assertEqual(['AAD18037.1', 'AAD38154.1', 'AFS49951.1', 'NP_001075529.1', 'Q9XSZ4.1', 'XP_008265676.1', 'XP_017202733.1'], result)