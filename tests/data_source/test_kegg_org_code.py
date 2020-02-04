import unittest
import shutil
import tempfile
from datanator.data_source import kegg_org_code
import datanator.config.core


class TestKeggOrgCode(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        db = 'test'
        username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
        cls.src = kegg_org_code.KeggOrgCode(MongoDB, db, max_entries=100, username=username, password=password,
                                            readPreference='nearest', authSource='admin', verbose=True, collection_str='kegg_organisms_code')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.src.db.drop_collection(cls.src.collection_str)

    @unittest.skip('passed')
    def test_parse_html_iter(self):
        results = self.src.parse_html_iter()
        for i, result in enumerate(results):
            if i == self.src.max_entries:
                break
            print(result)

    @unittest.skip('passed')
    def test_make_bulk(self):
        result = self.src.make_bulk(offset=500)
        print(result)
        self.assertEqual(len(result), 100)

    @unittest.skip('passed')
    def test_get_ncbi_id_rest(self):
        name = "homo sapiens (human)"
        self.assertEqual(self.src.get_ncbi_id_rest(name), 9606)

    def test_get_ncbi_id(self):
        name = 'Mus musculus'
        print(self.src.get_ncbi_id(name))