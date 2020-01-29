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
        cls.src = kegg_org_code.KeggOrgCode(MongoDB, db, max_entries=10, username=username, password=password,
                                            readPreference='nearest', authSource='admin', verbose=True, collection_str='kegg_organisms_code')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.src.db.drop_collection(self.src.collection_str)

    @unittest.skip('passed')
    def test_parse_ids(self):
        _ids = self.src.parse_ids()
        result = []
        i = 0
        for i, _id in enumerate(_ids):
            if i == self.src.max_entries:
                return (result)
            result.append(_id)
            i += 1
        self.assertEqual(result[0], 'hsa')
        self.assertEqual(len(result), self.src.max_entries)

    @unittest.skip('passed')
    def test_parse_names(self):
        names = self.src.parse_names()
        result = []
        i = 0
        for i, name in enumerate(names):
            if i == self.src.max_entries:
                return (result)
            result.append(name)
            i += 1
        self.assertEqual(result[0], 'Homo sapiens (human)')
        self.assertEqual(len(result), self.src.max_entries)

    @unittest.skip('passed')
    def test_make_bulk(self):
        result = self.src.make_bulk(offset=500)
        print(result)
        self.assertEqual(len(result), 100)

    @unittest.skip('passed')
    def test_bulk_load(self):
        self.src.bulk_load()

    @unittest.skip('passed')
    def test_parse_species_names(self):
        names = self.src.parse_species_name()
        result = []
        i = 0
        for i, name in enumerate(names):
            if i == self.src.max_entries:
                return (result)
            result.append(name)
            i += 1
        self.assertEqual(result[0], 'Homo sapiens')
        self.assertEqual(len(result), self.src.max_entries)

    @unittest.skip('passed')
    def test_parse_species_ids(self):
        names = self.src.parse_species_id()
        result = []
        i = 0
        for i, name in enumerate(names):
            if i == self.src.max_entries:
                return (result)
            result.append(name)
            i += 1
        self.assertEqual(result[0], '9606')
        self.assertEqual(len(result), self.src.max_entries)

    @unittest.skip('passed')
    def test_matching_species_name_id(self):
        names = self.src.parse_species_name()
        ids = self.src.parse_species_id()
        i = 0
        result = ''
        for i, (name, _id) in enumerate(zip(names, ids)):
            result = {name: _id}
        self.assertEqual(result, {'Candidatus Bathyarchaeota archaeon BA2': 1700836})

    # @unittest.skip('passed')
    def test_get_ncbi_id(self):
        name = 'Ornithobacterium rhinotracheale ORT-UMN 88'
        self.assertEqual(self.src.get_ncbi_id(name), 1401325)
        name = 'latoieruwerwe'
        self.assertEqual(self.src.get_ncbi_id(name), None)
        name = 'Pan troglodytes'
        self.assertEqual(self.src.get_ncbi_id(name), 9598)