import unittest
from datanator.schema_2 import migrate_uniprot


class TestMU(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.src = migrate_uniprot.MigrateUniprot()

    @classmethod
    def tearDownClass(cls):
        cls.src.client.close()

    def test_get_multi_cursor(self):
        results = self.src.get_multi_cursor(4)
        self.assertEqual(len(results), 4)