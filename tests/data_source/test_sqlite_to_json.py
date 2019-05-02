'''Tests of sqlite_to_json
'''

import unittest
import shutil
import tempfile
from datanator.data_source import sqlite_to_json
import json


class TestSQLToJSON(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        cls.query = "select * from "
        cls.quilt_package = 'datanator'
        cls.system_path = 'SabioRk.sqlite'
        cls.src = sqlite_to_json.SQLToJSON(
            cls.query, cls.cache_dirname, cls.quilt_package, cls.system_path)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    def test_table(self):
        table_names = self.src.table()
        self.assertEqual(table_names, ['synonym', 'resource', 'entry', 'compound_structure', 'entry_synonym', 'entry_resource', 'compound', 'enzyme',
                                       'compartment', 'compound_compound_structure', 'kinetic_law', 'enzyme_subunit', 'kinetic_law_resource', 'reaction_participant', 'parameter'])

    def test_query_table(self):
        tables = ['kinetic_law', 'resource']
        data = {}
        for table in tables:
            result = self.src.query_table(table)
            data[table] = result
        null = None
        self.assertEqual(data['kinetic_law'][0],
                         {"_id": 5,
                          "enzyme_id": 3,
                          "enzyme_compartment_id": null,
                          "enzyme_type": "",
                          "tissue": null,
                          "mechanism": null,
                          "equation": null,
                          "taxon": 1467,
                          "taxon_wildtype": 1,
                          "taxon_variant": "variant DSAI (N76D/N87S/S103A/V104I)",
                          "temperature": 25.0,
                          "ph": 7.5,
                          "media": "50 mM potassium phosphate, 4 % DMSO"}
                         )
        self.assertEqual(data['kinetic_law'][1],
                         {
            "_id": 11,
            "enzyme_id": 10,
            "enzyme_compartment_id": null,
            "enzyme_type": "",
            "tissue": null,
            "mechanism": null,
            "equation": null,
            "taxon": 1467,
            "taxon_wildtype": 0,
            "taxon_variant": "S156E of subtilisin DSAI (N76D/N87S/S103A/V104I)",
            "temperature": 25.0,
            "ph": 7.5,
            "media": "50 mM potassium phosphate, 4 % DMSO"
        }
        )
        self.assertEqual(data['resource'][0],
                         {
            "_id": 1,
            "namespace": "chebi",
            "id": "CHEBI:16670"
        }
        )
