""" Tests API

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-03-12
:Copyright: 2018, Karr Lab
:License: MIT
"""

import datanator
import types
import unittest


class ApiTestCase(unittest.TestCase):
    def test(self):
        self.assertIsInstance(datanator, types.ModuleType)
        self.assertIsInstance(datanator.config, types.ModuleType)
        self.assertIsInstance(datanator.config.get_config, types.FunctionType)
        self.assertIsInstance(datanator.core, types.ModuleType)
        self.assertIsInstance(datanator.core.data_model, types.ModuleType)
        self.assertIsInstance(datanator.api.query, types.ModuleType)
        self.assertIsInstance(datanator.api.query.dna_protein_interactions, types.ModuleType)
        self.assertIsInstance(datanator.data_source, types.ModuleType)
        self.assertIsInstance(datanator.data_source.array_express_tools, types.ModuleType)
        self.assertIsInstance(datanator.data_source.array_express_tools.ensembl_tools, types.ModuleType)
        self.assertIsInstance(datanator.data_source.process_rna_seq, types.ModuleType)
        self.assertIsInstance(datanator.data_source.process_rna_seq.get_processed_data_samples, types.FunctionType)
        self.assertIsInstance(datanator.data_source.process_rna_seq.download_cdna, types.ModuleType)
        self.assertIsInstance(datanator.util, types.ModuleType)
        self.assertIsInstance(datanator.util.molecule_util, types.ModuleType)
