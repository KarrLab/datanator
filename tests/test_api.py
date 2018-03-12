""" Tests API

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-03-12
:Copyright: 2018, Karr Lab
:License: MIT
"""

import kinetic_datanator
import types
import unittest


class ApiTestCase(unittest.TestCase):
    def test(self):
        self.assertIsInstance(kinetic_datanator, types.ModuleType)
        self.assertIsInstance(kinetic_datanator.config, types.ModuleType)
        self.assertIsInstance(kinetic_datanator.config.get_config, types.FunctionType)
        self.assertIsInstance(kinetic_datanator.core, types.ModuleType)
        self.assertIsInstance(kinetic_datanator.core.data_model, types.ModuleType)
        self.assertIsInstance(kinetic_datanator.data_query, types.ModuleType)
        self.assertIsInstance(kinetic_datanator.data_query.dna_protein_interactions, types.ModuleType)
        self.assertIsInstance(kinetic_datanator.data_source, types.ModuleType)
        self.assertIsInstance(kinetic_datanator.data_source.array_express_tools, types.ModuleType)
        self.assertIsInstance(kinetic_datanator.data_source.array_express_tools.ensembl_tools, types.ModuleType)
        self.assertIsInstance(kinetic_datanator.data_source.process_rna_seq, types.ModuleType)
        self.assertIsInstance(kinetic_datanator.data_source.process_rna_seq.get_processed_data, types.FunctionType)
        self.assertIsInstance(kinetic_datanator.data_source.process_rna_seq.download_cdna, types.ModuleType)
        self.assertIsInstance(kinetic_datanator.util, types.ModuleType)
        self.assertIsInstance(kinetic_datanator.util.molecule_util, types.ModuleType)
