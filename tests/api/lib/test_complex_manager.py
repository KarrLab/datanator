""" Test of text metabolite manager

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2018-08-08
:Copyright: 2018, Karr Lab
:License: MIT
"""

from kinetic_datanator.api.lib.complex.manager import complex_manager
from kinetic_datanator.core import common_schema, models
import unittest


class TestProteinComplexManager(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.signal_peptidase = complex_manager.data_source.session.query(models.ProteinComplex).filter_by(complex_name = 'signal peptidase I').first()

    def test_get_complex_by_id(self):
        prot = complex_manager.get_complex_by_id(self.signal_peptidase.id)
        self.assertEqual(prot.complex_name, self.signal_peptidase.complex_name)

    def test__search(self):
        result = complex_manager._search('signal peptidase I')
        self.assertIn(self.signal_peptidase, result)

    def test__port(self):
        ported_prot = complex_manager._port(self.signal_peptidase)
        self.assertEqual(ported_prot.name, 'signal peptidase I')
        self.assertEqual(ported_prot.cross_references, [])
