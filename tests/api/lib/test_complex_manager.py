""" Test of text metabolite manager

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2018-08-08
:Copyright: 2018, Karr Lab
:License: MIT
"""

from datanator.api.lib.complex.manager import complex_manager
from datanator.core import common_schema, models
import unittest


class TestProteinComplexManager(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.signal_peptidase = complex_manager.data_source.session.query(models.ProteinComplex).filter_by(complex_name = 'signal peptidase I').first()
        cls.rhino_complex = complex_manager.data_source.session.query(models.ProteinComplex).filter_by(complex_name = 'Rhino-Deadlock-Cutoff Complex').first()



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


    def test_get_observable_subunits(self):
        observed_subunits = complex_manager.get_observable_subunits(self.rhino_complex)
        self.assertEqual(set(['Q9VIF5', 'Q9V629', 'Q7JXA8']), set([subunit.specie.uniprot_id for subunit in observed_subunits]))

    def test_get_subunits_by_known_complex(self):
        subunits = complex_manager.get_subunits_by_known_complex(self.rhino_complex.complex_name)
        self.assertEqual(set(['Q9VIF5', 'Q9V629', 'Q7JXA8']), set([subunit.uniprot_id for subunit in subunits]))

    def test_get_known_complex_by_subunit(self):
        complex = complex_manager.get_known_complex_by_subunit('Q9CWF2').all()

        self.assertEqual(len(complex), 1)
        self.assertEqual(complex[0].complex_name, 'Profilin 1 complex')
        self.assertEqual(complex[0].go_id, 'GO:0006461;GO:0006897;GO:0030036')
