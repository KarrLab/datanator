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
    def setUpClass(self):
        self.protein_p53622 = complex_manager.data_source.session.query(models.ProteinSubunit).filter_by(uniprot_id = 'P53622').first()

    def test_get_observable_complex(self):
        complex = complex_manager.get_observable_complex(self.protein_p53622)
        self.assertEqual(complex[0].specie.name,'COPI vesicle coat complex')
