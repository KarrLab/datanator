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
        self.protein_Q9CWF2 = flk.session.query(models.ProteinSubunit).filter_by(uniprot_id = 'Q9CWF2').first()
        self.protein_p53622 = flk.session.query(models.ProteinSubunit).filter_by(uniprot_id = 'P53622').first()
        self.protein_p49418 = flk.session.query(models.ProteinSubunit).filter_by(uniprot_id = 'P49418').first()
        self.rhino_complex = flk.session.query(models.ProteinComplex).filter_by(complex_name = 'Rhino-Deadlock-Cutoff Complex').first()
        self.collagen = flk.session.query(models.ProteinComplex).filter_by(complex_name = 'Collagen type III trimer').first()

    def test_get_observable_complex(self):
        complex = self.q.get_observable_complex(self.protein_p53622)
        self.assertEqual(complex[0].specie.name,'COPI vesicle coat complex')
