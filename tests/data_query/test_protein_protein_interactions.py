
import unittest
from kinetic_datanator.core import data_model, common_schema
from kinetic_datanator.data_query import protein_protein_interactions as ppi
import unittest


class TestProteinProteinInteractionGenerator(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.protein_Q9CWF2 = data_model.ProteinSpecie(uniprot_id = 'Q9CWF2')

    def test_get_interaction_by_subunit(self):
        pass

    def test_get_known_complex_by_subunit(self):
        q= ppi.ProteintoProteinInteractionQueryGenerator()

        ans = q.get_known_complex_by_subunit(self.protein_Q9CWF2.uniprot_id).all()

        self.assertEqual(len(ans), 1)
        self.assertEqual(ans[0].complex_name, 'Profilin 1 complex')
        self.assertEqual(ans[0].go_id, 'GO:0006461;GO:0006897;GO:0030036')

    def test_get_subunits_by_known_complex(self):
        pass
