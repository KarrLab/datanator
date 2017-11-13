
import unittest
from kinetic_datanator.core import data_model, common_schema
from kinetic_datanator.data_query import protein_protein_interactions as ppi
import unittest


class TestProteinProteinInteractionGenerator(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.protein_Q9CWF2 = data_model.ProteinSpecie(uniprot_id = 'Q9CWF2')
        self.protein = data_model.ProteinSpecie(uniprot_id = 'P53622')

    def get_observed_values(self):
        pass

    def test_get_observable_interactions_and_complex(self):
        q= ppi.ProteintoProteinInteractionQueryGenerator()

        interaction, complex = q.get_observable_interactions_and_complex(self.protein)

        self.assertEqual(set([i.participant_a for i in interaction]),
            set([u'uniprotkb:P41810', u'uniprotkb:P41811', u'uniprotkb:P53622',
            u'uniprotkb:P33767', u'uniprotkb:P40509']))
        self.assertEqual(complex, [])


    def test_get_interaction_by_subunit(self):
        q= ppi.ProteintoProteinInteractionQueryGenerator()

        ans = q.get_interaction_by_subunit('P49418').all()

        self.assertEqual(set([c.site_a for c in ans]),set([u'binding-associated region:545-695(MINT-1523780)',
        u'binding-associated region:846-854(MINT-8096500)', u'binding-associated region:626-695(MINT-8094602)',
        u'binding-associated region:626-695(MINT-8094834)',
        u'glutathione s tranferase tag:?-?(MINT-8096286)|binding-associated region:620-695(MINT-8096289)',
        u'glutathione s tranferase tag:?-?(MINT-8096392)|binding-associated region:620-695(MINT-8096394)',
        u'binding-associated region:626-695(MINT-8094740)', u'-',
        u'binding-associated region:620-695(MINT-8095970)|glutathione s tranferase tag:?-?(MINT-8095968)',
        u'binding-associated region:626-695(MINT-376295)', u'binding-associated region:626-695(MINT-8094630)',
        u'binding-associated region:620-695(MINT-8096476)|glutathione s tranferase tag:?-?(MINT-8096479)',
        u'binding-associated region:545-695(MINT-1523797)',
        u'binding-associated region:620-695(MINT-8095931)|glutathione s tranferase tag:?-?(MINT-8095934)',
        u'glutathione s tranferase tag:?-?(MINT-8095915)|binding-associated region:620-695(MINT-8095918)',
        u'glutathione s tranferase tag:?-?(MINT-8096369)|binding-associated region:620-695(MINT-8096366)',
        u'binding-associated region:626-695(MINT-8094654)', u'binding-associated region:626-695(MINT-8094709)']))

    def test_get_known_complex_by_subunit(self):
        q= ppi.ProteintoProteinInteractionQueryGenerator()

        ans = q.get_known_complex_by_subunit(self.protein_Q9CWF2.uniprot_id).all()

        self.assertEqual(len(ans), 1)
        self.assertEqual(ans[0].complex_name, 'Profilin 1 complex')
        self.assertEqual(ans[0].go_id, 'GO:0006461;GO:0006897;GO:0030036')

    def test_get_subunits_by_known_complex(self):
        q= ppi.ProteintoProteinInteractionQueryGenerator()

        name = 'Succinyl-CoA synthetase, GDP-forming'

        ans = q.get_subunits_by_known_complex(name).all()

        self.assertEqual([c.uniprot_id for c in ans],
            [u'O19069', u'P53590', u'P53597', u'Q96I99', u'Q9WUM5', u'Q9Z2I8'])
