
import unittest
from kinetic_datanator.core import data_model
from kinetic_datanator.api.query import protein_protein_interactions as ppi
from kinetic_datanator.core import models, common_schema
import tempfile
import shutil
import unittest

@unittest.skip('skip')
class ProteinInteractionandComplexQuery(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        flk = common_schema.CommonSchema(cache_dirname=self.cache_dirname)

        self.protein_Q9CWF2 = flk.session.query(models.ProteinSubunit).filter_by(uniprot_id = 'Q9CWF2').first()
        self.protein_p53622 = flk.session.query(models.ProteinSubunit).filter_by(uniprot_id = 'P53622').first()
        self.protein_p49418 = flk.session.query(models.ProteinSubunit).filter_by(uniprot_id = 'P49418').first()
        self.rhino_complex = flk.session.query(models.ProteinComplex).filter_by(complex_name = 'Rhino-Deadlock-Cutoff Complex').first()
        self.collagen = flk.session.query(models.ProteinComplex).filter_by(complex_name = 'Collagen type III trimer').first()


        self.q= ppi.ProteinInteractionandComplexQuery(cache_dirname=self.cache_dirname)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)

    def test_get_observed_result_subunit(self):
        observed_interaction, observed_complex = self.q.get_observed_result(self.protein_p49418)

        self.assertEqual(set([o.interaction.specie_a.uniprot_id for o in observed_interaction if o.interaction.specie_a.__class__.__name__ == 'ProteinSpecie']),
            set(['Q05193', 'P49418', 'Q9UQ16', 'O43426', 'P17427']))
        self.assertEqual(set([o.interaction.specie_a.sequence for o in observed_interaction if o.interaction.specie_a.__class__.__name__ == 'PolymerSpecie']),
            set(['hrpvrraap', 'qrplrrprm', 'grpprnarv', 'vrparrvlw', 'vrptraada', 'hrptrskla']))

        # self.assertEqual(set([o.metadata.method.name for o in observed_interaction]), set(['peptide array', 'phage display', 'pull down', 'enzyme linked immunosorbent assay']))
        #TODO: Need to be more encompassing of all the cross_references for order issues
        # self.assertEqual(set([o.interaction.cross_references[0].namespace for o in observed_interaction]), set(['pubmed', 'paper']))

        observed_interaction, observed_complex = self.q.get_observed_result(self.protein_p53622)

        self.assertIn(set(['P33748', 'P38682', 'Q04651',]), set([o.interaction.specie_a.uniprot_id for o in observed_interaction if o.interaction.specie_a.__class__.__name__ == 'ProteinSpecie']))

        self.assertEqual(len(observed_complex), 1)
        self.assertEqual(observed_complex[0].specie.name,'COPI vesicle coat complex')


    def test_get_observed_result_complex(self):
        observed_subunits = self.q.get_observed_result(self.rhino_complex)

        self.assertEqual(set([c.specie.uniprot_id for c in observed_subunits]), set(['Q9VIF5', 'Q7JXA8', 'Q9V629']))

        observed_subunits = self.q.get_observed_result(self.collagen)

        self.assertEqual(set([c.specie.uniprot_id for c in observed_subunits]), set(['P08121', 'P02461', 'P04258']))


    def test_get_observable_complex(self):
        complex = self.q.get_observable_complex(self.protein_p53622)

        self.assertEqual(complex[0].specie.name,'COPI vesicle coat complex')

    def test_get_observable_interactions(self):
        observed_interaction = self.q.get_observable_interactions(self.protein_p49418)

        self.assertEqual(set([o.interaction.specie_a.uniprot_id for o in observed_interaction if o.interaction.specie_a.__class__.__name__ == 'ProteinSpecie']),
            set(['Q05193', 'P49418', 'Q9UQ16', 'O43426', 'P17427']))
        self.assertEqual(set([o.interaction.specie_a.sequence for o in observed_interaction if o.interaction.specie_a.__class__.__name__ == 'PolymerSpecie']),
            set(['hrpvrraap', 'qrplrrprm', 'grpprnarv', 'vrparrvlw', 'vrptraada', 'hrptrskla']))

        # self.assertIn( set(['peptide array', 'phage display', 'pull down', 'enzyme linked immunosorbent assay']), set([o.metadata.method.name for o in observed_interaction]))
        #TODO: Need to be more encompassing of all the cross_references for order issues
        # self.assertEqual(set([o.interaction.cross_references[0].namespace for o in observed_interaction]), set(['pubmed', 'paper']))

    def test_get_observable_subunits(self):
        subunits = self.q.get_observable_subunits(self.rhino_complex)
        self.assertEqual(set([c.specie.uniprot_id for c in subunits]), set(['Q9VIF5', 'Q7JXA8', 'Q9V629']))

    def test_get_interaction_by_subunit(self):
        interactions = self.q.get_interaction_by_subunit(self.protein_p49418.uniprot_id).all()
        self.assertEqual(set([i.protein_a for i in interactions]), set(['qrplrrprm', 'hrptrskla', 'grpprnarv', 'O43426', 'Q9UQ16', 'hrpvrraap', 'Q05193', 'P49418', 'vrptraada', 'vrparrvlw', 'P17427']))
        self.assertEqual(set([i.protein_b for i in interactions]), set(['O43426', 'xrprrghal', 'rrprrslpe', 'P50570', 'krpprdgtl', 'tpkrppril', 'Q05193', 'P49418', 'hrpwrrwre', 'trpyrpthp', 'P17427']))

        interactions = self.q.get_interaction_by_subunit(self.protein_p53622.uniprot_id).all()
        self.assertEqual(set([i.protein_a for i in interactions]), set(['P33767', 'P53622', 'P41811']))
        self.assertEqual(set([i.protein_b for i in interactions]), set(['P32074', 'P41811', 'P41810', 'P43621', 'P53622', 'P40509', 'P40070']))


    def test_get_known_complex_by_subunit(self):
        complex = self.q.get_known_complex_by_subunit(self.protein_Q9CWF2.uniprot_id).all()

        self.assertEqual(len(complex), 1)
        self.assertEqual(complex[0].complex_name, 'Profilin 1 complex')
        self.assertEqual(complex[0].go_id, 'GO:0006461;GO:0006897;GO:0030036')

    def test_get_subunits_by_known_complex(self):
        subunits = self.q.get_subunits_by_known_complex(self.rhino_complex.complex_name).all()

        self.assertEqual(set([s.uniprot_id for s in subunits]), set(['Q9VIF5', 'Q7JXA8', 'Q9V629']))
