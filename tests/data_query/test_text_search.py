from kinetic_datanator.data_query import text_search
import unittest


class TestTextSearchSession(unittest.TestCase):

    # @classmethod
    # def setUpClass(self):
    #     self.sesh = text_search.TextSearchSession()

    def test_compound_object_search(self):
        # search_dict = self.sesh.text_search('2-Oxopentanoate')
        # self.assertEqual(set([c.compound_name for c in search_dict['Compound']]),
        #     set([u'4-Hydroxy-2-oxopentanoate', u'(S)-4-Amino-5-oxopentanoate',
        #     u'(R)-3-Hydroxy-3-methyl-2-oxopentanoate', u'(R) 2,3-Dihydroxy-3-methylvalerate',
        #     u'2-Oxopentanoate', u'Ketoleucine', u'4-Methyl-2-oxopentanoate',
        #     u'2-Isopropyl-3-oxosuccinate']))
        pass

    def test_compound_concentration_query(self):
        pass

    def test_reaction_query(self):
        pass

    def test_protein_object_search(self):
        pass

    def test_protein_abundance_query(self):
        pass

    def test_interaction_object_search(self):
        pass

    def test_interaction_query(self):
        pass
