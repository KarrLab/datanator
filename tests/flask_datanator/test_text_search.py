from kinetic_datanator.flask_datanator import text_search
import unittest


class TestTextSearchSession(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.sesh = text_search.TextSearchSession()

    def test_return_all_queries(self):
        """
        Test ability to collect all relevant objects and run quries from related objects
        """

        pass

    def test_collect_objects(self):
        """
        Tests ability to collect objects from full text search of database
        """
        search_dict= self.sesh.collect_objects('2-Oxopentanoate')
        for c in search_dict['Compound']:
            self.assertIn(c.compound_name,
                         [u'4-Hydroxy-2-oxopentanoate', u'(S)-4-Amino-5-oxopentanoate',
                              u'(R)-3-Hydroxy-3-methyl-2-oxopentanoate', u'(R) 2,3-Dihydroxy-3-methylvalerate',
                              u'2-Oxopentanoate', u'Ketoleucine', u'4-Methyl-2-oxopentanoate',
                              u'2-Isopropyl-3-oxosuccinate'])

        search_dict= self.sesh.collect_objects('MCM complex')
        self.assertEqual(set([c.su_cmt for c in search_dict['ProteinComplex']]),
                         set([u'None', u'Q91876(1)|P55862(1)|Q7ZY18(1)|P55861(1)|P30664(1)|P49739(1)',
                         u'P29469(1)|P53091(1)|P30665(1)|P24279(1)|P38132(1)|P29496(1)',
                         u'O75001(1)|P40377(1)|P41389(1)|P49731(1)|P29458(1)|P30666(1)',
                         u'Q9VGW6(1)|P49735(1)|Q26454(1)|Q9V461(1)|Q9XYU0(1)|Q9XYU1(1)',
                         u'P25206(1)|Q61881(1)|P49717(1)|P97310(1)|P97311(1)|P49718(1)',
                         u'P33991(1)|P33993(1)|P33992(1)|P25205(1)|Q14566(1)|P49736(1)']))

        search_dict = self.sesh.collect_objects('P49418')
        for c in search_dict['ProteinInteractions']:
            self.assertIn(c.interaction,
                        set([u'intact:EBI-7121760|mint:MINT-8094677', u'intact:EBI-7121870|mint:MINT-8094737',
                        u'intact:EBI-7121552|mint:MINT-16056', u'intact:EBI-7122020|mint:MINT-8094831',
                        u'intact:EBI-7121659|mint:MINT-8094627', u'intact:EBI-7121975|mint:MINT-8094817',
                        u'intact:EBI-7122056|mint:MINT-8094848', u'intact:EBI-7121816|mint:MINT-8094722',
                        u'intact:EBI-7121780|mint:MINT-8094706', u'intact:EBI-7121634|mint:MINT-8094596',
                        u'intact:EBI-7121926|mint:MINT-8094793', u'intact:EBI-7121710|mint:MINT-8094651',
                        u'intact:EBI-7121911|mint:MINT-8094755']))

    def test_text_return_metabolite_concentration(self):
        search_dict = self.sesh.collect_objects('Xylulose')
        concentrations = self.sesh.text_return_metabolite_concentration()
        for items in concentrations:
            self.assertIn(items[0].value, [1020.0, 686.0, 1320.0])
            break

    def test_text_return_protein_concentration(self):
        # search_dict = self.sesh.collect_objects('Q72D86')
        # print(search_dict)
        # abundances = self.sesh.text_return_protein_concentration()
        # print(abundances)
        pass

    def test_text_return_reaction_kinetics(self):
        pass

    def test_text_return_protein_protein_interactions(self):
        search_dict = self.sesh.collect_objects('O43426')
        interactions, plex = self.sesh.text_return_protein_protein_interactions()
        print(interactions, plex)

    def test_text_return_protein_dna_interactions(self):
        pass

    def test_text_return_dna_protein_interactions(self):
        pass
