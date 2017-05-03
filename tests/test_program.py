""" Tests of kinetic_datanator

:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-06
:Copyright: 2017, Karr Lab
:License: MIT
"""

from os import path
from kinetic_datanator import datanator
from kinetic_datanator import inchi_generator
from kinetic_datanator import query_string_manipulator
from kinetic_datanator import reaction_queries
from kinetic_datanator import sabio_interface
from kinetic_datanator import translator_for_sabio
from kinetic_datanator.util import reaction_util
from kinetic_datanator.util import molecule_util
from kinetic_datanator.util import warning_util
import openpyxl
import os
import sys
import unittest

warning_util.disable_warnings()


class TestProgram(unittest.TestCase):


    def test_getTaxonomicDistance(self):
        distance = TaxonFinder.getTaxonomicDistance('Mycoplasma pneumoniae', 'homo sapiens')
        self.assertTrue(distance == 8)

        distance = TaxonFinder.getTaxonomicDistance('homo sapiens', 'Mycoplasma pneumoniae')
        self.assertTrue(distance == 30)

        distance = TaxonFinder.getTaxonomicDistance('Saccharomyces cerevisiae Fleischmanns baking yeast', 'Mycoplasma pneumoniae')
        self.assertTrue(distance == 13)

        distance = TaxonFinder.getTaxonomicDistance('Escherichia coli', 'Bacillus algicola')
        self.assertTrue(distance == 6)

        #try what happens if it doesn't recognize the species
        distance = TaxonFinder.getTaxonomicDistance('mycoplasma pneumoniae', 'Streptococcus SomeLatinusNamus')
        self.assertTrue(distance == 6)

        #if it doesn't recognize the genus, it should return an empty string
        distance = TaxonFinder.getTaxonomicDistance('mycoplasma pneumoniae', 'SomeLatinusNamus speciesus')
        self.assertTrue(distance == "")



    def test_get_query_search_string(self):
        #test query_string_manipulator

        #test a regular case
        a = ['dGMP', 'GMP', 'Lactose 6-phosphate', 'dGDP', "Orotidine 5'-phosphate", "Guanosine 3'-phosphate", "2',3'-Cyclic GMP", 'L-Arogenate', 'N-Acylneuraminate 9-phosphate', "Maltose 6'-phosphate", "5-Amino-6-(5'-phosphoribitylamino)uracil", '6-Phospho-beta-D-glucosyl-(1,4)-D-glucose', '2-Amino-4-hydroxy-6-(D-erythro-1,2,3-trihydroxypropyl)-7,8- dihydropteridine', '2-Amino-4-hydroxy-6-(erythro-1,2,3-trihydroxypropyl)dihydropteridine triphosphate', 'Dihydroneopterin phosphate', 'Ganciclovir', '8-Br-cGMP', "2'-Deoxyguanosine 3'-phosphate", "8-Azaguanosine-5'-monophosphate", '8-oxo-dGMP', 'Dihydroneopterin triphosphate', '8-oxo-dGTP', "2'-Deoxy-8-hydroxyguanosine"]
        b = ['dGMP', 'GMP', 'Lactose 6-phosphate', 'dGDP', "Orotidine 5'-phosphate", "Guanosine 3'-phosphate", "2',3'-Cyclic GMP", 'L-Arogenate', 'N-Acylneuraminate 9-phosphate', "Maltose 6'-phosphate", "5-Amino-6-(5'-phosphoribitylamino)uracil", '6-Phospho-beta-D-glucosyl-(1,4)-D-glucose', '2-Amino-4-hydroxy-6-(D-erythro-1,2,3-trihydroxypropyl)-7,8- dihydropteridine', '2-Amino-4-hydroxy-6-(erythro-1,2,3-trihydroxypropyl)dihydropteridine triphosphate', 'Dihydroneopterin phosphate', 'Ganciclovir', '8-Br-cGMP', "2'-Deoxyguanosine 3'-phosphate", "8-Azaguanosine-5'-monophosphate", '8-oxo-dGMP', 'Dihydroneopterin triphosphate', '8-oxo-dGTP', "2'-Deoxy-8-hydroxyguanosine"]
        c = ["H2O"]
        d = ["phosphate"]
        participants = []
        participants.append([a, c])
        participants.append([b, d])
        response = query_string_manipulator.getQuerySearchString(participants)
        expectedAnswer = """((Substrate:"dGMP" OR Substrate:"GMP" OR Substrate:"Lactose 6-phosphate" OR Substrate:"dGDP" OR Substrate:"Orotidine 5'-phosphate" OR Substrate:"Guanosine 3'-phosphate" OR Substrate:"2',3'-Cyclic GMP" OR Substrate:"L-Arogenate" OR Substrate:"N-Acylneuraminate 9-phosphate" OR Substrate:"Maltose 6'-phosphate" OR Substrate:"5-Amino-6-(5'-phosphoribitylamino)uracil" OR Substrate:"6-Phospho-beta-D-glucosyl-(1,4)-D-glucose" OR Substrate:"2-Amino-4-hydroxy-6-(D-erythro-1,2,3-trihydroxypropyl)-7,8- dihydropteridine" OR Substrate:"2-Amino-4-hydroxy-6-(erythro-1,2,3-trihydroxypropyl)dihydropteridine triphosphate" OR Substrate:"Dihydroneopterin phosphate" OR Substrate:"Ganciclovir" OR Substrate:"8-Br-cGMP" OR Substrate:"2'-Deoxyguanosine 3'-phosphate" OR Substrate:"8-Azaguanosine-5'-monophosphate" OR Substrate:"8-oxo-dGMP" OR Substrate:"Dihydroneopterin triphosphate" OR Substrate:"8-oxo-dGTP" OR Substrate:"2'-Deoxy-8-hydroxyguanosine") AND (Substrate:"H2O")) AND ((Product:"dGMP" OR Product:"GMP" OR Product:"Lactose 6-phosphate" OR Product:"dGDP" OR Product:"Orotidine 5'-phosphate" OR Product:"Guanosine 3'-phosphate" OR Product:"2',3'-Cyclic GMP" OR Product:"L-Arogenate" OR Product:"N-Acylneuraminate 9-phosphate" OR Product:"Maltose 6'-phosphate" OR Product:"5-Amino-6-(5'-phosphoribitylamino)uracil" OR Product:"6-Phospho-beta-D-glucosyl-(1,4)-D-glucose" OR Product:"2-Amino-4-hydroxy-6-(D-erythro-1,2,3-trihydroxypropyl)-7,8- dihydropteridine" OR Product:"2-Amino-4-hydroxy-6-(erythro-1,2,3-trihydroxypropyl)dihydropteridine triphosphate" OR Product:"Dihydroneopterin phosphate" OR Product:"Ganciclovir" OR Product:"8-Br-cGMP" OR Product:"2'-Deoxyguanosine 3'-phosphate" OR Product:"8-Azaguanosine-5'-monophosphate" OR Product:"8-oxo-dGMP" OR Product:"Dihydroneopterin triphosphate" OR Product:"8-oxo-dGTP" OR Product:"2'-Deoxy-8-hydroxyguanosine") AND (Product:"phosphate"))"""
        self.assertEqual(response, expectedAnswer)

        #test a blank:
        blank = [[],[]]
        response = query_string_manipulator.getQuerySearchString(blank)
        expectedAnswer = ""
        self.assertEqual(response, "")

        a = []
        b = ['dGMP', 'GMP', 'Lactose 6-phosphate', 'dGDP', "Orotidine 5'-phosphate", "Guanosine 3'-phosphate", "2',3'-Cyclic GMP", 'L-Arogenate', 'N-Acylneuraminate 9-phosphate', "Maltose 6'-phosphate", "5-Amino-6-(5'-phosphoribitylamino)uracil", '6-Phospho-beta-D-glucosyl-(1,4)-D-glucose', '2-Amino-4-hydroxy-6-(D-erythro-1,2,3-trihydroxypropyl)-7,8- dihydropteridine', '2-Amino-4-hydroxy-6-(erythro-1,2,3-trihydroxypropyl)dihydropteridine triphosphate', 'Dihydroneopterin phosphate', 'Ganciclovir', '8-Br-cGMP', "2'-Deoxyguanosine 3'-phosphate", "8-Azaguanosine-5'-monophosphate", '8-oxo-dGMP', 'Dihydroneopterin triphosphate', '8-oxo-dGTP', "2'-Deoxy-8-hydroxyguanosine"]
        c = []
        d = ["phosphate"]
        participants = []
        participants.append([a, c])
        participants.append([b, d])
        response = query_string_manipulator.getQuerySearchString(participants)
        expectedAnswer = """((Product:"dGMP" OR Product:"GMP" OR Product:"Lactose 6-phosphate" OR Product:"dGDP" OR Product:"Orotidine 5'-phosphate" OR Product:"Guanosine 3'-phosphate" OR Product:"2',3'-Cyclic GMP" OR Product:"L-Arogenate" OR Product:"N-Acylneuraminate 9-phosphate" OR Product:"Maltose 6'-phosphate" OR Product:"5-Amino-6-(5'-phosphoribitylamino)uracil" OR Product:"6-Phospho-beta-D-glucosyl-(1,4)-D-glucose" OR Product:"2-Amino-4-hydroxy-6-(D-erythro-1,2,3-trihydroxypropyl)-7,8- dihydropteridine" OR Product:"2-Amino-4-hydroxy-6-(erythro-1,2,3-trihydroxypropyl)dihydropteridine triphosphate" OR Product:"Dihydroneopterin phosphate" OR Product:"Ganciclovir" OR Product:"8-Br-cGMP" OR Product:"2'-Deoxyguanosine 3'-phosphate" OR Product:"8-Azaguanosine-5'-monophosphate" OR Product:"8-oxo-dGMP" OR Product:"Dihydroneopterin triphosphate" OR Product:"8-oxo-dGTP" OR Product:"2'-Deoxy-8-hydroxyguanosine") AND (Product:"phosphate"))"""
        self.assertEqual(response, expectedAnswer)

    def test_get_parsed_reaction(self):
        reaction_string = "[c]: Complex_Acp + ATP + HDCA ==> PPI + AMP + Complex_Acp_hdc"
        response = query_string_manipulator.getParsedReaction(reaction_string)
        expectedAnswer = [['Complex_Acp', 'ATP', 'HDCA'], ['PPI', 'AMP', 'Complex_Acp_hdc']]
        self.assertEqual(response, expectedAnswer)

        reaction_string = "A + B + C ==> D + E + F"
        response = query_string_manipulator.getParsedReaction(reaction_string)
        expectedAnswer = [['A', 'B', 'C'], ['D', 'E', 'F']]
        self.assertEqual(response, expectedAnswer)

        #make sure it removes hyrdrogens
        reaction_string = "[c]: Complex_Acp + H2O ==> ProtMon_MPN406 + H + Pantetheine4Phosphate"
        response = query_string_manipulator.getParsedReaction(reaction_string)
        expectedAnswer = [['Complex_Acp', 'H2O'], ['ProtMon_MPN406', 'Pantetheine4Phosphate']]
        self.assertEqual(response, expectedAnswer)

        reaction_string = "cyclobutane_dTdC[c] ==> cyclobutane_dTdC[e]"
        response = query_string_manipulator.getParsedReaction(reaction_string)
        expectedAnswer = [['cyclobutane_dTdC'], ['cyclobutane_dTdC']]
        self.assertEqual(response, expectedAnswer)

        reaction_string = "UDPG[c] + DAG161[m] ==> UDP[c] + H[c] + GlcDAG161[m]"
        response = query_string_manipulator.getParsedReaction(reaction_string)
        expectedAnswer = [['UDPG', 'DAG161'], ['UDP', 'GlcDAG161']]
        self.assertEqual(response, expectedAnswer)

        reaction_string = 'h_m + 1a25dhvitd2_m + o2_m + nadph_m --> h2o_m + nadp_m + 1a2425thvitd2_m'
        response = query_string_manipulator.getParsedReaction(reaction_string)
        expectedAnswer = [['1a25dhvitd2_m', 'o2_m', 'nadph_m'], ['h2o_m', 'nadp_m', '1a2425thvitd2_m']]
        self.assertEqual(response, expectedAnswer)

        reaction_string = 'h_m + 1a25dhvitd2_m + o2_m + nadph_m <-> h2o_m + nadp_m + 1a2425thvitd2_m'
        response = query_string_manipulator.getParsedReaction(reaction_string)
        expectedAnswer = [['1a25dhvitd2_m', 'o2_m', 'nadph_m'], ['h2o_m', 'nadp_m', '1a2425thvitd2_m']]
        self.assertEqual(response, expectedAnswer)

        reaction_string = 'h_m + 1a25dhvitd2_m + o2_m + nadph_m <=> h2o_m + nadp_m + 1a2425thvitd2_m'
        response = query_string_manipulator.getParsedReaction(reaction_string)
        expectedAnswer = [['1a25dhvitd2_m', 'o2_m', 'nadph_m'], ['h2o_m', 'nadp_m', '1a2425thvitd2_m']]
        self.assertEqual(response, expectedAnswer)

    def test_generate_generic_inchi(self):
        #test inchi_generator

        #this is a dictionary of smiles strings to corresponding inchi:
        smilesAndInchi = {
        "C1=CN2C=NC3=C(C2=N1)N=CN3[C@@H](C=O)O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C=O":"""InChI=1S/C12H14N5O13P3/c18-3-8(5-27-32(23,24)30-33(25,26)29-31(20,21)22)28-9(4-19)17-7-14-10-11-13-1-2-16(11)6-15-12(10)17/h1-4,6-9H,5H2,(H,23,24)(H,25,26)(H2,20,21,22)/t8-,9+/m0/s1""",
        "CC=C/C=N/NC(=O)C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F":"""InChI=1S/C11H7F13N2O/c1-2-3-4-25-26-5(27)6(12,13)7(14,15)8(16,17)9(18,19)10(20,21)11(22,23)24/h2-4H,1H3,(H,26,27)/b3-2?,25-4+""",
        "C1=CC(=C(C=C1Cl)[N+](=O)[O-])SC2=NN=C(S2)N":"""InChI=1S/C8H5ClN4O2S2/c9-4-1-2-6(5(3-4)13(14)15)16-8-12-11-7(10)17-8/h1-3H,(H2,10,11)""",
        "O.O.O.O.O.O.O.O.O.O.O.O.[Fe].[Fe].[Fe]": "InChI=1S/3Fe.12H2O/h;;;12*1H2"
        }

        testWorked = True
        for smiles in smilesAndInchi:
            if inchi_generator.generateGenericInchi(smiles) != inchi_generator.trimInchi(smilesAndInchi[smiles]):
                testWorked = False
        self.assertTrue(testWorked)

        #now make sure if the strings are different, it wont recognize it. 
        #some of the sttrings have gibberish in them
        smilesAndInchi = {
        "greenC1=CN2C=NC3=C(C2=N1)N=CN3[C@@H](C=O)O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C=O":"""InChI=1S/C12H14N5O13P3/c18-3-8(5-27-32(23,24)30-33(25,26)29-31(20,21)22)28-9(4-19)17-7-14-10-11-13-1-2-16(11)6-15-12(10)17/h1-4,6-9H,5H2,(H,23,24)(H,25,26)(H2,20,21,22)/t8-,9+/m0/s1""",
        "blueCC=C/C=N/NC(=O)C(C(C(C(C(C(F)(F)F)IlluminatiWasHere(F)F)(F)F)(F)F)(F)F)(F)F":"""InChI=1S/C11H7F13N2O/c1-2-3-4-25-26-5(27)6(12,13)7(14,15)8(16,17)9(18,19)10(20,21)11(22,23)24/h2-4H,1H3,(H,26,27)/b3-2?,25-4+""",
        "yelloC1=CC(=C(C=C1Cl)[N+](=O)[O-])SC2=NN=C(S2)N":"""InChI=1S/C8H5ClN4O2S2/c9-4-1-2-6(5(3-4)13(14)15)16-8-12-11-7(10)17-8/h1-3H,(H2,10,11)""",
        "O.O.O.O.O.O.O.O.O.O.O.O.O.[Fe].[Fe].[Fe]": "InChI=1S/3Fe.12H2O/h;;;12*1H2"
        }
        testWorked = False
        for smiles in smilesAndInchi:
            if inchi_generator.generateGenericInchi(smiles) == inchi_generator.trimInchi(smilesAndInchi[smiles]):
                testWorked = True
        self.assertFalse(testWorked)

        #now we need to make sure that it can also take in inchi strings, and return "generic" strings
        inchiToGenericInchi = {
        "InChI=1/C23H19NO5/c1-12(25)7-17-21-14(5-6-18-23(21)29-11-26-18)15-4-3-13-8-19-20(28-10-27-19)9-16(13)22(15)24(17)2/h3-6,8-9,17H,7,10-11H2,1-2H3":"InChI=1/C23H19NO5/c1-12(25)7-17-21-14(5-6-18-23(21)29-11-26-18)15-4-3-13-8-19-20(28-10-27-19)9-16(13)22(15)24(17)2",
        "InChI=1/C20H32O6/c1-5-15-11-13(2)16(21)10-8-6-7-9-14(3)26-18(23)12-17(22)20(25-4)19(15)24/h6-8,10,13-15,17,19-20,22,24H,5,9,11-12H2,1-4H3/b7-6+,10-8+/t13-,14-,15+,17-,19+,20+/m1/s1":"InChI=1/C20H32O6/c1-5-15-11-13(2)16(21)10-8-6-7-9-14(3)26-18(23)12-17(22)20(25-4)19(15)24",
        "InChI=1S/3CN.Fe/c3*1-2":"InChI=1S/3CN.Fe/c3*1-2",
        "InChI=1/C9H8O4.C2H5NO2.CH2O3.Al.Mg.2H2O/c1-6(10)13-8-5-3-2-4-7(8)9(11)12;3-1-2(4)5;2-1(3)4;;;;/h2-5H,1H3,(H,11,12);1,3H2,(H,4,5);(H2,2,3,4);;;2*1H2/q;;;+3;+2;;/p-5":"InChI=1/C9H8O4.C2H5NO2.CH2O3.Al.Mg.2H2O/c1-6(10)13-8-5-3-2-4-7(8)9(11)12;3-1-2(4)5;2-1(3)4;;;;"
        }
        testWorked = True
        for inchi in inchiToGenericInchi:
            if inchi_generator.generateGenericInchi(inchi) != inchi_generator.trimInchi(inchiToGenericInchi[inchi]):
                testWorked = False
        self.assertTrue(testWorked)
   	
    def test_format_ec_number_for_sabio(self):
        #this simply takes an EC number as string input, and makes a generic search string
        #by returing a string that searches for N.N.N.01, N.N.N.02, etc

        # four digit EC numbers
        self.assertEqual(translator_for_sabio.format_ec_number_for_sabio('3.1.3.4'), 'ECNumber: "3.1.3.4"')

        # three digit EC numbers should be expanded into a family of four digit EC numbers
        expectedString = 'ECNumber: ("3.1.3.1" OR "3.1.3.2" OR "3.1.3.3" OR "3.1.3.4" OR "3.1.3.5")'
        self.assertEqual(translator_for_sabio.format_ec_number_for_sabio("3.1.3", number_four_digit_ec_numbers=5), expectedString)

        #if an empty string is sent, it should return an empty string
        self.assertEqual(translator_for_sabio.format_ec_number_for_sabio(""), "")


    #testingtranslator_for_sabio

    def test_get_substrate_product_query_string(self):
        id = "Example Reaction 1"
        #[c]: ATP + Pantetheine4Phosphate ==> DPCOA + PPI
        reaction = reaction_queries.ReactionQuery(id)
        reaction.substrates = [reaction_queries.Compound("ATP", sabioNames = ['ATP']), reaction_queries.Compound("Pantetheine4Phosphate", sabioNames = ["4'-Phosphopantetheine"])]
        reaction.products = [reaction_queries.Compound("DPCOA", sabioNames = ['Dephospho-CoA', "3'-Dephospho-CoA"]), reaction_queries.Compound("PPI", sabioNames = ['Diphosphate'])]
        searchString = translator_for_sabio.getSubstrateProductQueryString(reaction)
        expectedString = """((Substrate:"ATP") AND (Substrate:"4'-Phosphopantetheine")) AND ((Product:"Dephospho-CoA" OR Product:"3'-Dephospho-CoA") AND (Product:"Diphosphate"))"""
        self.assertEqual(searchString, expectedString)

        #if sabio recognizes all of the participants except for one, it should still return a string
        reaction = reaction_queries.ReactionQuery(id)
        reaction.substrates = [reaction_queries.Compound("ATP", sabioNames = ['ATP']), reaction_queries.Compound("Pantetheine4Phosphate", sabioNames = [])]
        reaction.products = [reaction_queries.Compound("DPCOA", sabioNames = ['Dephospho-CoA', "3'-Dephospho-CoA"]), reaction_queries.Compound("PPI", sabioNames = ['Diphosphate'])]
        searchString = translator_for_sabio.getSubstrateProductQueryString(reaction)
        expectedString = """((Substrate:"ATP")) AND ((Product:"Dephospho-CoA" OR Product:"3'-Dephospho-CoA") AND (Product:"Diphosphate"))"""
        self.assertEqual(searchString, expectedString)

        #if sabio does not recognize two or more, the code should return an empty string
        reaction = reaction_queries.ReactionQuery(id)
        reaction.substrates = [reaction_queries.Compound("ATP", sabioNames = ['ATP']), reaction_queries.Compound("Pantetheine4Phosphate", sabioNames = [])]
        reaction.products = [reaction_queries.Compound("DPCOA", sabioNames = []), reaction_queries.Compound("PPI", sabioNames = ['Diphosphate'])]
        searchString = translator_for_sabio.getSubstrateProductQueryString(reaction)
        expectedString = ""
        self.assertEqual(searchString, expectedString)


        id = "Example Reaction 2"
        #[c]: PAP + H2O ==> AMP + PI
        reaction = reaction_queries.ReactionQuery(id)
        reaction.substrates = [reaction_queries.Compound("PAP", sabioNames = ["Adenosine 3',5'-bisphosphate"]), reaction_queries.Compound("H2O", sabioNames = ['H2O', 'OH-'])]
        reaction.products = [reaction_queries.Compound("AMP", sabioNames = ['AMP', "Adenine-9-beta-D-arabinofuranoside 5'-monophosphate"]), reaction_queries.Compound("PI", sabioNames = ['Dihydrogen phosphate', 'Phosphate'])]
        searchString = translator_for_sabio.getSubstrateProductQueryString(reaction)
        expectedString ="""((Substrate:"Adenosine 3',5'-bisphosphate") AND (Substrate:"H2O" OR Substrate:"OH-")) AND ((Product:"AMP" OR Product:"Adenine-9-beta-D-arabinofuranoside 5'-monophosphate") AND (Product:"Dihydrogen phosphate" OR Product:"Phosphate"))"""
        self.assertEqual(searchString, expectedString)



    #####################################
    #Tesst the reaction_queries module
    def test_reaction_queries(self):

        inputFileName = path.join(path.dirname(__file__), "fixtures", "five_reactions.xlsx")
        #turn Excel sheet into openpyxl workbook
        if not path.isdir(path.join(path.dirname(__file__), "output")):
            os.makedirs(path.join(path.dirname(__file__), "output"))
        outputFilename = path.join(path.dirname(__file__), "output", "five_reactions.xlsx")
        species = 'mycoplasma pneumoniae'


    def test_generate_reaction_queries(self):
        inputFileName = path.join(path.dirname(__file__), "fixtures", "five_reactions.xlsx")
        #turn Excel sheet into openpyxl workbook
        if not path.isdir(path.join(path.dirname(__file__), "output")):
            os.makedirs(path.join(path.dirname(__file__), "output"))
        outputFilename = path.join(path.dirname(__file__), "output", "five_reactions.xlsx")
        species = 'mycoplasma pneumoniae'
        wb = openpyxl.load_workbook(filename=inputFileName)
        #generate_reaction_queries is given an openpyxl workbook as an arg. It outputs a list
        
        rxn_queries = reaction_queries.generate_reaction_queries(wb)

        self.assertEqual(set([rxn.id for rxn in rxn_queries]), set([
            'ump_kinase',
            'gmp_kinase',
            'oligonucleotidase_dcmp_dtmp',
            'fmn_reductase',
            'nucleotidase_7gmp',
            ]))
        for rxn in rxn_queries:
            if rxn.id == 'ump_kinase':
                self.assertEqual(rxn.ec_number, ""),
                self.assertEqual(rxn.keggID, ""),
                self.assertEqual(rxn.num_participants, [2,2]),
                self.assertEqual(rxn.reaction_string, '[c]: ATP + UMP <==> UDP + ADP')
                self.assertEqual(set([comp.id for comp in rxn.substrates]), set([
                    'ATP',
                    'UMP'
                    ]))
                self.assertEqual(set([comp.id for comp in rxn.products]), set([
                    'ADP',
                    'UDP'
                    ]))

                for comp in rxn.substrates:
                    if comp.id == 'ATP':
                        self.assertEqual(comp.inchi_smiles, 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1')
                        self.assertEqual(comp.sabioNames, ['ATP'])
                    if comp.id == 'UMP':
                        self.assertEqual(comp.inchi_smiles, 'OC1C(O)C(OC1COP([O-])([O-])=O)N1C=CC(=O)NC1=O')
                        #UMP has two generic inchi structure in the Sabio Database that match the structural information we  provided.
                        #Therefore UMP has two sabio names
                        self.assertEqual(comp.sabioNames, ['UMP', "Uridine 5'-phosphate"])
                for comp in rxn.products:
                    if comp.id == 'ADP':
                        self.assertEqual(comp.inchi_smiles, 'C1=NC2=C(C(=N1)N)N=CN2C3C(C(C(O3)COP(=O)(O)OP(=O)(O)O)O)O')
                        self.assertEqual(comp.sabioNames, ['ADP'])
                    if comp.id == 'UDP':
                        self.assertEqual(comp.inchi_smiles, 'OC1C(O)C(OC1COP([O-])(=O)OP([O-])([O-])=O)N1C=CC(=O)NC1=O')
                        #UMP has two generic inchi structure in the Sabio Database that match the structural information we  provided.
                        #Therefore UMP has two sabio names
                        self.assertEqual(comp.sabioNames, ['UDP'])


    def test_generateCompounds(self):
        pass






        
    







    def test_get_sabio_data(self):
        #test sabio_interface

        #[c]: PAP + H2O ==> AMP + PI
        base_species = 'mycoplasma pneumoniae'
        searchString = """((Substrate:"Adenosine 3',5'-bisphosphate") AND (Substrate:"H2O" OR Substrate:"OH-")) AND ((Product:"AMP" OR Product:"Adenine-9-beta-D-arabinofuranoside 5'-monophosphate") AND (Product:"Dihydrogen phosphate" OR Product:"Phosphate"))"""
        results =  sabio_interface.get_sabio_data(searchString, base_species)
        self.assertEqual(len(results.entry_list), 3)

        #one of the entries should have a km of 2.0E-6 M. 
        containsKm = False
        for entry in results.entry_list:
            if entry.km == "2.0E-6":
                containsKm = True
        self.assertTrue(containsKm)

        containsVmax = False
        for entry in results.entry_list:
            if entry.vmax == "1.33333333E-5":
                containsVmax = True
        self.assertTrue(containsVmax)


        base_species = 'mycoplasma pneumoniae'
        searchString = """Product:ADP AND Substrate:AMP AND ADP"""
        results =  sabio_interface.get_sabio_data(searchString, base_species)
        self.assertEqual(len(results.entry_list), 77)


    def test_get_EC_number(self):

        def get_EC_from_structures(sub_array, prod_array):
            """ Return EC prediction from list of structures in reaction
            Args:
                sub_array(:obj:`list` of 'str'): list of structural strings (inchi or smiles) of the substrates in a reaction
                sub_array(:obj:`list` of 'str'): list of structural strings (inchi or smiles) of the products in a reaction

            Returns:
                :'str': highest scored ec number
            """

            participants = []
            for part in sub_array:
                participants.append(reaction_util.ReactionParticipant(molecule_util.Molecule(part), coefficient=-1))
            for part in prod_array:
                participants.append(reaction_util.ReactionParticipant(molecule_util.Molecule(part), coefficient=1))
            reaction = reaction_util.Reaction(participants)
            results = reaction_util.Ezyme().run(reaction)
            if len(results)>0:
            #    for thing in results:
            #        print "{}:    {}".format(thing.ec_number, thing.score)
                return reaction_util.Ezyme().run(reaction)[0].ec_number
            else:
                return "No EC Found"


        #[c]: T3P1 + S7P <==> X5P + R5P
        sub_array = ['OC(COP([O-])([O-])=O)C=O', 'OCC(=O)C(O)C(O)C(O)C(O)COP([O-])([O-])=O']
        prod_array = ['OCC(=O)C(O)C(O)COP([O-])([O-])=O', 'OC(COP([O-])([O-])=O)C(O)C(O)C=O']
        self.assertEqual(get_EC_from_structures(sub_array, prod_array), '4.1.2')

        #[c]: PI + INS <==> R1P + HYXN
        sub_array = ['OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)N1C=NC2=C(O)N=CN=C12', 'OP([O-])([O-])=O']
        prod_array = ['OCC1OC(OP([O-])([O-])=O)C(O)C1O', 'O=C1NC=NC2=C1NC=N2']
        self.assertEqual(get_EC_from_structures(sub_array, prod_array), "2.4.2")

        #[c]: H2O + ho5hydantoin_dRiboseMP ==> ho5hydantoin_dRibose + PI
        sub_array = ['O', 'OC1CC(OC1COP([O-])([O-])=O)N1C(O)C(=O)NC1=O']
        prod_array = ['OCC1OC(CC1O)N1C(O)C(=O)NC1=O', 'OP([O-])([O-])=O']
        self.assertEqual(get_EC_from_structures(sub_array, prod_array), "3.1.4")


        #[c]: e1dGMP + H2O ==> e1dG + PI
        sub_array = ['CCN1C(N)=NC2=C(N=CN2C2CC(O)C(COP([O-])([O-])=O)O2)C1=O', 'O']
        prod_array = ['CCN1C(N)=NC2=C(N=CN2C2CC(O)C(CO)O2)C1=O', 'OP([O-])([O-])=O']
        self.assertEqual(get_EC_from_structures(sub_array, prod_array), "3.1.3")

        #The following is an artifact of the previous way of searching, where the first match was recorded. Therefore 3.1.4 was possible
        sub_array = ['O', 'CCN1C(N)=NC2=C(N=CN2C2CC(O)C(COP([O-])([O-])=O)O2)C1=O']
        prod_array = ['CCN1C(N)=NC2=C(N=CN2C2CC(O)C(CO)O2)C1=O', 'OP([O-])([O-])=O']
        #self.assertEqual(get_EC_from_structures(sub_array, prod_array), "3.1.4")

        #[c]: e3dCMP + H2O ==> e3dC + PI
        #The old method was not able to find this reaction. However, the new method works!
        sub_array = ['CCNC1=NC(=O)N(C=C1)C1CC(O)C(COP([O-])([O-])=O)O1', 'O']
        prod_array = ['CCNC1=NC(=O)N(C=C1)C1CC(O)C(CO)O1', 'OP([O-])([O-])=O']
        self.assertEqual(get_EC_from_structures(sub_array, prod_array), "3.1.3")


        #try it with three
        #this fails because the normalize method in the reaction class requires a name present in each molecule. 
        #this should be changed because names are only keyword args in the molecule objects
        #IleIle[e] + ATP[c] + H2O[c] ==> IleIle[c] + PI[c] + H[c] + ADP[c]
        sub_array = ['CCC(C)C([NH3+])C(=O)NC(C(C)CC)C([O-])=O', 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1', 'O']
        prod_array = ['CCC(C)C([NH3+])C(=O)NC(C(C)CC)C([O-])=O', 'OP([O-])([O-])=O', 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1']
        self.assertEqual(get_EC_from_structures(sub_array, prod_array), "3.6.1")
        
        #things to test further:
        #1) cases where the reaction has 3 substrates and 2 products, or vice versa. I don't think kegg supports that



    
    def test_datanator(self):
        #Find Excel sheet with reaction data
        inputFileName = path.join(path.dirname(__file__), "fixtures", "five_reactions.xlsx")
        #turn Excel sheet into openpyxl workbook
        if not path.isdir(path.join(path.dirname(__file__), "output")):
            os.makedirs(path.join(path.dirname(__file__), "output"))
        outputFilename = path.join(path.dirname(__file__), "output", "five_reactions.xlsx")
        species = 'mycoplasma pneumoniae'

        #formatted data list 
        #The formatted data list is a list of FormattedData objects
        #FormattedData objects store and organize the final relevant expiremental data
        #FormattedData is the central Class that should contain all the data needed to respond to a query
        #Ideally, a FormattedData objective should be comprehensive enough to pass to a user, and the user should be able 
        #to derive all the information he wants from that object

        rxns = datanator.get_kinetic_data(inputFileName, outputFilename, species)
        self.assertEqual(set([rxn.id for rxn in rxns]), set([
            'ump_kinase',
            'gmp_kinase',
            'oligonucleotidase_dcmp_dtmp',
            'fmn_reductase',
            'nucleotidase_7gmp',
            ]))

        ##################################################
        ##################################################
        """ UMP kinase """
        rxn = next((rxn for rxn in rxns if rxn.id == 'ump_kinase'), None)
        self.assertNotEqual(rxn, None)

        #reaction_ids is there to classify how many distinct sabio-RK Ids reference the desired reaction
        #there are two reasons a reaction ID list may contain more than a single entry. Either if there is an error
        #or if there are two reactions that differ only in stereochemistry (the inchi strings and smiles cannot detect
        #stereochemical distinctions).
        #the reaction IDs for this reaction should just be 201
        self.assertEqual(rxn.reaction_ids, ['201'])
        #both km and vmax should have KineticInfo objects in their fields
        self.assertNotEqual(rxn.km_data, None)
        self.assertNotEqual(rxn.vmax_data, None)

        #Now the KineticInfo fields will be tested
        #First is some summaries of the expiremental data. 
        km_data = rxn.km_data
        self.assertEqual(km_data.ec_numbers, ['2.7.4.14', '2.7.4.22']) #it has two EC numbers for the same reaction. EC numbers are not a perfect system
        self.assertEqual(km_data.sabio_reaction_ids, ['201'])
        self.assertEqual(km_data.lift_info, "Lift Not Used")
        self.assertEqual(km_data.reaction_list, ['ATP + UMP = UDP + ADP'])

        #Each kinetic info class contains Entry objects from the sabio_interface module. 
        #The following tests will tests ensure the KineticInfo fields contains the right Entry objects
        
        closest_entry_ids = km_data.closest_entry_ids
        self.assertTrue(len(closest_entry_ids)>12 and len(closest_entry_ids)<25)
        closest_entries = km_data.closest_entries #this is a list of the most closely related expiremental entries
        self.assertTrue(len(closest_entries)>12 and len(closest_entries) < 25) #when I made this test, there were 13 entries. If there are more than 25 entries, something probably went wrong 

        #from the set of closest_entries, the entry with median km value is the median_entry
        #the median_entry is an Entry object defined in sabio_interface
        median_entry = km_data.median_entry
        self.assertEqual(median_entry.ec_number, '2.7.4.22')
        self.assertEqual(median_entry.entry_id, '17927') #this is the specific Entry Number sabio assigns to each expiremental entry
        self.assertEqual(median_entry.km, '1.0E-4') #this is the km value for this expiremental entry
        self.assertEqual(median_entry.vmax, '0.00665') #this is the vmax for this expiremental entry
        self.assertEqual(median_entry.num_participants, [2,2]) #this is a list to record the number of substrates and products
        self.assertEqual(median_entry.species, 'Streptococcus pneumoniae') #this is the species the expirement was done in
        self.assertEqual(median_entry.proximity, 6) #this is closely related the expiremental species is to the modeler's species
        self.assertEqual(median_entry.reaction_id, '201') #this is the Sabio assigned ID for this reaction in general
        
        #check the entry_id of the min and the max within km_data
        self.assertEqual(km_data.min_entry.entry_id, '51904')
        self.assertEqual(km_data.max_entry.entry_id, '51903')


        #now the vmax infor fields will  be tested
        vmaxData = rxn.vmax_data
        self.assertEqual(vmaxData.ec_numbers, ['2.7.4.14', '2.7.4.22']) #it has two EC numbers for the same reaction. EC numbers are not a perfect system
        self.assertEqual(vmaxData.sabio_reaction_ids, ['201'])
        self.assertEqual(vmaxData.lift_info, "Lift Not Used")
        self.assertEqual(vmaxData.reaction_list, ['ATP + UMP = UDP + ADP'])

        #Each kinetic info class contains Entry objects from the sabio_interface module. 
        #The following tests will tests ensure the KineticInfo fields contains the right Entry objects
        
        closest_entry_ids = vmaxData.closest_entry_ids
        self.assertTrue(len(closest_entry_ids)>44 and len(closest_entry_ids)<50)
        closest_entries = vmaxData.closest_entries #this is a list of the most closely related expiremental entries
        self.assertTrue(len(closest_entries)>44 and len(closest_entries) < 50) #when I made this test, there were 45 entries. If there are more than 50 entries, something probably went wrong 

        #from the set of closest_entries, the entry with median km value is the median_entry
        #the median_entry is an Entry object defined in sabio_interface
        median_entry = vmaxData.median_entry
        self.assertEqual(median_entry.ec_number, '2.7.4.22')
        self.assertEqual(median_entry.entry_id, '17907') #this is the specific Entry Number sabio assigns to each expiremental entry
        self.assertEqual(median_entry.km, '') #this is the km value for this expiremental entry
        self.assertEqual(median_entry.vmax, '0.00456666667') #this is the vmax for this expiremental entry
        self.assertEqual(median_entry.num_participants, [2,2]) #this is a list to record the number of substrates and products
        self.assertEqual(median_entry.species, 'Streptococcus pneumoniae') #this is the species the expirement was done in
        self.assertEqual(median_entry.proximity, 6) #this is closely related the expiremental species is to the modeler's species
        self.assertEqual(median_entry.reaction_id, '201') #this is the Sabio assigned ID for this reaction in general
        
        #check the entry_id of the min and the max within km_data
        self.assertEqual(vmaxData.min_entry.entry_id, '40434')
        self.assertEqual(vmaxData.max_entry.entry_id, '17934')

        ##################################################
        ##################################################
        rxn = next((rxn for rxn in rxns if rxn.id == 'oligonucleotidase_dcmp_dtmp'), None)
        self.assertNotEqual(rxn, None)        

        #the reaction IDs for this reaction should just be empty becuase no exact mathches were found in the Sabio Database
        self.assertEqual(rxn.reaction_ids, [])
        #both km and vmax should have KineticInfo objects in their fields
        self.assertNotEqual(rxn.km_data, None)
        self.assertNotEqual(rxn.vmax_data, None)

        #Now the KineticInfo fields will be tested
        #First is some summaries of the expiremental data. 
        km_data = rxn.km_data
        #No Sabio reaction should match the query, and therefore Sabio will gather data from any reaction in 
        #the same EC classification subclass (the first three digits of the four digit EC number)
        self.assertEqual(km_data.lift_info, 'Lifted From 3.1.3')
        #The query gathers data from many different reactions, and therefore many the reactions will have different EC numbers. 
        self.assertEqual(km_data.ec_numbers, ['3.1.3.1', '3.1.3.11', '3.1.3.13', '3.1.3.16', '3.1.3.17', '3.1.3.18', '3.1.3.2', '3.1.3.22', '3.1.3.23', '3.1.3.24', '3.1.3.25', '3.1.3.29', '3.1.3.3', '3.1.3.31', '3.1.3.33', '3.1.3.34', '3.1.3.36', '3.1.3.37', '3.1.3.4', '3.1.3.43', '3.1.3.46', '3.1.3.48', '3.1.3.5', '3.1.3.56', '3.1.3.57', '3.1.3.6', '3.1.3.64', '3.1.3.7', '3.1.3.74', '3.1.3.76', '3.1.3.8', '3.1.3.9', '3.1.3.91', '3.1.3.93', '3.1.3.95']) #it has two EC numbers for the same reaction. EC numbers are not a perfect system
        #Similarly, the reactions gatheres will have different Sabio reaction IDs
        self.assertEqual(km_data.sabio_reaction_ids, ['10261', '10263', '10269', '10270', '10271', '10330', '10377', '10378', '10388', '10394', '10399', '10400', '10401', '10631', '10700', '10707', '10711', '10726', '10734', '1117', '1118', '11265', '11266', '11267', '11288', '11318', '1150', '1182', '12091', '12166', '124', '12404', '12405', '129', '13054', '13055', '1325', '13296', '13453', '13635', '13655', '13970', '14184', '14198', '14219', '1422', '1423', '1424', '1683', '1686', '1851', '186', '1891', '1923', '1965', '200', '207', '2084', '209', '210', '211', '2219', '2424', '246', '2493', '2503', '27', '2717', '2730', '282', '2874', '295', '304', '305', '309', '3140', '3170', '3177', '3197', '3198', '3224', '3227', '339', '4095', '433', '444', '479', '480', '487', '501', '508', '581', '6096', '6098', '6445', '6446', '664', '6734', '6769', '697', '7185', '7186', '7187', '7188', '7189', '7190', '7191', '7192', '7193', '7194', '7195', '7196', '7278', '728', '736', '742', '75', '76', '7712', '7713', '7952', '7953', '7954', '796', '797', '7970', '7971', '7972', '7973', '7974', '7975', '7976', '7977', '8072', '8082', '816', '8203', '8204', '8205', '8206', '8749', '8902', '8953', '9216', '937', '9556', '9571', '9583', '9586', '9599', '9603', '964', '983', '9864', '9910', '9945', '9953', '9954', '9955', '9956', '9957', '9958', '9959', '9960', '9963'])
        self.assertEqual(km_data.reaction_list, ["1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3',4',5'-trisphosphate) + H2O = 1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3',4'-bisphosphate) + Phosphate", "1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3',5'-bisphosphate) + H2O = 1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3'-phosphate) + Phosphate", "1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 4',5'-bisphosphate) + H2O = 1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 4'-phosphate) + Phosphate", '1-Phosphatidyl-D-myo-inositol 4,5-bisphosphate + H2O = Phosphate + 1-Phosphatidyl-1D-myo-inositol 4-phosphate', '1D-myo-Inositol 1,4,5,6-tetrakisphosphate + H2O = Inositol 1,4,6-trisphosphate + Phosphate', '1D-myo-Inositol 2-phosphate + H2O = myo-Inositol + Phosphate', "2'-Deoxyadenosine 3'-phosphate + H2O = Deoxyadenosine + Phosphate", "2'-Deoxycytidine 3'-phosphate + H2O = 2'-Deoxycytidine + Phosphate", "2'-Deoxyguanosine 3'-phosphate + H2O = Deoxyguanosine + Phosphate", "2'-Deoxyuridine 3'-phosphate + H2O = 2'-Deoxyuridine + Phosphate", '2,3-Diphosphoglycerate + H2O = Unknown + Phosphate', '2-Deoxy-D-glucose 6-phosphate + H2O = 2-Deoxy-D-glucose + Phosphate', '2-Deoxyglucose 6-phosphate + H2O = 2-Deoxyglucose + Phosphate', "3'-Phosphoadenylyl sulfate + H2O = Phosphate + Adenylylsulfate", '4-Pyridoxic acid 5-phosphate + H2O = 4-Pyridoxate + Phosphate', "5'-Ribonucleotide + H2O = Nucleoside + Phosphate", 'ADP + H2O = Diphosphate + Adenosine', 'ADP + H2O = Phosphate + AMP', 'ATP + H2O = Adenosine + Triphosphate', 'ATP + H2O = Phosphate + ADP', "Adenosine 3',5'-bisphosphate + H2O = Phosphate + AMP", 'Arabinose 5-phosphate + H2O = Arabinose + Phosphate', 'CDP + H2O = CMP + Phosphate', "Cytidine 3'-phosphate + H2O = Cytidine + Phosphate", 'D-Fructose 1,6-bisphosphate + H2O = D-Fructose monophosphate + Phosphate', 'D-Glucose 1-phosphate + H2O = D-Glucose + Phosphate', 'D-Glucose 6-phosphate + H2O = D-Glucose + Phosphate', 'D-Mannitol 1-phosphate + H2O = Phosphate + D-Mannitol', 'D-myo-Inositol 1,2,4,5,6-pentakisphosphate + H2O = Inositol 1,2,4,6-tetrakisphosphate + Phosphate', 'Diphosphate + H2O = Phosphate', 'Fructose 1,6-bisphosphate + H2O = Unknown + Phosphate', 'Fructose 1-phosphate + H2O = Fructose + Phosphate', 'Fructose 6-phosphate + H2O = Fructose + Phosphate', 'GDP + H2O = Phosphate + GMP', 'Galactose 1-phosphate + H2O = Galactose + Phosphate', 'Glucose + 3-Phospho-D-glyceroyl phosphate = Glucose 6-phosphate + 3-Phospho-D-glycerate', 'Glucose 1-phosphate + H2O = Glucose + Phosphate', 'Glycerate 3-phosphate + H2O = Glycerate + Phosphate', 'Glycerate 3-phosphate + H2O = Phosphate + D-Glycerate', 'Glycerol 3-phosphate + H2O = Glycerol + Phosphate', 'H2O + 1-Phosphatidyl-1D-myo-inositol 3,5-bisphosphate = Phosphate + 1-Phosphatidyl-1D-myo-inositol 5-phosphate', 'H2O + 1-Phosphatidyl-1D-myo-inositol 3-phosphate = Phosphate + 1-Phosphatidyl-D-myo-inositol', 'H2O + 1-Phospho-D-glycerate = D-Glycerate + Phosphate', 'H2O + 1D-myo-Inositol 1,3,4,5-tetrakisphosphate = Phosphate + 1D-myo-Inositol 1,3,4-trisphosphate', 'H2O + 1D-myo-Inositol 1,3,4-trisphosphate = Phosphate + D-myo-Inositol 3,4-bisphosphate', 'H2O + 1D-myo-Inositol 1,4-bisphosphate = Phosphate + myo-Inositol 4-phosphate', 'H2O + 1D-myo-Inositol 1-phosphate = Phosphate + myo-Inositol', 'H2O + 1D-myo-Inositol 3-phosphate = Phosphate + myo-Inositol', 'H2O + 2-Chloro-4-nitrophenyl phosphate = Phosphate + 2-Chloro-4-nitrophenol', 'H2O + 2-Phospho-D-glycerate = Phosphate + D-Glycerate', 'H2O + 2-Phosphoglycolate = Phosphate + Glycolate', 'H2O + 3-O-Methylfluorescein phosphate = Phosphate + 3-O-Methylfluorescein', 'H2O + 4-Chlorophenyl phosphate = Phosphate + 4-Chlorophenol', 'H2O + 4-Cyanophenyl phosphate = Phosphate + 4-Cyanophenol', 'H2O + 4-Methylumbelliferyl phosphate = Phosphate + 4-Methylumbelliferone', 'H2O + 4-Nitrophenyl phenyl phosphonate = p-Nitrophenol + Phenyl phosphonate', 'H2O + 4-Nitrophenyl phosphate = Phosphate + p-Nitrophenol', 'H2O + 4-Trifluoromethylphenyl phosphate = Phosphate + 4-Trifluoromethylphenol', "H2O + 5'-Phosphopolynucleotide = Phosphate + Polynucleotide", 'H2O + 5-Fluoro-4-methylumbelliferyl phosphate = Phosphate + 5-Fluoro-4-methylumbelliferone', 'H2O + 6,8-Difluoro-4-methylumbelliferyl phosphate = Phosphate + 6,8-Difluoro-4-methylumbelliferone', 'H2O + 6-Fluoro-4-methylumbelliferyl phosphate = Phosphate + 6-Fluoro-4-methylumbelliferone', 'H2O + 6-Phosphogluconate = Phosphate + Gluconate', "H2O + 7-Methylguanosine 5'-phosphate = Phosphate + 7-Methylguanosine", 'H2O + 8-Fluoro-4-methylumbelliferyl phosphate = Phosphate + 8-Fluoro-4-methylumbelliferone', 'H2O + AMP = Phosphate + Adenosine', "H2O + Adenosine 2'-phosphate = Adenosine + Phosphate", "H2O + Adenosine 3'-phosphate = Phosphate + Adenosine", 'H2O + Bis-4-nitrophenyl phosphate = p-Nitrophenol + 4-Nitrophenyl phosphate', 'H2O + CMP = Cytidine + Phosphate', 'H2O + CTP = CDP + Phosphate', 'H2O + Casein kinase I epsilon phosphorylated = Phosphate + Casein kinase I epsilon', 'H2O + Ceramide 1-phosphate = Phosphate + N-Acylsphingosine', 'H2O + Choline phosphate = Phosphate + Choline', 'H2O + D-Fructose 1,6-bisphosphate = Phosphate + D-Fructose 6-phosphate', 'H2O + D-Fructose 1,6-bisphosphate = Phosphate + beta-D-Fructose 6-phosphate', 'H2O + D-Fructose 2,6-bisphosphate = Phosphate + D-Fructose 6-phosphate', 'H2O + D-Fructose 2,6-bisphosphate = Phosphate + beta-D-Fructose 6-phosphate', 'H2O + D-Galactose 1-phosphate = Phosphate + D-Galactose', 'H2O + D-Mannose 6-phosphate = D-Mannose + Phosphate', 'H2O + D-O-Phosphoserine = Phosphate + D-Serine', 'H2O + D-myo-Inositol 1,3-bisphosphate = Phosphate + 1D-myo-Inositol 1-phosphate', 'H2O + D-myo-Inositol 1,4,5-trisphosphate = Phosphate + 1D-myo-Inositol 1,4-bisphosphate', "H2O + Deoxythymidine 3',5'-diphosphate = Phosphate + dTMP", "H2O + Deoxythymidine 3'-phosphate = Phosphate + Thymidine", 'H2O + Ethanolamine phosphate = Phosphate + Ethanolamine', 'H2O + Fructose 1,6-bisphosphate = Phosphate + Fructose 6-phosphate', 'H2O + GMP = Phosphate + Guanosine', 'H2O + GTP = GDP + Phosphate', 'H2O + Glucosamine 6-phosphate = Phosphate + Glucosamine', 'H2O + Glucose 6-phosphate = Glucose + Phosphate', 'H2O + Glycerate 2,3-bisphosphate = Glycerate 3-phosphate + Phosphate', 'H2O + Glycerate 2,3-bisphosphate = Phosphate + 3-Phospho-D-glycerate', 'H2O + Glycerol 2-phosphate = Phosphate + Glycerol', "H2O + Guanosine 3'-phosphate = Phosphate + Guanosine", 'H2O + IMP = Phosphate + Inosine', 'H2O + Inositol 1,2,3,4,5-pentakisphosphate = Phosphate + Inositol 1,2,3,4-tetrakisphosphate', 'H2O + Inositol 1,2,4,5-tetrakisphosphate = Inositol 1,2,4-trisphosphate + Phosphate', 'H2O + Inositol 4,5-bisphosphate = Phosphate + myo-Inositol 4-phosphate', 'H2O + L-Galactose 1-phosphate = Phosphate + L-Galactose', 'H2O + L-Phosphotyrosine = Phosphate + L-Tyrosine', 'H2O + L-Threonine O-3-phosphate = L-Threonine + Phosphate', 'H2O + N-(5-Phospho-4-pyridoxyl)glycine = 4-Pyridoxylglycine + Phosphate', 'H2O + N-Acetylneuraminate 9-phosphate = Phosphate + N-Acetylneuraminate', 'H2O + N-Acylneuraminate 9-phosphate = Phosphate + N-Acylneuraminate', 'H2O + NADP+ = Phosphate + NAD+', 'H2O + NADPH = NADH + Phosphate', 'H2O + O-Phospho-L-serine = L-Serine + Phosphate', 'H2O + O-Phospho-L-serine = Serine + Phosphate', 'H2O + O-Phospho-tau-protein = Phosphate + tau-Protein', 'H2O + Phenolic phosphate = Phosphate + Phenol', 'H2O + Phenolphthalein diphosphate = Phosphate + Phenolphthalein', 'H2O + Phosphoenolpyruvate = Phosphate + Pyruvate', 'H2O + Phosphotyrosine = Phosphate + Tyrosine', 'H2O + Propan-1-ol 2-phosphate = Phosphate + 1-Propanol', 'H2O + Pyridoxal phosphate = Phosphate + Pyridoxal', 'H2O + Ribulose 5-phosphate = Ribulose + Phosphate', 'H2O + Sorbitol 6-phosphate = Phosphate + Sorbitol', 'H2O + Sphinganine 1-phosphate = Phosphate + Sphinganine', 'H2O + SpoIIAA-phosphorylated = Phosphate + SpoIIAA', 'H2O + Sucrose 6-phosphate = Sucrose + Phosphate', 'H2O + TDP = Phosphate + Thiamine monophosphate', 'H2O + Tetrapolyphosphate = Phosphate + Triphosphate', 'H2O + UDP = UMP + Phosphate', "H2O + Uridine 3'-phosphate = Uridine + Phosphate", 'H2O + XMP = Phosphate + Xanthosine', 'H2O + beta-D-Fructose 2,6-bisphosphate = Phosphate + D-Fructose 6-phosphate', 'H2O + beta-Naphthyl phosphate = Phosphate + beta-Naphthol', 'H2O + dGDP = Phosphate + dGMP', 'H2O + dIMP = Phosphate + Deoxyinosine', 'H2O + dTMP = Phosphate + Thymidine', 'H2O + dTTP = Phosphate + TDP', 'H2O + myo-Inositol 4-phosphate = Phosphate + myo-Inositol', 'H2O + myo-Inositol hexakisphosphate = Phosphate + D-myo-Inositol 1,2,4,5,6-pentakisphosphate', 'H2O + myo-Inositol phosphate = myo-Inositol + Phosphate', 'H2O + o-Carboxyphenyl phosphate = Phosphate + Salicylate', 'H2O + p-Nitrophenylthymidine phosphate = p-Nitrophenol + dTMP', 'Inositol 2,4,5,6-tetrakisphosphate + H2O = Phosphate + Inositol 2,4,6-trisphosphate', 'Inositol 2,4,5-trisphosphate + H2O = Phosphate + Inositol 2,4-bisphosphate', 'Inositol 4,5,6-trisphosphate + H2O = Inositol 4,6-bisphosphate + Phosphate', 'Mannitol 1-phosphate + H2O = Mannitol + Phosphate', 'Mannose 6-phosphate + H2O = D-Mannose + Phosphate', 'N-(5-Phospho-4-pyridoxyl)benzylamine + H2O = 4-Pyridoxylbenzylamine + Phosphate', 'N-(5-Phospho-4-pyridoxyl)ethanolamine + H2O = 4-Pyridoxylethanolamine + Phosphate', 'N-(5-Phospho-4-pyridoxyl)phenylalanine + H2O = 4-Pyridoxylphenylalanine + Phosphate', 'N-Acetyl-D-mannosamine 6-phosphate + H2O = N-Acetylmannosamine + Phosphate', 'N-Acetylglucosamine 6-phosphate + H2O = N-Acetylglucosamine + Phosphate', 'Phosphate + Pyridoxamine = H2O + Pyridoxamine phosphate', 'Phosphate + Pyridoxine = H2O + Pyridoxine phosphate', 'Phosphorylase a + H2O = Phosphorylase b + Phosphate', 'Proteine tyrosine phosphate + H2O = Phosphate + Protein tyrosine', 'Riboflavin-5-phosphate + H2O = Phosphate + Riboflavin', 'Ribose 5-phosphate + H2O = Phosphate + beta-D-Ribopyranose', 'Sedoheptulose 1,7-bisphosphate + H2O = Phosphate + Sedoheptulose 7-phosphate', 'Thymolphthalein monophosphate + H2O = Thymolphthalein + Phosphate', 'UMP + H2O = Uridine + Phosphate', 'UTP + H2O = Phosphate + UDP', 'alpha-Naphthyl phosphate + H2O = alpha-Naphthol + Phosphate', 'beta-D-Thiogalactopyranoside 6-phosphate + H2O = beta-D-Thiogalactopyranoside + Phosphate', 'dAMP + H2O = Phosphate + Deoxyadenosine', "dCMP + H2O = Phosphate + 2'-Deoxycytidine", 'dGMP + H2O = Phosphate + Deoxyguanosine', "dUMP + H2O = Phosphate + 2'-Deoxyuridine", 'sn-Glycerol 1-phosphate + H2O = Phosphate + Glycerol', 'sn-Glycerol 3-phosphate + H2O = Phosphate + Glycerol'])

        #Each kinetic info class contains Entry objects from the sabio_interface module. 
        #The following tests will tests ensure the KineticInfo fields contains the right Entry objects
        
        closest_entry_ids = km_data.closest_entry_ids
        self.assertTrue(len(closest_entry_ids)>9 and len(closest_entry_ids)<15)
        closest_entries = km_data.closest_entries #this is a list of the most closely related expiremental entries
        self.assertTrue(len(closest_entries)>9 and len(closest_entries) < 15) #when I made this test, there were 13 entries. If there are more than 25 entries, something probably went wrong 

        #from the set of closest_entries, the entry with median km value is the median_entry
        #the median_entry is an Entry object defined in sabio_interface
        median_entry = km_data.median_entry
        self.assertEqual(median_entry.ec_number, '3.1.3.5')
        self.assertEqual(median_entry.entry_id, '30219') #this is the specific Entry Number sabio assigns to each expiremental entry
        self.assertEqual(median_entry.km, '2.6E-6') #this is the km value for this expiremental entry
        self.assertEqual(median_entry.vmax, '') #this is the vmax for this expiremental entry
        self.assertEqual(median_entry.num_participants, [2,2]) #this is a list to record the number of substrates and products
        self.assertEqual(median_entry.species, 'Mycoplasma fermentans') #this is the species the expirement was done in
        self.assertEqual(median_entry.proximity, 1) #this is closely related the expiremental species is to the modeler's species
        self.assertEqual(median_entry.reaction_id, '295') #this is the Sabio assigned ID for this reaction in general
        
        #check the entry_id of the min and the max within km_data
        self.assertEqual(km_data.min_entry.entry_id, '30226')
        self.assertEqual(km_data.max_entry.entry_id, '30224')

        #now the vmax info fields will  be tested
        vmaxData = rxn.vmax_data
        self.assertEqual(vmaxData.ec_numbers, ['3.1.3.1', '3.1.3.11', '3.1.3.13', '3.1.3.16', '3.1.3.17', '3.1.3.18', '3.1.3.2', '3.1.3.22', '3.1.3.23', '3.1.3.24', '3.1.3.25', '3.1.3.29', '3.1.3.3', '3.1.3.31', '3.1.3.33', '3.1.3.34', '3.1.3.36', '3.1.3.37', '3.1.3.4', '3.1.3.43', '3.1.3.46', '3.1.3.48', '3.1.3.5', '3.1.3.56', '3.1.3.57', '3.1.3.6', '3.1.3.64', '3.1.3.7', '3.1.3.74', '3.1.3.76', '3.1.3.8', '3.1.3.9', '3.1.3.91', '3.1.3.93', '3.1.3.95']) #it has two EC numbers for the same reaction. EC numbers are not a perfect system
        self.assertEqual(vmaxData.sabio_reaction_ids, ['10261', '10263', '10269', '10270', '10271', '10330', '10377', '10378', '10388', '10394', '10399', '10400', '10401', '10631', '10700', '10707', '10711', '10726', '10734', '1117', '1118', '11265', '11266', '11267', '11288', '11318', '1150', '1182', '12091', '12166', '124', '12404', '12405', '129', '13054', '13055', '1325', '13296', '13453', '13635', '13655', '13970', '14184', '14198', '14219', '1422', '1423', '1424', '1683', '1686', '1851', '186', '1891', '1923', '1965', '200', '207', '2084', '209', '210', '211', '2219', '2424', '246', '2493', '2503', '27', '2717', '2730', '282', '2874', '295', '304', '305', '309', '3140', '3170', '3177', '3197', '3198', '3224', '3227', '339', '4095', '433', '444', '479', '480', '487', '501', '508', '581', '6096', '6098', '6445', '6446', '664', '6734', '6769', '697', '7185', '7186', '7187', '7188', '7189', '7190', '7191', '7192', '7193', '7194', '7195', '7196', '7278', '728', '736', '742', '75', '76', '7712', '7713', '7952', '7953', '7954', '796', '797', '7970', '7971', '7972', '7973', '7974', '7975', '7976', '7977', '8072', '8082', '816', '8203', '8204', '8205', '8206', '8749', '8902', '8953', '9216', '937', '9556', '9571', '9583', '9586', '9599', '9603', '964', '983', '9864', '9910', '9945', '9953', '9954', '9955', '9956', '9957', '9958', '9959', '9960', '9963'])
        self.assertEqual(vmaxData.lift_info, 'Lifted From 3.1.3')
        self.assertEqual(vmaxData.reaction_list, ["1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3',4',5'-trisphosphate) + H2O = 1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3',4'-bisphosphate) + Phosphate", "1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3',5'-bisphosphate) + H2O = 1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3'-phosphate) + Phosphate", "1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 4',5'-bisphosphate) + H2O = 1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 4'-phosphate) + Phosphate", '1-Phosphatidyl-D-myo-inositol 4,5-bisphosphate + H2O = Phosphate + 1-Phosphatidyl-1D-myo-inositol 4-phosphate', '1D-myo-Inositol 1,4,5,6-tetrakisphosphate + H2O = Inositol 1,4,6-trisphosphate + Phosphate', '1D-myo-Inositol 2-phosphate + H2O = myo-Inositol + Phosphate', "2'-Deoxyadenosine 3'-phosphate + H2O = Deoxyadenosine + Phosphate", "2'-Deoxycytidine 3'-phosphate + H2O = 2'-Deoxycytidine + Phosphate", "2'-Deoxyguanosine 3'-phosphate + H2O = Deoxyguanosine + Phosphate", "2'-Deoxyuridine 3'-phosphate + H2O = 2'-Deoxyuridine + Phosphate", '2,3-Diphosphoglycerate + H2O = Unknown + Phosphate', '2-Deoxy-D-glucose 6-phosphate + H2O = 2-Deoxy-D-glucose + Phosphate', '2-Deoxyglucose 6-phosphate + H2O = 2-Deoxyglucose + Phosphate', "3'-Phosphoadenylyl sulfate + H2O = Phosphate + Adenylylsulfate", '4-Pyridoxic acid 5-phosphate + H2O = 4-Pyridoxate + Phosphate', "5'-Ribonucleotide + H2O = Nucleoside + Phosphate", 'ADP + H2O = Diphosphate + Adenosine', 'ADP + H2O = Phosphate + AMP', 'ATP + H2O = Adenosine + Triphosphate', 'ATP + H2O = Phosphate + ADP', "Adenosine 3',5'-bisphosphate + H2O = Phosphate + AMP", 'Arabinose 5-phosphate + H2O = Arabinose + Phosphate', 'CDP + H2O = CMP + Phosphate', "Cytidine 3'-phosphate + H2O = Cytidine + Phosphate", 'D-Fructose 1,6-bisphosphate + H2O = D-Fructose monophosphate + Phosphate', 'D-Glucose 1-phosphate + H2O = D-Glucose + Phosphate', 'D-Glucose 6-phosphate + H2O = D-Glucose + Phosphate', 'D-Mannitol 1-phosphate + H2O = Phosphate + D-Mannitol', 'D-myo-Inositol 1,2,4,5,6-pentakisphosphate + H2O = Inositol 1,2,4,6-tetrakisphosphate + Phosphate', 'Diphosphate + H2O = Phosphate', 'Fructose 1,6-bisphosphate + H2O = Unknown + Phosphate', 'Fructose 1-phosphate + H2O = Fructose + Phosphate', 'Fructose 6-phosphate + H2O = Fructose + Phosphate', 'GDP + H2O = Phosphate + GMP', 'Galactose 1-phosphate + H2O = Galactose + Phosphate', 'Glucose + 3-Phospho-D-glyceroyl phosphate = Glucose 6-phosphate + 3-Phospho-D-glycerate', 'Glucose 1-phosphate + H2O = Glucose + Phosphate', 'Glycerate 3-phosphate + H2O = Glycerate + Phosphate', 'Glycerate 3-phosphate + H2O = Phosphate + D-Glycerate', 'Glycerol 3-phosphate + H2O = Glycerol + Phosphate', 'H2O + 1-Phosphatidyl-1D-myo-inositol 3,5-bisphosphate = Phosphate + 1-Phosphatidyl-1D-myo-inositol 5-phosphate', 'H2O + 1-Phosphatidyl-1D-myo-inositol 3-phosphate = Phosphate + 1-Phosphatidyl-D-myo-inositol', 'H2O + 1-Phospho-D-glycerate = D-Glycerate + Phosphate', 'H2O + 1D-myo-Inositol 1,3,4,5-tetrakisphosphate = Phosphate + 1D-myo-Inositol 1,3,4-trisphosphate', 'H2O + 1D-myo-Inositol 1,3,4-trisphosphate = Phosphate + D-myo-Inositol 3,4-bisphosphate', 'H2O + 1D-myo-Inositol 1,4-bisphosphate = Phosphate + myo-Inositol 4-phosphate', 'H2O + 1D-myo-Inositol 1-phosphate = Phosphate + myo-Inositol', 'H2O + 1D-myo-Inositol 3-phosphate = Phosphate + myo-Inositol', 'H2O + 2-Chloro-4-nitrophenyl phosphate = Phosphate + 2-Chloro-4-nitrophenol', 'H2O + 2-Phospho-D-glycerate = Phosphate + D-Glycerate', 'H2O + 2-Phosphoglycolate = Phosphate + Glycolate', 'H2O + 3-O-Methylfluorescein phosphate = Phosphate + 3-O-Methylfluorescein', 'H2O + 4-Chlorophenyl phosphate = Phosphate + 4-Chlorophenol', 'H2O + 4-Cyanophenyl phosphate = Phosphate + 4-Cyanophenol', 'H2O + 4-Methylumbelliferyl phosphate = Phosphate + 4-Methylumbelliferone', 'H2O + 4-Nitrophenyl phenyl phosphonate = p-Nitrophenol + Phenyl phosphonate', 'H2O + 4-Nitrophenyl phosphate = Phosphate + p-Nitrophenol', 'H2O + 4-Trifluoromethylphenyl phosphate = Phosphate + 4-Trifluoromethylphenol', "H2O + 5'-Phosphopolynucleotide = Phosphate + Polynucleotide", 'H2O + 5-Fluoro-4-methylumbelliferyl phosphate = Phosphate + 5-Fluoro-4-methylumbelliferone', 'H2O + 6,8-Difluoro-4-methylumbelliferyl phosphate = Phosphate + 6,8-Difluoro-4-methylumbelliferone', 'H2O + 6-Fluoro-4-methylumbelliferyl phosphate = Phosphate + 6-Fluoro-4-methylumbelliferone', 'H2O + 6-Phosphogluconate = Phosphate + Gluconate', "H2O + 7-Methylguanosine 5'-phosphate = Phosphate + 7-Methylguanosine", 'H2O + 8-Fluoro-4-methylumbelliferyl phosphate = Phosphate + 8-Fluoro-4-methylumbelliferone', 'H2O + AMP = Phosphate + Adenosine', "H2O + Adenosine 2'-phosphate = Adenosine + Phosphate", "H2O + Adenosine 3'-phosphate = Phosphate + Adenosine", 'H2O + Bis-4-nitrophenyl phosphate = p-Nitrophenol + 4-Nitrophenyl phosphate', 'H2O + CMP = Cytidine + Phosphate', 'H2O + CTP = CDP + Phosphate', 'H2O + Casein kinase I epsilon phosphorylated = Phosphate + Casein kinase I epsilon', 'H2O + Ceramide 1-phosphate = Phosphate + N-Acylsphingosine', 'H2O + Choline phosphate = Phosphate + Choline', 'H2O + D-Fructose 1,6-bisphosphate = Phosphate + D-Fructose 6-phosphate', 'H2O + D-Fructose 1,6-bisphosphate = Phosphate + beta-D-Fructose 6-phosphate', 'H2O + D-Fructose 2,6-bisphosphate = Phosphate + D-Fructose 6-phosphate', 'H2O + D-Fructose 2,6-bisphosphate = Phosphate + beta-D-Fructose 6-phosphate', 'H2O + D-Galactose 1-phosphate = Phosphate + D-Galactose', 'H2O + D-Mannose 6-phosphate = D-Mannose + Phosphate', 'H2O + D-O-Phosphoserine = Phosphate + D-Serine', 'H2O + D-myo-Inositol 1,3-bisphosphate = Phosphate + 1D-myo-Inositol 1-phosphate', 'H2O + D-myo-Inositol 1,4,5-trisphosphate = Phosphate + 1D-myo-Inositol 1,4-bisphosphate', "H2O + Deoxythymidine 3',5'-diphosphate = Phosphate + dTMP", "H2O + Deoxythymidine 3'-phosphate = Phosphate + Thymidine", 'H2O + Ethanolamine phosphate = Phosphate + Ethanolamine', 'H2O + Fructose 1,6-bisphosphate = Phosphate + Fructose 6-phosphate', 'H2O + GMP = Phosphate + Guanosine', 'H2O + GTP = GDP + Phosphate', 'H2O + Glucosamine 6-phosphate = Phosphate + Glucosamine', 'H2O + Glucose 6-phosphate = Glucose + Phosphate', 'H2O + Glycerate 2,3-bisphosphate = Glycerate 3-phosphate + Phosphate', 'H2O + Glycerate 2,3-bisphosphate = Phosphate + 3-Phospho-D-glycerate', 'H2O + Glycerol 2-phosphate = Phosphate + Glycerol', "H2O + Guanosine 3'-phosphate = Phosphate + Guanosine", 'H2O + IMP = Phosphate + Inosine', 'H2O + Inositol 1,2,3,4,5-pentakisphosphate = Phosphate + Inositol 1,2,3,4-tetrakisphosphate', 'H2O + Inositol 1,2,4,5-tetrakisphosphate = Inositol 1,2,4-trisphosphate + Phosphate', 'H2O + Inositol 4,5-bisphosphate = Phosphate + myo-Inositol 4-phosphate', 'H2O + L-Galactose 1-phosphate = Phosphate + L-Galactose', 'H2O + L-Phosphotyrosine = Phosphate + L-Tyrosine', 'H2O + L-Threonine O-3-phosphate = L-Threonine + Phosphate', 'H2O + N-(5-Phospho-4-pyridoxyl)glycine = 4-Pyridoxylglycine + Phosphate', 'H2O + N-Acetylneuraminate 9-phosphate = Phosphate + N-Acetylneuraminate', 'H2O + N-Acylneuraminate 9-phosphate = Phosphate + N-Acylneuraminate', 'H2O + NADP+ = Phosphate + NAD+', 'H2O + NADPH = NADH + Phosphate', 'H2O + O-Phospho-L-serine = L-Serine + Phosphate', 'H2O + O-Phospho-L-serine = Serine + Phosphate', 'H2O + O-Phospho-tau-protein = Phosphate + tau-Protein', 'H2O + Phenolic phosphate = Phosphate + Phenol', 'H2O + Phenolphthalein diphosphate = Phosphate + Phenolphthalein', 'H2O + Phosphoenolpyruvate = Phosphate + Pyruvate', 'H2O + Phosphotyrosine = Phosphate + Tyrosine', 'H2O + Propan-1-ol 2-phosphate = Phosphate + 1-Propanol', 'H2O + Pyridoxal phosphate = Phosphate + Pyridoxal', 'H2O + Ribulose 5-phosphate = Ribulose + Phosphate', 'H2O + Sorbitol 6-phosphate = Phosphate + Sorbitol', 'H2O + Sphinganine 1-phosphate = Phosphate + Sphinganine', 'H2O + SpoIIAA-phosphorylated = Phosphate + SpoIIAA', 'H2O + Sucrose 6-phosphate = Sucrose + Phosphate', 'H2O + TDP = Phosphate + Thiamine monophosphate', 'H2O + Tetrapolyphosphate = Phosphate + Triphosphate', 'H2O + UDP = UMP + Phosphate', "H2O + Uridine 3'-phosphate = Uridine + Phosphate", 'H2O + XMP = Phosphate + Xanthosine', 'H2O + beta-D-Fructose 2,6-bisphosphate = Phosphate + D-Fructose 6-phosphate', 'H2O + beta-Naphthyl phosphate = Phosphate + beta-Naphthol', 'H2O + dGDP = Phosphate + dGMP', 'H2O + dIMP = Phosphate + Deoxyinosine', 'H2O + dTMP = Phosphate + Thymidine', 'H2O + dTTP = Phosphate + TDP', 'H2O + myo-Inositol 4-phosphate = Phosphate + myo-Inositol', 'H2O + myo-Inositol hexakisphosphate = Phosphate + D-myo-Inositol 1,2,4,5,6-pentakisphosphate', 'H2O + myo-Inositol phosphate = myo-Inositol + Phosphate', 'H2O + o-Carboxyphenyl phosphate = Phosphate + Salicylate', 'H2O + p-Nitrophenylthymidine phosphate = p-Nitrophenol + dTMP', 'Inositol 2,4,5,6-tetrakisphosphate + H2O = Phosphate + Inositol 2,4,6-trisphosphate', 'Inositol 2,4,5-trisphosphate + H2O = Phosphate + Inositol 2,4-bisphosphate', 'Inositol 4,5,6-trisphosphate + H2O = Inositol 4,6-bisphosphate + Phosphate', 'Mannitol 1-phosphate + H2O = Mannitol + Phosphate', 'Mannose 6-phosphate + H2O = D-Mannose + Phosphate', 'N-(5-Phospho-4-pyridoxyl)benzylamine + H2O = 4-Pyridoxylbenzylamine + Phosphate', 'N-(5-Phospho-4-pyridoxyl)ethanolamine + H2O = 4-Pyridoxylethanolamine + Phosphate', 'N-(5-Phospho-4-pyridoxyl)phenylalanine + H2O = 4-Pyridoxylphenylalanine + Phosphate', 'N-Acetyl-D-mannosamine 6-phosphate + H2O = N-Acetylmannosamine + Phosphate', 'N-Acetylglucosamine 6-phosphate + H2O = N-Acetylglucosamine + Phosphate', 'Phosphate + Pyridoxamine = H2O + Pyridoxamine phosphate', 'Phosphate + Pyridoxine = H2O + Pyridoxine phosphate', 'Phosphorylase a + H2O = Phosphorylase b + Phosphate', 'Proteine tyrosine phosphate + H2O = Phosphate + Protein tyrosine', 'Riboflavin-5-phosphate + H2O = Phosphate + Riboflavin', 'Ribose 5-phosphate + H2O = Phosphate + beta-D-Ribopyranose', 'Sedoheptulose 1,7-bisphosphate + H2O = Phosphate + Sedoheptulose 7-phosphate', 'Thymolphthalein monophosphate + H2O = Thymolphthalein + Phosphate', 'UMP + H2O = Uridine + Phosphate', 'UTP + H2O = Phosphate + UDP', 'alpha-Naphthyl phosphate + H2O = alpha-Naphthol + Phosphate', 'beta-D-Thiogalactopyranoside 6-phosphate + H2O = beta-D-Thiogalactopyranoside + Phosphate', 'dAMP + H2O = Phosphate + Deoxyadenosine', "dCMP + H2O = Phosphate + 2'-Deoxycytidine", 'dGMP + H2O = Phosphate + Deoxyguanosine', "dUMP + H2O = Phosphate + 2'-Deoxyuridine", 'sn-Glycerol 1-phosphate + H2O = Phosphate + Glycerol', 'sn-Glycerol 3-phosphate + H2O = Phosphate + Glycerol'])
        #Each kinetic info class contains Entry objects from the sabio_interface module. 
        #The following tests will tests ensure the KineticInfo fields contains the right Entry objects
        
        closest_entry_ids = vmaxData.closest_entry_ids
        self.assertTrue(len(closest_entry_ids)>6 and len(closest_entry_ids)<10)
        closest_entries = vmaxData.closest_entries #this is a list of the most closely related expiremental entries
        self.assertTrue(len(closest_entries)>6 and len(closest_entries) < 10) #when I made this test, there were 45 entries. If there are more than 50 entries, something probably went wrong 

        #from the set of closest_entries, the entry with median km value is the median_entry
        #the median_entry is an Entry object defined in sabio_interface
        median_entry = vmaxData.median_entry
        self.assertEqual(median_entry.ec_number, '3.1.3.22')
        self.assertEqual(median_entry.entry_id, '22920') #this is the specific Entry Number sabio assigns to each expiremental entry
        self.assertEqual(median_entry.km, '9.0E-6') #this is the km value for this expiremental entry
        self.assertEqual(median_entry.vmax, '4.33333333E-5') #this is the vmax for this expiremental entry
        self.assertEqual(median_entry.num_participants, [2,2]) #this is a list to record the number of substrates and products
        self.assertEqual(median_entry.species, 'Streptococcus bovis') #this is the species the expirement was done in
        self.assertEqual(median_entry.proximity, 6) #this is closely related the expiremental species is to the modeler's species
        self.assertEqual(median_entry.reaction_id, '581') #this is the Sabio assigned ID for this reaction in general
        
        #check the entry_id of the min and the max within km_data
        self.assertEqual(vmaxData.min_entry.entry_id, '22921')
        self.assertEqual(vmaxData.max_entry.entry_id, '16039')

        ##################################################
        ##################################################        
        #this Entry is an example of where less than 3 relevant entries are found
        #If two relevant entries are found the min and max are filled in but the median is blank
        #if one relevant entry is found, the median is filled in, but the in and max are blank

        rxn = next((rxn for rxn in rxns if rxn.id == 'fmn_reductase'), None)
        self.assertNotEqual(rxn, None)

        self.assertEqual(rxn.reaction_ids, ['5301'])
        #both km and vmax should have KineticInfo objects in their fields
        self.assertNotEqual(rxn.km_data, None)
        self.assertNotEqual(rxn.vmax_data, None)

        #Now the KineticInfo fields will be tested
        #First is some summaries of the expiremental data. 
        km_data = rxn.km_data
        
        self.assertEqual(km_data.lift_info, 'Lift Not Used')

        #Even though all the reactions are the same, the reaction has many EC numbers (because EC numbers are an imperfect system)
        self.assertEqual(km_data.ec_numbers, ['1.14.13', '1.14.13.7', '1.5.1', '1.5.1.30', '1.5.1.39'])
        #only one reaction ID should be present
        self.assertEqual(km_data.sabio_reaction_ids, ['5301'])
        self.assertEqual(km_data.reaction_list, ['NAD+ + Reduced FMN = NADH + H+ + Riboflavin-5-phosphate'])
        #Each kinetic info class contains Entry objects from the sabio_interface module. 
        #The following tests will tests ensure the KineticInfo fields contains the right Entry objects
        
        closest_entries = km_data.closest_entries
        self.assertEqual(len(closest_entries), 2)
        #make sure that median is blank but min and max are filled in
        self.assertEqual(km_data.median_entry, None)
        self.assertFalse(km_data.min_entry==None or km_data.max_entry==None)
        
        vmaxData = rxn.vmax_data
        
        self.assertEqual(vmaxData.lift_info, 'Lift Not Used')

        #Even though all the reactions are the same, the reaction has many EC numbers (because EC numbers are an imperfect system)
        self.assertEqual(vmaxData.ec_numbers, ['1.14.13', '1.14.13.7', '1.5.1', '1.5.1.30', '1.5.1.39'])
        #only one reaction ID should be present
        self.assertEqual(vmaxData.sabio_reaction_ids, ['5301'])
        self.assertEqual(vmaxData.reaction_list, ['NAD+ + Reduced FMN = NADH + H+ + Riboflavin-5-phosphate'])
        #Each kinetic info class contains Entry objects from the sabio_interface module. 
        #The following tests will tests ensure the KineticInfo fields contains the right Entry objects
        
        closest_entries = vmaxData.closest_entries
        self.assertEqual(len(closest_entries), 1)
        #make sure that median is blank but min and max are filled in
        self.assertNotEqual(vmaxData.median_entry, None)
        self.assertTrue(vmaxData.min_entry==None or vmaxData.max_entry==None)

        
        ##################################################
        ##################################################
        #todo the next entry is a case where Sabio did find the queried reaction, however it did not find any vmax infromation
        #it only found km




class TestsCollectedFromMain(unittest.TestCase):
    @unittest.skip('Too long for typical testing')
    def test_datanator(self):
        species = 'mycoplasma pneumoniae'
        
        input_filename = path.join(path.dirname(__file__), 'fixtures', 'Mycoplasma_pneumoniae.xlsx')

        out_dir = path.join(path.dirname(__file__), 'output')
        if not path.isdir(out_dir):
            os.makedirs(out_dir)
        output_filename = path.join(out_dir, 'Mycoplasma_pneumoniae.xlsx')
        
        datanator.get_kinetic_data(input_filename, output_filename, species, proxim_limit = 8)

    def test_reaction_queries(self):
        filename = path.join(path.dirname(__file__), 'fixtures', 'twenty_reactions.xlsx')
        wb = openpyxl.load_workbook(filename=filename)
        rxns = reaction_queries.generate_reaction_queries(wb)

        for rxn in rxns:
            print(rxn.__dict__)

            for comp in rxn.substrates:
                print(comp.id)
                print(comp.inchi_smiles)
                print(comp.sabioNames)

            for comp in rxn.products:
                print(comp.id)
                print(comp.inchi_smiles)
                print(comp.sabioNames)
                
    def test_sabio_interface(self):
        query_dict = {
                    #"Organism":'"Homo sapiens"',
                    "Substrate": "AMP AND ADP", 
                    "Product": "ADP",
                    #"Enzymename":"Adk",
                    #"enzymeType":"wildtype",
                    #"Sabioreaction_id" : "5305",

                    #"Organism":'"Homo sapiens"',
                    #"Substrate": "nad",
                    #"Product": "nadh",
                    #"Enzymename":"Adk",
                    #"enzymeType":"wildtype",
                    }
        #answer =  get_sabio_data(query_dict)

        #print(answer)
        #blue = TotalResult(answer)
        
        searchString = """((Substrate:"Adenosine 3',5'-bisphosphate") AND (Substrate:"H2O" OR Substrate:"OH-")) AND ((Product:"AMP" OR Product:"Adenine-9-beta-D-arabinofuranoside 5'-monophosphate") AND (Product:"Dihydrogen phosphate" OR Product:"Phosphate"))"""
        searchString = """enzymeType:wildtype AND TemperatureRange:[30 TO 40] AND pHValueRange:[5 TO 9] AND ((Substrate:"Glyceraldehyde 3-phosphate" OR Substrate:"L-Glyceraldehyde 3-phosphate" OR Substrate:"Glycerone phosphate" OR Substrate:"D-Glyceraldehyde 3-phosphate") AND (Substrate:"D-Sedoheptulose 7-phosphate" OR Substrate:"Sedoheptulose 1-phosphate" OR Substrate:"Sedoheptulose 7-phosphate")) AND ((Product:"L-Xylulose 1-phosphate" OR Product:"D-Ribulose 5-phosphate" OR Product:"D-Xylose 5-phosphate" OR Product:"Ribose 5-phosphate" OR Product:"D-Arabinose 5-phosphate" OR Product:"D-Ribose 5-phosphate" OR Product:"D-Xylulose 1-phosphate" OR Product:"L-Xylulose 5-phosphate" OR Product:"Ribulose 5-phosphate" OR Product:"L-Ribulose 5-phosphate" OR Product:"Arabinose 5-phosphate" OR Product:"D-Xylulose 5-phosphate") AND (Product:"L-Xylulose 1-phosphate" OR Product:"D-Ribulose 5-phosphate" OR Product:"D-Xylose 5-phosphate" OR Product:"Ribose 5-phosphate" OR Product:"D-Arabinose 5-phosphate" OR Product:"D-Ribose 5-phosphate" OR Product:"D-Xylulose 1-phosphate" OR Product:"L-Xylulose 5-phosphate" OR Product:"Ribulose 5-phosphate" OR Product:"L-Ribulose 5-phosphate" OR Product:"Arabinose 5-phosphate" OR Product:"D-Xylulose 5-phosphate"))"""

        searchString = """ECNumber: ("1.2.7.0" OR "1.2.7.1" OR "1.2.7.2" OR "1.2.7.3" OR "1.2.7.4" OR "1.2.7.5" OR "1.2.7.6" OR "1.2.7.7" OR "1.2.7.8" OR "1.2.7.9" OR "1.2.7.10" OR "1.2.7.11" OR "1.2.7.12" OR "1.2.7.13" OR "1.2.7.14" OR "1.2.7.15" OR "1.2.7.16" OR "1.2.7.17" OR "1.2.7.18" OR "1.2.7.19" OR "1.2.7.20" OR "1.2.7.21" OR "1.2.7.22" OR "1.2.7.23" OR "1.2.7.24" OR "1.2.7.25" OR "1.2.7.26" OR "1.2.7.27" OR "1.2.7.28" OR "1.2.7.29" OR "1.2.7.30" OR "1.2.7.31" OR "1.2.7.32" OR "1.2.7.33" OR "1.2.7.34" OR "1.2.7.35" OR "1.2.7.36" OR "1.2.7.37" OR "1.2.7.38" OR "1.2.7.39" OR "1.2.7.40" OR "1.2.7.41" OR "1.2.7.42" OR "1.2.7.43" OR "1.2.7.44" OR "1.2.7.45" OR "1.2.7.46" OR "1.2.7.47" OR "1.2.7.48" OR "1.2.7.49" OR "1.2.7.50" OR "1.2.7.51" OR "1.2.7.52" OR "1.2.7.53" OR "1.2.7.54" OR "1.2.7.55" OR "1.2.7.56" OR "1.2.7.57" OR "1.2.7.58" OR "1.2.7.59" OR "1.2.7.60" OR "1.2.7.61" OR "1.2.7.62" OR "1.2.7.63" OR "1.2.7.64" OR "1.2.7.65" OR "1.2.7.66" OR "1.2.7.67" OR "1.2.7.68" OR "1.2.7.69" OR "1.2.7.70" OR "1.2.7.71" OR "1.2.7.72" OR "1.2.7.73" OR "1.2.7.74" OR "1.2.7.75" OR "1.2.7.76" OR "1.2.7.77" OR "1.2.7.78" OR "1.2.7.79" OR "1.2.7.80" OR "1.2.7.81" OR "1.2.7.82" OR "1.2.7.83" OR "1.2.7.84" OR "1.2.7.85" OR "1.2.7.86" OR "1.2.7.87" OR "1.2.7.88" OR "1.2.7.89" OR "1.2.7.90" OR "1.2.7.91" OR "1.2.7.92" OR "1.2.7.93" OR "1.2.7.94" OR "1.2.7.95" OR "1.2.7.96" OR "1.2.7.97" OR "1.2.7.98" OR "1.2.7.99" OR "1.2.7.100")"""
        base_species = 'mycoplasma pneumoniae'
        searchString = """((Substrate:"ATP") AND (Substrate:"UMP" OR Substrate:"Uridine 5'-phosphate")) AND ((Product:"UDP"))"""
        results =  sabio_interface.get_sabio_data(searchString, base_species)
        print(len(results.entry_list))
        
        """
        for entry in results.entry_list:
            print(entry.entry_id)
            print(entry.km)
        """
        
        """
        array = [1,2,3]
        blue = []
        blue.append(array)
        print(blue)
        """

    def test_query_string_manipulator(self):
        a = []#['dGMP', 'GMP', 'Lactose 6-phosphate', 'dGDP', "Orotidine 5'-phosphate", "Guanosine 3'-phosphate", "2',3'-Cyclic GMP", 'L-Arogenate', 'N-Acylneuraminate 9-phosphate', "Maltose 6'-phosphate", "5-Amino-6-(5'-phosphoribitylamino)uracil", '6-Phospho-beta-D-glucosyl-(1,4)-D-glucose', '2-Amino-4-hydroxy-6-(D-erythro-1,2,3-trihydroxypropyl)-7,8- dihydropteridine', '2-Amino-4-hydroxy-6-(erythro-1,2,3-trihydroxypropyl)dihydropteridine triphosphate', 'Dihydroneopterin phosphate', 'Ganciclovir', '8-Br-cGMP', "2'-Deoxyguanosine 3'-phosphate", "8-Azaguanosine-5'-monophosphate", '8-oxo-dGMP', 'Dihydroneopterin triphosphate', '8-oxo-dGTP', "2'-Deoxy-8-hydroxyguanosine"]
        b = ['dGMP', 'GMP', 'Lactose 6-phosphate', 'dGDP', "Orotidine 5'-phosphate", "Guanosine 3'-phosphate", "2',3'-Cyclic GMP", 'L-Arogenate', 'N-Acylneuraminate 9-phosphate', "Maltose 6'-phosphate", "5-Amino-6-(5'-phosphoribitylamino)uracil", '6-Phospho-beta-D-glucosyl-(1,4)-D-glucose', '2-Amino-4-hydroxy-6-(D-erythro-1,2,3-trihydroxypropyl)-7,8- dihydropteridine', '2-Amino-4-hydroxy-6-(erythro-1,2,3-trihydroxypropyl)dihydropteridine triphosphate', 'Dihydroneopterin phosphate', 'Ganciclovir', '8-Br-cGMP', "2'-Deoxyguanosine 3'-phosphate", "8-Azaguanosine-5'-monophosphate", '8-oxo-dGMP', 'Dihydroneopterin triphosphate', '8-oxo-dGTP', "2'-Deoxy-8-hydroxyguanosine"]
        c = []#["H2O"]
        d = ["phosphate"]
        stuff = []
        stuff.append([a, c])
        stuff.append([b, d])

        blank = [[], []]
        #print(stuff)


        string = 'h_m + 1a25dhvitd2_m + o2_m + nadph_m --> h2o_m + nadp_m + 1a2425thvitd2_m'
        print(query_string_manipulator.getParsedReaction(string))

    def test_inchi_generator(self):
        ourAMP = "NC1=C2N=CN(C3OC(COP([O-])([O-])=O)C(O)C3O)C2=NC=N1"
        ourATP = "NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1"
        aspartylAMP = "NC1=C2N=CN(C3OC(COP([O-])([O-])=O)C(OC(=O)C([NH3+])CC([O-])=O)C3O)C2=NC=N1"
        choline = "C[N+](C)(C)CCO"
        phosphate = "CC1=NC=C(COP([O-])([O-])=O)C(C=O)=C1O"
        water = "InChI=1S/H2O/h1H2"
        #water = "O"
        weird = "CC1NC2=C(NC(N)=NC2=O)N1C1CC(O)C(COP([O-])([O-])=O)O1"
        m8dg = "CC1NC2=C(NC(N)=NC2=O)N1C1CC(O)C(CO)O1"
        mgdGMP = "NC1=NC(=O)C2=C(N1)N(C1CC(O)C(COP([O-])([O-])=O)O1)C(=O)[N-]2"
        e3dCMP = "CCNC1=NC(=O)N(C=C1)C1CC(O)C(COP([O-])([O-])=O)O1"
        answer = inchi_generator.generateGenericInchi("InChI=1S/C10H14N5O7P/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(22-10)1-21-23(18,19)20")
        #print(getSabioNameToInchiDict())

        print(inchi_generator.getSabioNameToInchiDict())
