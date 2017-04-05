import sys
import os
import openpyxl
sys.path.insert(0, '/home/yosef/Desktop/KineticDatanator/mypackage')

import unittest
from KineticDatanator import TaxonFinder
#import ReactionQueries
from KineticDatanator import QueryStringManipulator
from KineticDatanator import InchiGenerator
from KineticDatanator import ECNumberFinder

from KineticDatanator import ReactionQueries
from KineticDatanator import TranslatorForSabio

from KineticDatanator import SabioInterface
from KineticDatanator import Datanator


class TestProgram(unittest.TestCase):
	def test_1(self):
		self.assertEqual(1, 1)

	#TaxonFinder
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


	#QueryStringManipulator

	def test_getQuerySearchString(self):

		#test a regular case
		a = ['dGMP', 'GMP', 'Lactose 6-phosphate', 'dGDP', "Orotidine 5'-phosphate", "Guanosine 3'-phosphate", "2',3'-Cyclic GMP", 'L-Arogenate', 'N-Acylneuraminate 9-phosphate', "Maltose 6'-phosphate", "5-Amino-6-(5'-phosphoribitylamino)uracil", '6-Phospho-beta-D-glucosyl-(1,4)-D-glucose', '2-Amino-4-hydroxy-6-(D-erythro-1,2,3-trihydroxypropyl)-7,8- dihydropteridine', '2-Amino-4-hydroxy-6-(erythro-1,2,3-trihydroxypropyl)dihydropteridine triphosphate', 'Dihydroneopterin phosphate', 'Ganciclovir', '8-Br-cGMP', "2'-Deoxyguanosine 3'-phosphate", "8-Azaguanosine-5'-monophosphate", '8-oxo-dGMP', 'Dihydroneopterin triphosphate', '8-oxo-dGTP', "2'-Deoxy-8-hydroxyguanosine"]
		b = ['dGMP', 'GMP', 'Lactose 6-phosphate', 'dGDP', "Orotidine 5'-phosphate", "Guanosine 3'-phosphate", "2',3'-Cyclic GMP", 'L-Arogenate', 'N-Acylneuraminate 9-phosphate', "Maltose 6'-phosphate", "5-Amino-6-(5'-phosphoribitylamino)uracil", '6-Phospho-beta-D-glucosyl-(1,4)-D-glucose', '2-Amino-4-hydroxy-6-(D-erythro-1,2,3-trihydroxypropyl)-7,8- dihydropteridine', '2-Amino-4-hydroxy-6-(erythro-1,2,3-trihydroxypropyl)dihydropteridine triphosphate', 'Dihydroneopterin phosphate', 'Ganciclovir', '8-Br-cGMP', "2'-Deoxyguanosine 3'-phosphate", "8-Azaguanosine-5'-monophosphate", '8-oxo-dGMP', 'Dihydroneopterin triphosphate', '8-oxo-dGTP', "2'-Deoxy-8-hydroxyguanosine"]
		c = ["H2O"]
		d = ["phosphate"]
		participants = []
		participants.append([a, c])
		participants.append([b, d])
		response = QueryStringManipulator.getQuerySearchString(participants)
		expectedAnswer = """((Substrate:"dGMP" OR Substrate:"GMP" OR Substrate:"Lactose 6-phosphate" OR Substrate:"dGDP" OR Substrate:"Orotidine 5'-phosphate" OR Substrate:"Guanosine 3'-phosphate" OR Substrate:"2',3'-Cyclic GMP" OR Substrate:"L-Arogenate" OR Substrate:"N-Acylneuraminate 9-phosphate" OR Substrate:"Maltose 6'-phosphate" OR Substrate:"5-Amino-6-(5'-phosphoribitylamino)uracil" OR Substrate:"6-Phospho-beta-D-glucosyl-(1,4)-D-glucose" OR Substrate:"2-Amino-4-hydroxy-6-(D-erythro-1,2,3-trihydroxypropyl)-7,8- dihydropteridine" OR Substrate:"2-Amino-4-hydroxy-6-(erythro-1,2,3-trihydroxypropyl)dihydropteridine triphosphate" OR Substrate:"Dihydroneopterin phosphate" OR Substrate:"Ganciclovir" OR Substrate:"8-Br-cGMP" OR Substrate:"2'-Deoxyguanosine 3'-phosphate" OR Substrate:"8-Azaguanosine-5'-monophosphate" OR Substrate:"8-oxo-dGMP" OR Substrate:"Dihydroneopterin triphosphate" OR Substrate:"8-oxo-dGTP" OR Substrate:"2'-Deoxy-8-hydroxyguanosine") AND (Substrate:"H2O")) AND ((Product:"dGMP" OR Product:"GMP" OR Product:"Lactose 6-phosphate" OR Product:"dGDP" OR Product:"Orotidine 5'-phosphate" OR Product:"Guanosine 3'-phosphate" OR Product:"2',3'-Cyclic GMP" OR Product:"L-Arogenate" OR Product:"N-Acylneuraminate 9-phosphate" OR Product:"Maltose 6'-phosphate" OR Product:"5-Amino-6-(5'-phosphoribitylamino)uracil" OR Product:"6-Phospho-beta-D-glucosyl-(1,4)-D-glucose" OR Product:"2-Amino-4-hydroxy-6-(D-erythro-1,2,3-trihydroxypropyl)-7,8- dihydropteridine" OR Product:"2-Amino-4-hydroxy-6-(erythro-1,2,3-trihydroxypropyl)dihydropteridine triphosphate" OR Product:"Dihydroneopterin phosphate" OR Product:"Ganciclovir" OR Product:"8-Br-cGMP" OR Product:"2'-Deoxyguanosine 3'-phosphate" OR Product:"8-Azaguanosine-5'-monophosphate" OR Product:"8-oxo-dGMP" OR Product:"Dihydroneopterin triphosphate" OR Product:"8-oxo-dGTP" OR Product:"2'-Deoxy-8-hydroxyguanosine") AND (Product:"phosphate"))"""
		self.assertTrue(response == expectedAnswer)

		#test a blank:
		blank = [[],[]]
		response = QueryStringManipulator.getQuerySearchString(blank)
		expectedAnswer = ""
		self.assertTrue(response == "")

		a = []
		b = ['dGMP', 'GMP', 'Lactose 6-phosphate', 'dGDP', "Orotidine 5'-phosphate", "Guanosine 3'-phosphate", "2',3'-Cyclic GMP", 'L-Arogenate', 'N-Acylneuraminate 9-phosphate', "Maltose 6'-phosphate", "5-Amino-6-(5'-phosphoribitylamino)uracil", '6-Phospho-beta-D-glucosyl-(1,4)-D-glucose', '2-Amino-4-hydroxy-6-(D-erythro-1,2,3-trihydroxypropyl)-7,8- dihydropteridine', '2-Amino-4-hydroxy-6-(erythro-1,2,3-trihydroxypropyl)dihydropteridine triphosphate', 'Dihydroneopterin phosphate', 'Ganciclovir', '8-Br-cGMP', "2'-Deoxyguanosine 3'-phosphate", "8-Azaguanosine-5'-monophosphate", '8-oxo-dGMP', 'Dihydroneopterin triphosphate', '8-oxo-dGTP', "2'-Deoxy-8-hydroxyguanosine"]
		c = []
		d = ["phosphate"]
		participants = []
		participants.append([a, c])
		participants.append([b, d])
		response = QueryStringManipulator.getQuerySearchString(participants)
		expectedAnswer = """((Product:"dGMP" OR Product:"GMP" OR Product:"Lactose 6-phosphate" OR Product:"dGDP" OR Product:"Orotidine 5'-phosphate" OR Product:"Guanosine 3'-phosphate" OR Product:"2',3'-Cyclic GMP" OR Product:"L-Arogenate" OR Product:"N-Acylneuraminate 9-phosphate" OR Product:"Maltose 6'-phosphate" OR Product:"5-Amino-6-(5'-phosphoribitylamino)uracil" OR Product:"6-Phospho-beta-D-glucosyl-(1,4)-D-glucose" OR Product:"2-Amino-4-hydroxy-6-(D-erythro-1,2,3-trihydroxypropyl)-7,8- dihydropteridine" OR Product:"2-Amino-4-hydroxy-6-(erythro-1,2,3-trihydroxypropyl)dihydropteridine triphosphate" OR Product:"Dihydroneopterin phosphate" OR Product:"Ganciclovir" OR Product:"8-Br-cGMP" OR Product:"2'-Deoxyguanosine 3'-phosphate" OR Product:"8-Azaguanosine-5'-monophosphate" OR Product:"8-oxo-dGMP" OR Product:"Dihydroneopterin triphosphate" OR Product:"8-oxo-dGTP" OR Product:"2'-Deoxy-8-hydroxyguanosine") AND (Product:"phosphate"))"""
		self.assertTrue(response == expectedAnswer)

	def test_getParsedReaction(self):
		reactionString = "[c]: Complex_Acp + ATP + HDCA ==> PPI + AMP + Complex_Acp_hdc"
		response = QueryStringManipulator.getParsedReaction(reactionString)
		expectedAnswer = [['Complex_Acp', 'ATP', 'HDCA'], ['PPI', 'AMP', 'Complex_Acp_hdc']]
		self.assertTrue(response==expectedAnswer)

		reactionString = "A + B + C ==> D + E + F"
		response = QueryStringManipulator.getParsedReaction(reactionString)
		expectedAnswer = [['A', 'B', 'C'], ['D', 'E', 'F']]
		self.assertTrue(response==expectedAnswer)

		#make sure it removes hyrdrogens
		reactionString = "[c]: Complex_Acp + H2O ==> ProtMon_MPN406 + H + Pantetheine4Phosphate"
		response = QueryStringManipulator.getParsedReaction(reactionString)
		expectedAnswer = [['Complex_Acp', 'H2O'], ['ProtMon_MPN406', 'Pantetheine4Phosphate']]
		self.assertTrue(response==expectedAnswer)


		reactionString = "cyclobutane_dTdC[c] ==> cyclobutane_dTdC[e]"
		response = QueryStringManipulator.getParsedReaction(reactionString)
		expectedAnswer = [['cyclobutane_dTdC'], ['cyclobutane_dTdC']]
		self.assertTrue(response==expectedAnswer)


		reactionString = "UDPG[c] + DAG161[m] ==> UDP[c] + H[c] + GlcDAG161[m]"
		response = QueryStringManipulator.getParsedReaction(reactionString)
		expectedAnswer = [['UDPG', 'DAG161'], ['UDP', 'GlcDAG161']]
		self.assertTrue(response==expectedAnswer)

		reactionString = 'h_m + 1a25dhvitd2_m + o2_m + nadph_m --> h2o_m + nadp_m + 1a2425thvitd2_m'
		response = QueryStringManipulator.getParsedReaction(reactionString)
		expectedAnswer = [['1a25dhvitd2_m', 'o2_m', 'nadph_m'], ['h2o_m', 'nadp_m', '1a2425thvitd2_m']]
		self.assertTrue(response==expectedAnswer)

		reactionString = 'h_m + 1a25dhvitd2_m + o2_m + nadph_m <-> h2o_m + nadp_m + 1a2425thvitd2_m'
		response = QueryStringManipulator.getParsedReaction(reactionString)
		expectedAnswer = [['1a25dhvitd2_m', 'o2_m', 'nadph_m'], ['h2o_m', 'nadp_m', '1a2425thvitd2_m']]
		self.assertTrue(response==expectedAnswer)

		reactionString = 'h_m + 1a25dhvitd2_m + o2_m + nadph_m <=> h2o_m + nadp_m + 1a2425thvitd2_m'
		response = QueryStringManipulator.getParsedReaction(reactionString)
		expectedAnswer = [['1a25dhvitd2_m', 'o2_m', 'nadph_m'], ['h2o_m', 'nadp_m', '1a2425thvitd2_m']]
		self.assertTrue(response==expectedAnswer)


		#test InchiGenerator

	def test_generateGenericInchi(self):
		#this is a dictionary of smiles strings to corresponding inchi:
		smilesAndInchi = {
		"C1=CN2C=NC3=C(C2=N1)N=CN3[C@@H](C=O)O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C=O":"""InChI=1S/C12H14N5O13P3/c18-3-8(5-27-32(23,24)30-33(25,26)29-31(20,21)22)28-9(4-19)17-7-14-10-11-13-1-2-16(11)6-15-12(10)17/h1-4,6-9H,5H2,(H,23,24)(H,25,26)(H2,20,21,22)/t8-,9+/m0/s1""",
		"CC=C/C=N/NC(=O)C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F":"""InChI=1S/C11H7F13N2O/c1-2-3-4-25-26-5(27)6(12,13)7(14,15)8(16,17)9(18,19)10(20,21)11(22,23)24/h2-4H,1H3,(H,26,27)/b3-2?,25-4+""",
		"C1=CC(=C(C=C1Cl)[N+](=O)[O-])SC2=NN=C(S2)N":"""InChI=1S/C8H5ClN4O2S2/c9-4-1-2-6(5(3-4)13(14)15)16-8-12-11-7(10)17-8/h1-3H,(H2,10,11)""",
		"O.O.O.O.O.O.O.O.O.O.O.O.[Fe].[Fe].[Fe]": "InChI=1S/3Fe.12H2O/h;;;12*1H2"
		}

		testWorked = True
		for smiles in smilesAndInchi:
			if InchiGenerator.generateGenericInchi(smiles) != InchiGenerator.trimInchi(smilesAndInchi[smiles]):
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
			if InchiGenerator.generateGenericInchi(smiles) == InchiGenerator.trimInchi(smilesAndInchi[smiles]):
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
			if InchiGenerator.generateGenericInchi(inchi) != InchiGenerator.trimInchi(inchiToGenericInchi[inchi]):
				testWorked = False
		self.assertTrue(testWorked)

	#testECNumber
	def test_getECNumber(self):


		#[c]: T3P1 + S7P <==> X5P + R5P
		subArray = [u'OC(COP([O-])([O-])=O)C=O', u'OCC(=O)C(O)C(O)C(O)C(O)COP([O-])([O-])=O']
		prodArray = [u'OCC(=O)C(O)C(O)COP([O-])([O-])=O', u'OC(COP([O-])([O-])=O)C(O)C(O)C=O']
		expectedECNum = "4.1.2"
		self.assertEqual(ECNumberFinder.getECNumber(subArray, prodArray), expectedECNum)

		#becuase the substrates and proucts have to be correctly matched, the EC number search
		#should only find results for one pairing of subArray and prodArray. Therefore if the items
		#in subArray are swapped, the result should be an empty string
		subArray = [u'OCC(=O)C(O)C(O)C(O)C(O)COP([O-])([O-])=O', u'OC(COP([O-])([O-])=O)C=O']
		prodArray = [u'OCC(=O)C(O)C(O)COP([O-])([O-])=O', u'OC(COP([O-])([O-])=O)C(O)C(O)C=O']
		expectedECNum = ""
		self.assertEqual(ECNumberFinder.getECNumber(subArray, prodArray), expectedECNum)

		#[c]: PI + INS <==> R1P + HYXN
		subArray = [u'OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)N1C=NC2=C(O)N=CN=C12', u'OP([O-])([O-])=O']
		prodArray = [u'OCC1OC(OP([O-])([O-])=O)C(O)C1O', u'O=C1NC=NC2=C1NC=N2']
		expectedECNum = "2.4.2"
		self.assertEqual(ECNumberFinder.getECNumber(subArray, prodArray), expectedECNum)


		#[c]: H2O + ho5hydantoin_dRiboseMP ==> ho5hydantoin_dRibose + PI
		subArray = [u'O', u'OC1CC(OC1COP([O-])([O-])=O)N1C(O)C(=O)NC1=O']
		prodArray = [u'OCC1OC(CC1O)N1C(O)C(=O)NC1=O', u'OP([O-])([O-])=O']
		expectedECNum = "3.1.4"
		self.assertEqual(ECNumberFinder.getECNumber(subArray, prodArray), expectedECNum)


		#THE FOLLOWING IS A BIG PROBLEM. TWO REACTIONS CAN GET DIFFERENT EC NUMBERS 
		#JUST BY CHANGING THE ORDER OF THE ITEMS IN THE ARRAY
		#[c]: e1dGMP + H2O ==> e1dG + PI
		subArray = [u'CCN1C(N)=NC2=C(N=CN2C2CC(O)C(COP([O-])([O-])=O)O2)C1=O', u'O']
		prodArray = [u'CCN1C(N)=NC2=C(N=CN2C2CC(O)C(CO)O2)C1=O', u'OP([O-])([O-])=O']
		expectedECNum = "3.1.3"
		self.assertEqual(ECNumberFinder.getECNumber(subArray, prodArray), expectedECNum)

		subArray = [u'O', u'CCN1C(N)=NC2=C(N=CN2C2CC(O)C(COP([O-])([O-])=O)O2)C1=O']
		prodArray = [u'CCN1C(N)=NC2=C(N=CN2C2CC(O)C(CO)O2)C1=O', u'OP([O-])([O-])=O']
		expectedECNum = "3.1.4"
		self.assertEqual(ECNumberFinder.getECNumber(subArray, prodArray), expectedECNum)

		#[c]: e3dCMP + H2O ==> e3dC + PI
		#at the moment, kegg's algorithm doens't recognize this reaction. 
		#It is probably because it doens't actually have an official EC number
		#The official EC database itself doesn't have it
		#(it should be in 3.1.3)
		subArray = [u'CCNC1=NC(=O)N(C=C1)C1CC(O)C(COP([O-])([O-])=O)O1', u'O']
		prodArray = [u'CCNC1=NC(=O)N(C=C1)C1CC(O)C(CO)O1', u'OP([O-])([O-])=O']
		expectedECNum = "3.1.3"
		self.assertEqual(ECNumberFinder.getECNumber(subArray, prodArray), expectedECNum)

		#try it with three
		#IleIle[e] + ATP[c] + H2O[c] ==> IleIle[c] + PI[c] + H[c] + ADP[c]
		subArray = [u'CCC(C)C([NH3+])C(=O)NC(C(C)CC)C([O-])=O', u'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1', u'O']
		prodArray = [u'CCC(C)C([NH3+])C(=O)NC(C(C)CC)C([O-])=O', u'OP([O-])([O-])=O', u'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1']
		expectedECNum = "3.6.1"
		self.assertEqual(ECNumberFinder.getECNumber(subArray, prodArray), expectedECNum)

		#things to test further:
		#1) cases where the reaction is 3 and 2, or 2 and 3. 
		#2) cases where the reaction is 3 and 3, but simply looping through each three with subInchiSmiles = subInchiSmiles[-1:] + subInchiSmiles[:-1]
		#becuase the position of the 3 is wrong, and I would have to change the position of the middle, so teh side, or soemthign

		"""

		[c]: DR1P <==> DR5P
		[u'OCC1OC(CC1O)OP([O-])([O-])=O']
		[u'OC1CC(O)C(COP([O-])([O-])=O)O1']
		2.7.1
		self.assertEqual(ECNumberFinder.getECNumber(subArray, prodArray), expectedECNum)

		ATP[c] + SerSer[e] + H2O[c] ==> SerSer[c] + PI[c] + H[c] + ADP[c]
		[u'[NH3+]C(CO)C(=O)NC(CO)C([O-])=O', u'O', u'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1']
		[u'[NH3+]C(CO)C(=O)NC(CO)C([O-])=O', u'OP([O-])([O-])=O', u'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1']
		3.6.1
		self.assertEqual(ECNumberFinder.getECNumber(subArray, prodArray), expectedECNum)

		[c]: ATP + GDP ==> AMP + H + PPGPP
		[u'NC1=NC2=C(N=CN2C2OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C2O)C(=O)N1', u'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1']
		[u'NC1=C2N=CN(C3OC(COP([O-])([O-])=O)C(O)C3O)C2=NC=N1', u'NC1=NC2=C(N=CN2C2OC(COP([O-])(=O)OP([O-])([O-])=O)C(OP([O-])(=O)OP([O-])([O-])=O)C2O)C(=O)N1']
		2.7.1
		self.assertEqual(ECNumberFinder.getECNumber(subArray, prodArray), expectedECNum)

		Why is this the same as two above it?

		[c]: AHCYS + H2O ==> HCYS + ADN
		[u'O', u'NC1=C2N=CN(C3OC(CSCCC([NH3+])C([O-])=O)C(O)C3O)C2=NC=N1']
		[u'[NH3+]C(CCS)C([O-])=O', u'NC1=C2N=CN(C3OC(CO)C(O)C3O)C2=NC=N1']
		4.4.1
		self.assertEqual(ECNumberFinder.getECNumber(subArray, prodArray), expectedECNum)

		[c]: METHF + H2O <==> FTHF10 + H
		[u'O', u'NC1=NC(=O)C2=C(NCC3CN(C=[N+]23)C2=CC=C(C=C2)C(=O)NC(CCC([O-])=O)C([O-])=O)N1']
		[u'NC1=NC(=O)C2=C(NCC(CN(C=O)C3=CC=C(C=C3)C(=O)NC(CCC([O-])=O)C([O-])=O)N2)N1']
		3.5.4

		[c]: dCMP64dTMP + (2) H2O ==> dC64dT + (2) PI
		[u'O', u'CC1=CN(C2CC(O)C(COP([O-])([O-])=O)O2)C(=O)NC1C1=[N+](C2CC(O)C(COP([O-])([O-])=O)O2)C(=O)N=C(N)[C-]1O']
		[u'CC1=CN(C2CC(O)C(CO)O2)C(=O)NC1C1=[N+](C2CC(O)C(CO)O2)C(=O)N=C(N)[C-]1O', u'OP([O-])([O-])=O']
		3.1.3

		[c]: H2O + LeuLeu ==> (2) LEU
		[u'O', u'CC(C)CC([NH3+])C(=O)NC(CC(C)C)C([O-])=O']
		[u'CC(C)CC([NH3+])C([O-])=O']
		3.5.1

		UDPG[c] + GlcGlcDAG181[m] ==> UDP[c] + H[c] + GlcGlcGlcDAG181[m]
		[u'CCCCCCCC\\C=C\\CCCCCCCC(=O)OCC(COC1OC(COC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O)OC(=O)CCCCCCC\\C=C\\CCCCCCCC', u'OCC1OC(OP([O-])(=O)OP([O-])(=O)OCC2OC(C(O)C2O)N2C=CC(=O)NC2=O)C(O)C(O)C1O']
		[u'OC1C(O)C(OC1COP([O-])(=O)OP([O-])([O-])=O)N1C=CC(=O)NC1=O', u'CCCCCCCC\\C=C\\CCCCCCCC(=O)OCC(COC1OC(COC2OC(COC3OC(CO)C(O)C(O)C3O)C(O)C(O)C2O)C(O)C(O)C1O)OC(=O)CCCCCCC\\C=C\\CCCCCCCC']
		2.4.1

		[c]: ATP + URI ==> H + ADP + UMP
		[u'OCC1OC(C(O)C1O)N1C=CC(=O)NC1=O', u'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1']
		[u'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1', u'OC1C(O)C(OC1COP([O-])([O-])=O)N1C=CC(=O)NC1=O']
		2.7.1
		"""

	def test_formatECForSabio(self):
		#this simply takes an EC number as string input, and makes a generic search string
		#by returing a string that searches for N.N.N.01, N.N.N.02, etc

		expectedString = """ECNumber: ("3.1.3.0" OR "3.1.3.1" OR "3.1.3.2" OR "3.1.3.3" OR "3.1.3.4" OR "3.1.3.5" OR "3.1.3.6" OR "3.1.3.7" OR "3.1.3.8" OR "3.1.3.9" OR "3.1.3.10" OR "3.1.3.11" OR "3.1.3.12" OR "3.1.3.13" OR "3.1.3.14" OR "3.1.3.15" OR "3.1.3.16" OR "3.1.3.17" OR "3.1.3.18" OR "3.1.3.19" OR "3.1.3.20" OR "3.1.3.21" OR "3.1.3.22" OR "3.1.3.23" OR "3.1.3.24" OR "3.1.3.25" OR "3.1.3.26" OR "3.1.3.27" OR "3.1.3.28" OR "3.1.3.29" OR "3.1.3.30" OR "3.1.3.31" OR "3.1.3.32" OR "3.1.3.33" OR "3.1.3.34" OR "3.1.3.35" OR "3.1.3.36" OR "3.1.3.37" OR "3.1.3.38" OR "3.1.3.39" OR "3.1.3.40" OR "3.1.3.41" OR "3.1.3.42" OR "3.1.3.43" OR "3.1.3.44" OR "3.1.3.45" OR "3.1.3.46" OR "3.1.3.47" OR "3.1.3.48" OR "3.1.3.49" OR "3.1.3.50" OR "3.1.3.51" OR "3.1.3.52" OR "3.1.3.53" OR "3.1.3.54" OR "3.1.3.55" OR "3.1.3.56" OR "3.1.3.57" OR "3.1.3.58" OR "3.1.3.59" OR "3.1.3.60" OR "3.1.3.61" OR "3.1.3.62" OR "3.1.3.63" OR "3.1.3.64" OR "3.1.3.65" OR "3.1.3.66" OR "3.1.3.67" OR "3.1.3.68" OR "3.1.3.69" OR "3.1.3.70" OR "3.1.3.71" OR "3.1.3.72" OR "3.1.3.73" OR "3.1.3.74" OR "3.1.3.75" OR "3.1.3.76" OR "3.1.3.77" OR "3.1.3.78" OR "3.1.3.79" OR "3.1.3.80" OR "3.1.3.81" OR "3.1.3.82" OR "3.1.3.83" OR "3.1.3.84" OR "3.1.3.85" OR "3.1.3.86" OR "3.1.3.87" OR "3.1.3.88" OR "3.1.3.89" OR "3.1.3.90" OR "3.1.3.91" OR "3.1.3.92" OR "3.1.3.93" OR "3.1.3.94" OR "3.1.3.95" OR "3.1.3.96" OR "3.1.3.97" OR "3.1.3.98" OR "3.1.3.99" OR "3.1.3.100")"""
		self.assertEqual(ECNumberFinder.formatECForSabio("3.1.3"), expectedString)

		#if an empty string is sent, it should return an empty string
		expectedString = ""
		self.assertEqual(ECNumberFinder.formatECForSabio(""), expectedString)


	#testingTranslatorForSabio

	def test_getSubstrateProductQueryString(self):
		id = "Example Reaction 1"
		#[c]: ATP + Pantetheine4Phosphate ==> DPCOA + PPI
		reaction = ReactionQueries.ReactionQuery(id)
		reaction.substrates = [ReactionQueries.Compound("ATP", sabioNames = ['ATP']), ReactionQueries.Compound("Pantetheine4Phosphate", sabioNames = ["4'-Phosphopantetheine"])]
		reaction.products = [ReactionQueries.Compound("DPCOA", sabioNames = ['Dephospho-CoA', "3'-Dephospho-CoA"]), ReactionQueries.Compound("PPI", sabioNames = ['Diphosphate'])]
		searchString = TranslatorForSabio.getSubstrateProductQueryString(reaction)
		expectedString = """((Substrate:"ATP") AND (Substrate:"4'-Phosphopantetheine")) AND ((Product:"Dephospho-CoA" OR Product:"3'-Dephospho-CoA") AND (Product:"Diphosphate"))"""
		self.assertEqual(searchString, expectedString)

		#if sabio recognizes all of the participants except for one, it should still return a string
		reaction = ReactionQueries.ReactionQuery(id)
		reaction.substrates = [ReactionQueries.Compound("ATP", sabioNames = ['ATP']), ReactionQueries.Compound("Pantetheine4Phosphate", sabioNames = [])]
		reaction.products = [ReactionQueries.Compound("DPCOA", sabioNames = ['Dephospho-CoA', "3'-Dephospho-CoA"]), ReactionQueries.Compound("PPI", sabioNames = ['Diphosphate'])]
		searchString = TranslatorForSabio.getSubstrateProductQueryString(reaction)
		expectedString = """((Substrate:"ATP")) AND ((Product:"Dephospho-CoA" OR Product:"3'-Dephospho-CoA") AND (Product:"Diphosphate"))"""
		self.assertEqual(searchString, expectedString)

		#if sabio does not recognize two or more, the code should return an empty string
		reaction = ReactionQueries.ReactionQuery(id)
		reaction.substrates = [ReactionQueries.Compound("ATP", sabioNames = ['ATP']), ReactionQueries.Compound("Pantetheine4Phosphate", sabioNames = [])]
		reaction.products = [ReactionQueries.Compound("DPCOA", sabioNames = []), ReactionQueries.Compound("PPI", sabioNames = ['Diphosphate'])]
		searchString = TranslatorForSabio.getSubstrateProductQueryString(reaction)
		expectedString = ""
		self.assertEqual(searchString, expectedString)


		id = "Example Reaction 2"
		#[c]: PAP + H2O ==> AMP + PI
		reaction = ReactionQueries.ReactionQuery(id)
		reaction.substrates = [ReactionQueries.Compound("PAP", sabioNames = ["Adenosine 3',5'-bisphosphate"]), ReactionQueries.Compound("H2O", sabioNames = ['H2O', 'OH-'])]
		reaction.products = [ReactionQueries.Compound("AMP", sabioNames = ['AMP', "Adenine-9-beta-D-arabinofuranoside 5'-monophosphate"]), ReactionQueries.Compound("PI", sabioNames = ['Dihydrogen phosphate', 'Phosphate'])]
		searchString = TranslatorForSabio.getSubstrateProductQueryString(reaction)
		expectedString ="""((Substrate:"Adenosine 3',5'-bisphosphate") AND (Substrate:"H2O" OR Substrate:"OH-")) AND ((Product:"AMP" OR Product:"Adenine-9-beta-D-arabinofuranoside 5'-monophosphate") AND (Product:"Dihydrogen phosphate" OR Product:"Phosphate"))"""
		self.assertEqual(searchString, expectedString)


#SabioInterface
	def test_getSabioData(self):
		#[c]: PAP + H2O ==> AMP + PI
		baseSpecies = 'mycoplasma pneumoniae'
		searchString = """((Substrate:"Adenosine 3',5'-bisphosphate") AND (Substrate:"H2O" OR Substrate:"OH-")) AND ((Product:"AMP" OR Product:"Adenine-9-beta-D-arabinofuranoside 5'-monophosphate") AND (Product:"Dihydrogen phosphate" OR Product:"Phosphate"))"""
		results =  SabioInterface.getSabioData(searchString, baseSpecies)
		self.assertEqual(len(results.entryList), 3)

		#one of the entries should have a km of 2.0E-6 M. 
		containsKm = False
		for entry in results.entryList:
			if entry.km == "2.0E-6":
				containsKm = True
		self.assertTrue(containsKm)

		containsVmax = False
		for entry in results.entryList:
			if entry.vmax == "1.33333333E-5":
				containsVmax = True
		self.assertTrue(containsVmax)


		baseSpecies = 'mycoplasma pneumoniae'
		searchString = """Product:ADP AND Substrate:AMP AND ADP"""
		results =  SabioInterface.getSabioData(searchString, baseSpecies)
		self.assertEqual(len(results.entryList), 77)
		self.assertFalse(len(results.entryList)==75)


	
	def test_Datanator(self):
		#Find Excel sheet with reaction data
		inputFileName = os.path.join("tests/", "ExcelToTestDatanator.xlsx")
		#turn Excel sheet into openpyxl workbook
		outputFilename = "doesntmatter.xlsx"
		species = 'mycoplasma pneumoniae'

		#formatted data list 
		#The formatted data list is a list of FormattedData objects
		#FormattedData objects store and organize the final relevant expiremental data
		#FormattedData is the central Class that should contain all the data needed to respond to a query
		#Ideally, a FormattedData objective should be comprehensive enough to pass to a user, and the user should be able 
		#to derive all the information he wants from that object

		FormattedDataList = Datanator.getKineticData(inputFileName, outputFilename, species)

		found = False
		for formattedData in FormattedDataList:
			if formattedData.id == "A Reacion Name #1 (ATP + UMP <==> UDP + ADP)":
				found = True

				#test the FormattedData fields

				#reactionIDs is there to classify how many distinct sabio-RK Ids reference the desired reaction
				#there are two reasons a reaction ID list may contain more than a single entry. Either if there is an error
				#or if there are two reactions that differ only in stereochemistry (the inchi strings and smiles cannot detect
				#stereochemical distinctions).
				#the reaction IDs for this reaction should just be 201
				self.assertEqual(formattedData.reactionIDs, ['201'])
				#both km and vmax should have KineticInfo objects in their fields
				self.assertFalse(formattedData.KmData==None)
				self.assertFalse(formattedData.VmaxData==None)

				#Now the KineticInfo fields will be tested
				#First is some summaries of the expiremental data. 
				kmData = formattedData.KmData
				self.assertEqual(kmData.ECNumbers, ['2.7.4.14', '2.7.4.22']) #it has two EC numbers for the same reaction. EC numbers are not a perfect system
				self.assertEqual(kmData.SabioReactionIDs, ['201'])
				self.assertEqual(kmData.liftInfo, "Lift Not Used")
				self.assertEqual(kmData.reactionList, ['ATP + UMP = UDP + ADP'])

				#Each kinetic info class contains Entry objects from the SabioInterface module. 
				#The following tests will tests ensure the KineticInfo fields contains the right Entry objects
				
				closestEntryIDs = kmData.closestEntryIDs
				self.assertTrue(len(closestEntryIDs)>12 and len(closestEntryIDs)<25)
				closestEntries = kmData.closestEntries #this is a list of the most closely related expiremental entries
				self.assertTrue(len(closestEntries)>12 and len(closestEntries) < 25) #when I made this test, there were 13 entries. If there are more than 25 entries, something probably went wrong 

				#from the set of closestEntries, the entry with median km value is the medianEntry
				#the medianEntry is an Entry object defined in SabioInterface
				medianEntry = kmData.medianEntry
				self.assertEqual(medianEntry.ECNumber, '2.7.4.22')
				self.assertEqual(medianEntry.entryID, '17927') #this is the specific Entry Number sabio assigns to each expiremental entry
				self.assertEqual(medianEntry.km, '1.0E-4') #this is the km value for this expiremental entry
				self.assertEqual(medianEntry.vmax, '0.00665') #this is the vmax for this expiremental entry
				self.assertEqual(medianEntry.numParticipants, [2,2]) #this is a list to record the number of substrates and products
				self.assertEqual(medianEntry.species, 'Streptococcus pneumoniae') #this is the species the expirement was done in
				self.assertEqual(medianEntry.proximity, 6) #this is closely related the expiremental species is to the modeler's species
				self.assertEqual(medianEntry.reactionID, '201') #this is the Sabio assigned ID for this reaction in general
				
				#check the entryID of the min and the max within kmData
				self.assertEqual(kmData.minEntry.entryID, '51904')
				self.assertEqual(kmData.maxEntry.entryID, '51903')


				#now the vmax infor fields will  be tested
				vmaxData = formattedData.VmaxData
				self.assertEqual(vmaxData.ECNumbers, ['2.7.4.14', '2.7.4.22']) #it has two EC numbers for the same reaction. EC numbers are not a perfect system
				self.assertEqual(vmaxData.SabioReactionIDs, ['201'])
				self.assertEqual(vmaxData.liftInfo, "Lift Not Used")
				self.assertEqual(vmaxData.reactionList, ['ATP + UMP = UDP + ADP'])

				#Each kinetic info class contains Entry objects from the SabioInterface module. 
				#The following tests will tests ensure the KineticInfo fields contains the right Entry objects
				
				closestEntryIDs = vmaxData.closestEntryIDs
				self.assertTrue(len(closestEntryIDs)>44 and len(closestEntryIDs)<50)
				closestEntries = vmaxData.closestEntries #this is a list of the most closely related expiremental entries
				self.assertTrue(len(closestEntries)>44 and len(closestEntries) < 50) #when I made this test, there were 45 entries. If there are more than 50 entries, something probably went wrong 

				#from the set of closestEntries, the entry with median km value is the medianEntry
				#the medianEntry is an Entry object defined in SabioInterface
				medianEntry = vmaxData.medianEntry
				self.assertEqual(medianEntry.ECNumber, '2.7.4.22')
				self.assertEqual(medianEntry.entryID, '17907') #this is the specific Entry Number sabio assigns to each expiremental entry
				self.assertEqual(medianEntry.km, '') #this is the km value for this expiremental entry
				self.assertEqual(medianEntry.vmax, '0.00456666667') #this is the vmax for this expiremental entry
				self.assertEqual(medianEntry.numParticipants, [2,2]) #this is a list to record the number of substrates and products
				self.assertEqual(medianEntry.species, 'Streptococcus pneumoniae') #this is the species the expirement was done in
				self.assertEqual(medianEntry.proximity, 6) #this is closely related the expiremental species is to the modeler's species
				self.assertEqual(medianEntry.reactionID, '201') #this is the Sabio assigned ID for this reaction in general
				
				#check the entryID of the min and the max within kmData
				self.assertEqual(vmaxData.minEntry.entryID, '40434')
				self.assertEqual(vmaxData.maxEntry.entryID, '17934')
		self.assertTrue(found)




		#the next entry to test is one that has both the km and vmax lifted from other reactions
		found = False
		for formattedData in FormattedDataList:
			if formattedData.id == '[c]: dCMP64dTMP + (2) H2O ==> dC64dT + (2) PI':
				found = True

				#test the FormattedData fields

				#the reaction IDs for this reaction should just be empty becuase no exact mathches were found in the Sabio Database
				self.assertEqual(formattedData.reactionIDs, [])
				#both km and vmax should have KineticInfo objects in their fields
				self.assertFalse(formattedData.KmData==None)
				self.assertFalse(formattedData.VmaxData==None)

				#Now the KineticInfo fields will be tested
				#First is some summaries of the expiremental data. 
				kmData = formattedData.KmData
				#No Sabio reaction should match the query, and therefore Sabio will gather data from any reaction in 
				#the same EC classification subclass (the first three digits of the four digit EC number)
				self.assertEqual(kmData.liftInfo, 'Lifted From 3.1.3')
				#The query gathers data from many different reactions, and therefore many the reactions will have different EC numbers. 
				self.assertEqual(kmData.ECNumbers, [u'3.1.3.1', u'3.1.3.11', u'3.1.3.13', u'3.1.3.16', u'3.1.3.17', u'3.1.3.18', u'3.1.3.2', u'3.1.3.22', u'3.1.3.23', u'3.1.3.24', u'3.1.3.25', u'3.1.3.29', u'3.1.3.3', u'3.1.3.31', u'3.1.3.33', u'3.1.3.34', u'3.1.3.36', u'3.1.3.37', u'3.1.3.4', u'3.1.3.43', u'3.1.3.46', u'3.1.3.48', u'3.1.3.5', u'3.1.3.56', u'3.1.3.57', u'3.1.3.6', u'3.1.3.64', u'3.1.3.7', u'3.1.3.74', u'3.1.3.76', u'3.1.3.8', u'3.1.3.9', u'3.1.3.91', u'3.1.3.93', u'3.1.3.95']) #it has two EC numbers for the same reaction. EC numbers are not a perfect system
				#Similarly, the reactions gatheres will have different Sabio reaction IDs
				self.assertEqual(kmData.SabioReactionIDs, [u'10261', u'10263', u'10269', u'10270', u'10271', u'10330', u'10377', u'10378', u'10388', u'10394', u'10399', u'10400', u'10401', u'10631', u'10700', u'10707', u'10711', u'10726', u'10734', u'1117', u'1118', u'11265', u'11266', u'11267', u'11288', u'11318', u'1150', u'1182', u'12091', u'12166', u'124', u'12404', u'12405', u'129', u'13054', u'13055', u'1325', u'13296', u'13453', u'13635', u'13655', u'13970', u'14184', u'14198', u'14219', u'1422', u'1423', u'1424', u'1683', u'1686', u'1851', u'186', u'1891', u'1923', u'1965', u'200', u'207', u'2084', u'209', u'210', u'211', u'2219', u'2424', u'246', u'2493', u'2503', u'27', u'2717', u'2730', u'282', u'2874', u'295', u'304', u'305', u'309', u'3140', u'3170', u'3177', u'3197', u'3198', u'3224', u'3227', u'339', u'4095', u'433', u'444', u'479', u'480', u'487', u'501', u'508', u'581', u'6096', u'6098', u'6445', u'6446', u'664', u'6734', u'6769', u'697', u'7185', u'7186', u'7187', u'7188', u'7189', u'7190', u'7191', u'7192', u'7193', u'7194', u'7195', u'7196', u'7278', u'728', u'736', u'742', u'75', u'76', u'7712', u'7713', u'7952', u'7953', u'7954', u'796', u'797', u'7970', u'7971', u'7972', u'7973', u'7974', u'7975', u'7976', u'7977', u'8072', u'8082', u'816', u'8203', u'8204', u'8205', u'8206', u'8749', u'8902', u'8953', u'9216', u'937', u'9556', u'9571', u'9583', u'9586', u'9599', u'9603', u'964', u'983', u'9864', u'9910', u'9945', u'9953', u'9954', u'9955', u'9956', u'9957', u'9958', u'9959', u'9960', u'9963'])
				self.assertEqual(kmData.reactionList, [u"1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3',4',5'-trisphosphate) + H2O = 1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3',4'-bisphosphate) + Phosphate", u"1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3',5'-bisphosphate) + H2O = 1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3'-phosphate) + Phosphate", u"1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 4',5'-bisphosphate) + H2O = 1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 4'-phosphate) + Phosphate", u'1-Phosphatidyl-D-myo-inositol 4,5-bisphosphate + H2O = Phosphate + 1-Phosphatidyl-1D-myo-inositol 4-phosphate', u'1D-myo-Inositol 1,4,5,6-tetrakisphosphate + H2O = Inositol 1,4,6-trisphosphate + Phosphate', u'1D-myo-Inositol 2-phosphate + H2O = myo-Inositol + Phosphate', u"2'-Deoxyadenosine 3'-phosphate + H2O = Deoxyadenosine + Phosphate", u"2'-Deoxycytidine 3'-phosphate + H2O = 2'-Deoxycytidine + Phosphate", u"2'-Deoxyguanosine 3'-phosphate + H2O = Deoxyguanosine + Phosphate", u"2'-Deoxyuridine 3'-phosphate + H2O = 2'-Deoxyuridine + Phosphate", u'2,3-Diphosphoglycerate + H2O = Unknown + Phosphate', u'2-Deoxy-D-glucose 6-phosphate + H2O = 2-Deoxy-D-glucose + Phosphate', u'2-Deoxyglucose 6-phosphate + H2O = 2-Deoxyglucose + Phosphate', u"3'-Phosphoadenylyl sulfate + H2O = Phosphate + Adenylylsulfate", u'4-Pyridoxic acid 5-phosphate + H2O = 4-Pyridoxate + Phosphate', u"5'-Ribonucleotide + H2O = Nucleoside + Phosphate", u'ADP + H2O = Diphosphate + Adenosine', u'ADP + H2O = Phosphate + AMP', u'ATP + H2O = Adenosine + Triphosphate', u'ATP + H2O = Phosphate + ADP', u"Adenosine 3',5'-bisphosphate + H2O = Phosphate + AMP", u'Arabinose 5-phosphate + H2O = Arabinose + Phosphate', u'CDP + H2O = CMP + Phosphate', u"Cytidine 3'-phosphate + H2O = Cytidine + Phosphate", u'D-Fructose 1,6-bisphosphate + H2O = D-Fructose monophosphate + Phosphate', u'D-Glucose 1-phosphate + H2O = D-Glucose + Phosphate', u'D-Glucose 6-phosphate + H2O = D-Glucose + Phosphate', u'D-Mannitol 1-phosphate + H2O = Phosphate + D-Mannitol', u'D-myo-Inositol 1,2,4,5,6-pentakisphosphate + H2O = Inositol 1,2,4,6-tetrakisphosphate + Phosphate', u'Diphosphate + H2O = Phosphate', u'Fructose 1,6-bisphosphate + H2O = Unknown + Phosphate', u'Fructose 1-phosphate + H2O = Fructose + Phosphate', u'Fructose 6-phosphate + H2O = Fructose + Phosphate', u'GDP + H2O = Phosphate + GMP', u'Galactose 1-phosphate + H2O = Galactose + Phosphate', u'Glucose + 3-Phospho-D-glyceroyl phosphate = Glucose 6-phosphate + 3-Phospho-D-glycerate', u'Glucose 1-phosphate + H2O = Glucose + Phosphate', u'Glycerate 3-phosphate + H2O = Glycerate + Phosphate', u'Glycerate 3-phosphate + H2O = Phosphate + D-Glycerate', u'Glycerol 3-phosphate + H2O = Glycerol + Phosphate', u'H2O + 1-Phosphatidyl-1D-myo-inositol 3,5-bisphosphate = Phosphate + 1-Phosphatidyl-1D-myo-inositol 5-phosphate', u'H2O + 1-Phosphatidyl-1D-myo-inositol 3-phosphate = Phosphate + 1-Phosphatidyl-D-myo-inositol', u'H2O + 1-Phospho-D-glycerate = D-Glycerate + Phosphate', u'H2O + 1D-myo-Inositol 1,3,4,5-tetrakisphosphate = Phosphate + 1D-myo-Inositol 1,3,4-trisphosphate', u'H2O + 1D-myo-Inositol 1,3,4-trisphosphate = Phosphate + D-myo-Inositol 3,4-bisphosphate', u'H2O + 1D-myo-Inositol 1,4-bisphosphate = Phosphate + myo-Inositol 4-phosphate', u'H2O + 1D-myo-Inositol 1-phosphate = Phosphate + myo-Inositol', u'H2O + 1D-myo-Inositol 3-phosphate = Phosphate + myo-Inositol', u'H2O + 2-Chloro-4-nitrophenyl phosphate = Phosphate + 2-Chloro-4-nitrophenol', u'H2O + 2-Phospho-D-glycerate = Phosphate + D-Glycerate', u'H2O + 2-Phosphoglycolate = Phosphate + Glycolate', u'H2O + 3-O-Methylfluorescein phosphate = Phosphate + 3-O-Methylfluorescein', u'H2O + 4-Chlorophenyl phosphate = Phosphate + 4-Chlorophenol', u'H2O + 4-Cyanophenyl phosphate = Phosphate + 4-Cyanophenol', u'H2O + 4-Methylumbelliferyl phosphate = Phosphate + 4-Methylumbelliferone', u'H2O + 4-Nitrophenyl phenyl phosphonate = p-Nitrophenol + Phenyl phosphonate', u'H2O + 4-Nitrophenyl phosphate = Phosphate + p-Nitrophenol', u'H2O + 4-Trifluoromethylphenyl phosphate = Phosphate + 4-Trifluoromethylphenol', u"H2O + 5'-Phosphopolynucleotide = Phosphate + Polynucleotide", u'H2O + 5-Fluoro-4-methylumbelliferyl phosphate = Phosphate + 5-Fluoro-4-methylumbelliferone', u'H2O + 6,8-Difluoro-4-methylumbelliferyl phosphate = Phosphate + 6,8-Difluoro-4-methylumbelliferone', u'H2O + 6-Fluoro-4-methylumbelliferyl phosphate = Phosphate + 6-Fluoro-4-methylumbelliferone', u'H2O + 6-Phosphogluconate = Phosphate + Gluconate', u"H2O + 7-Methylguanosine 5'-phosphate = Phosphate + 7-Methylguanosine", u'H2O + 8-Fluoro-4-methylumbelliferyl phosphate = Phosphate + 8-Fluoro-4-methylumbelliferone', u'H2O + AMP = Phosphate + Adenosine', u"H2O + Adenosine 2'-phosphate = Adenosine + Phosphate", u"H2O + Adenosine 3'-phosphate = Phosphate + Adenosine", u'H2O + Bis-4-nitrophenyl phosphate = p-Nitrophenol + 4-Nitrophenyl phosphate', u'H2O + CMP = Cytidine + Phosphate', u'H2O + CTP = CDP + Phosphate', u'H2O + Casein kinase I epsilon phosphorylated = Phosphate + Casein kinase I epsilon', u'H2O + Ceramide 1-phosphate = Phosphate + N-Acylsphingosine', u'H2O + Choline phosphate = Phosphate + Choline', u'H2O + D-Fructose 1,6-bisphosphate = Phosphate + D-Fructose 6-phosphate', u'H2O + D-Fructose 1,6-bisphosphate = Phosphate + beta-D-Fructose 6-phosphate', u'H2O + D-Fructose 2,6-bisphosphate = Phosphate + D-Fructose 6-phosphate', u'H2O + D-Fructose 2,6-bisphosphate = Phosphate + beta-D-Fructose 6-phosphate', u'H2O + D-Galactose 1-phosphate = Phosphate + D-Galactose', u'H2O + D-Mannose 6-phosphate = D-Mannose + Phosphate', u'H2O + D-O-Phosphoserine = Phosphate + D-Serine', u'H2O + D-myo-Inositol 1,3-bisphosphate = Phosphate + 1D-myo-Inositol 1-phosphate', u'H2O + D-myo-Inositol 1,4,5-trisphosphate = Phosphate + 1D-myo-Inositol 1,4-bisphosphate', u"H2O + Deoxythymidine 3',5'-diphosphate = Phosphate + dTMP", u"H2O + Deoxythymidine 3'-phosphate = Phosphate + Thymidine", u'H2O + Ethanolamine phosphate = Phosphate + Ethanolamine', u'H2O + Fructose 1,6-bisphosphate = Phosphate + Fructose 6-phosphate', u'H2O + GMP = Phosphate + Guanosine', u'H2O + GTP = GDP + Phosphate', u'H2O + Glucosamine 6-phosphate = Phosphate + Glucosamine', u'H2O + Glucose 6-phosphate = Glucose + Phosphate', u'H2O + Glycerate 2,3-bisphosphate = Glycerate 3-phosphate + Phosphate', u'H2O + Glycerate 2,3-bisphosphate = Phosphate + 3-Phospho-D-glycerate', u'H2O + Glycerol 2-phosphate = Phosphate + Glycerol', u"H2O + Guanosine 3'-phosphate = Phosphate + Guanosine", u'H2O + IMP = Phosphate + Inosine', u'H2O + Inositol 1,2,3,4,5-pentakisphosphate = Phosphate + Inositol 1,2,3,4-tetrakisphosphate', u'H2O + Inositol 1,2,4,5-tetrakisphosphate = Inositol 1,2,4-trisphosphate + Phosphate', u'H2O + Inositol 4,5-bisphosphate = Phosphate + myo-Inositol 4-phosphate', u'H2O + L-Galactose 1-phosphate = Phosphate + L-Galactose', u'H2O + L-Phosphotyrosine = Phosphate + L-Tyrosine', u'H2O + L-Threonine O-3-phosphate = L-Threonine + Phosphate', u'H2O + N-(5-Phospho-4-pyridoxyl)glycine = 4-Pyridoxylglycine + Phosphate', u'H2O + N-Acetylneuraminate 9-phosphate = Phosphate + N-Acetylneuraminate', u'H2O + N-Acylneuraminate 9-phosphate = Phosphate + N-Acylneuraminate', u'H2O + NADP+ = Phosphate + NAD+', u'H2O + NADPH = NADH + Phosphate', u'H2O + O-Phospho-L-serine = L-Serine + Phosphate', u'H2O + O-Phospho-L-serine = Serine + Phosphate', u'H2O + O-Phospho-tau-protein = Phosphate + tau-Protein', u'H2O + Phenolic phosphate = Phosphate + Phenol', u'H2O + Phenolphthalein diphosphate = Phosphate + Phenolphthalein', u'H2O + Phosphoenolpyruvate = Phosphate + Pyruvate', u'H2O + Phosphotyrosine = Phosphate + Tyrosine', u'H2O + Propan-1-ol 2-phosphate = Phosphate + 1-Propanol', u'H2O + Pyridoxal phosphate = Phosphate + Pyridoxal', u'H2O + Ribulose 5-phosphate = Ribulose + Phosphate', u'H2O + Sorbitol 6-phosphate = Phosphate + Sorbitol', u'H2O + Sphinganine 1-phosphate = Phosphate + Sphinganine', u'H2O + SpoIIAA-phosphorylated = Phosphate + SpoIIAA', u'H2O + Sucrose 6-phosphate = Sucrose + Phosphate', u'H2O + TDP = Phosphate + Thiamine monophosphate', u'H2O + Tetrapolyphosphate = Phosphate + Triphosphate', u'H2O + UDP = UMP + Phosphate', u"H2O + Uridine 3'-phosphate = Uridine + Phosphate", u'H2O + XMP = Phosphate + Xanthosine', u'H2O + beta-D-Fructose 2,6-bisphosphate = Phosphate + D-Fructose 6-phosphate', u'H2O + beta-Naphthyl phosphate = Phosphate + beta-Naphthol', u'H2O + dGDP = Phosphate + dGMP', u'H2O + dIMP = Phosphate + Deoxyinosine', u'H2O + dTMP = Phosphate + Thymidine', u'H2O + dTTP = Phosphate + TDP', u'H2O + myo-Inositol 4-phosphate = Phosphate + myo-Inositol', u'H2O + myo-Inositol hexakisphosphate = Phosphate + D-myo-Inositol 1,2,4,5,6-pentakisphosphate', u'H2O + myo-Inositol phosphate = myo-Inositol + Phosphate', u'H2O + o-Carboxyphenyl phosphate = Phosphate + Salicylate', u'H2O + p-Nitrophenylthymidine phosphate = p-Nitrophenol + dTMP', u'Inositol 2,4,5,6-tetrakisphosphate + H2O = Phosphate + Inositol 2,4,6-trisphosphate', u'Inositol 2,4,5-trisphosphate + H2O = Phosphate + Inositol 2,4-bisphosphate', u'Inositol 4,5,6-trisphosphate + H2O = Inositol 4,6-bisphosphate + Phosphate', u'Mannitol 1-phosphate + H2O = Mannitol + Phosphate', u'Mannose 6-phosphate + H2O = D-Mannose + Phosphate', u'N-(5-Phospho-4-pyridoxyl)benzylamine + H2O = 4-Pyridoxylbenzylamine + Phosphate', u'N-(5-Phospho-4-pyridoxyl)ethanolamine + H2O = 4-Pyridoxylethanolamine + Phosphate', u'N-(5-Phospho-4-pyridoxyl)phenylalanine + H2O = 4-Pyridoxylphenylalanine + Phosphate', u'N-Acetyl-D-mannosamine 6-phosphate + H2O = N-Acetylmannosamine + Phosphate', u'N-Acetylglucosamine 6-phosphate + H2O = N-Acetylglucosamine + Phosphate', u'Phosphate + Pyridoxamine = H2O + Pyridoxamine phosphate', u'Phosphate + Pyridoxine = H2O + Pyridoxine phosphate', u'Phosphorylase a + H2O = Phosphorylase b + Phosphate', u'Proteine tyrosine phosphate + H2O = Phosphate + Protein tyrosine', u'Riboflavin-5-phosphate + H2O = Phosphate + Riboflavin', u'Ribose 5-phosphate + H2O = Phosphate + beta-D-Ribopyranose', u'Sedoheptulose 1,7-bisphosphate + H2O = Phosphate + Sedoheptulose 7-phosphate', u'Thymolphthalein monophosphate + H2O = Thymolphthalein + Phosphate', u'UMP + H2O = Uridine + Phosphate', u'UTP + H2O = Phosphate + UDP', u'alpha-Naphthyl phosphate + H2O = alpha-Naphthol + Phosphate', u'beta-D-Thiogalactopyranoside 6-phosphate + H2O = beta-D-Thiogalactopyranoside + Phosphate', u'dAMP + H2O = Phosphate + Deoxyadenosine', u"dCMP + H2O = Phosphate + 2'-Deoxycytidine", u'dGMP + H2O = Phosphate + Deoxyguanosine', u"dUMP + H2O = Phosphate + 2'-Deoxyuridine", u'sn-Glycerol 1-phosphate + H2O = Phosphate + Glycerol', u'sn-Glycerol 3-phosphate + H2O = Phosphate + Glycerol'])

				#Each kinetic info class contains Entry objects from the SabioInterface module. 
				#The following tests will tests ensure the KineticInfo fields contains the right Entry objects
				
				closestEntryIDs = kmData.closestEntryIDs
				self.assertTrue(len(closestEntryIDs)>9 and len(closestEntryIDs)<15)
				closestEntries = kmData.closestEntries #this is a list of the most closely related expiremental entries
				self.assertTrue(len(closestEntries)>9 and len(closestEntries) < 15) #when I made this test, there were 13 entries. If there are more than 25 entries, something probably went wrong 

				#from the set of closestEntries, the entry with median km value is the medianEntry
				#the medianEntry is an Entry object defined in SabioInterface
				medianEntry = kmData.medianEntry
				self.assertEqual(medianEntry.ECNumber, '3.1.3.5')
				self.assertEqual(medianEntry.entryID, '30219') #this is the specific Entry Number sabio assigns to each expiremental entry
				self.assertEqual(medianEntry.km, '2.6E-6') #this is the km value for this expiremental entry
				self.assertEqual(medianEntry.vmax, '') #this is the vmax for this expiremental entry
				self.assertEqual(medianEntry.numParticipants, [2,2]) #this is a list to record the number of substrates and products
				self.assertEqual(medianEntry.species, 'Mycoplasma fermentans') #this is the species the expirement was done in
				self.assertEqual(medianEntry.proximity, 1) #this is closely related the expiremental species is to the modeler's species
				self.assertEqual(medianEntry.reactionID, '295') #this is the Sabio assigned ID for this reaction in general
				
				#check the entryID of the min and the max within kmData
				self.assertEqual(kmData.minEntry.entryID, '30226')
				self.assertEqual(kmData.maxEntry.entryID, '30224')


				#now the vmax info fields will  be tested
				vmaxData = formattedData.VmaxData
				self.assertEqual(vmaxData.ECNumbers, [u'3.1.3.1', u'3.1.3.11', u'3.1.3.13', u'3.1.3.16', u'3.1.3.17', u'3.1.3.18', u'3.1.3.2', u'3.1.3.22', u'3.1.3.23', u'3.1.3.24', u'3.1.3.25', u'3.1.3.29', u'3.1.3.3', u'3.1.3.31', u'3.1.3.33', u'3.1.3.34', u'3.1.3.36', u'3.1.3.37', u'3.1.3.4', u'3.1.3.43', u'3.1.3.46', u'3.1.3.48', u'3.1.3.5', u'3.1.3.56', u'3.1.3.57', u'3.1.3.6', u'3.1.3.64', u'3.1.3.7', u'3.1.3.74', u'3.1.3.76', u'3.1.3.8', u'3.1.3.9', u'3.1.3.91', u'3.1.3.93', u'3.1.3.95']) #it has two EC numbers for the same reaction. EC numbers are not a perfect system
				self.assertEqual(vmaxData.SabioReactionIDs, [u'10261', u'10263', u'10269', u'10270', u'10271', u'10330', u'10377', u'10378', u'10388', u'10394', u'10399', u'10400', u'10401', u'10631', u'10700', u'10707', u'10711', u'10726', u'10734', u'1117', u'1118', u'11265', u'11266', u'11267', u'11288', u'11318', u'1150', u'1182', u'12091', u'12166', u'124', u'12404', u'12405', u'129', u'13054', u'13055', u'1325', u'13296', u'13453', u'13635', u'13655', u'13970', u'14184', u'14198', u'14219', u'1422', u'1423', u'1424', u'1683', u'1686', u'1851', u'186', u'1891', u'1923', u'1965', u'200', u'207', u'2084', u'209', u'210', u'211', u'2219', u'2424', u'246', u'2493', u'2503', u'27', u'2717', u'2730', u'282', u'2874', u'295', u'304', u'305', u'309', u'3140', u'3170', u'3177', u'3197', u'3198', u'3224', u'3227', u'339', u'4095', u'433', u'444', u'479', u'480', u'487', u'501', u'508', u'581', u'6096', u'6098', u'6445', u'6446', u'664', u'6734', u'6769', u'697', u'7185', u'7186', u'7187', u'7188', u'7189', u'7190', u'7191', u'7192', u'7193', u'7194', u'7195', u'7196', u'7278', u'728', u'736', u'742', u'75', u'76', u'7712', u'7713', u'7952', u'7953', u'7954', u'796', u'797', u'7970', u'7971', u'7972', u'7973', u'7974', u'7975', u'7976', u'7977', u'8072', u'8082', u'816', u'8203', u'8204', u'8205', u'8206', u'8749', u'8902', u'8953', u'9216', u'937', u'9556', u'9571', u'9583', u'9586', u'9599', u'9603', u'964', u'983', u'9864', u'9910', u'9945', u'9953', u'9954', u'9955', u'9956', u'9957', u'9958', u'9959', u'9960', u'9963'])
				self.assertEqual(vmaxData.liftInfo, 'Lifted From 3.1.3')
				self.assertEqual(vmaxData.reactionList, [u"1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3',4',5'-trisphosphate) + H2O = 1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3',4'-bisphosphate) + Phosphate", u"1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3',5'-bisphosphate) + H2O = 1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3'-phosphate) + Phosphate", u"1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 4',5'-bisphosphate) + H2O = 1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 4'-phosphate) + Phosphate", u'1-Phosphatidyl-D-myo-inositol 4,5-bisphosphate + H2O = Phosphate + 1-Phosphatidyl-1D-myo-inositol 4-phosphate', u'1D-myo-Inositol 1,4,5,6-tetrakisphosphate + H2O = Inositol 1,4,6-trisphosphate + Phosphate', u'1D-myo-Inositol 2-phosphate + H2O = myo-Inositol + Phosphate', u"2'-Deoxyadenosine 3'-phosphate + H2O = Deoxyadenosine + Phosphate", u"2'-Deoxycytidine 3'-phosphate + H2O = 2'-Deoxycytidine + Phosphate", u"2'-Deoxyguanosine 3'-phosphate + H2O = Deoxyguanosine + Phosphate", u"2'-Deoxyuridine 3'-phosphate + H2O = 2'-Deoxyuridine + Phosphate", u'2,3-Diphosphoglycerate + H2O = Unknown + Phosphate', u'2-Deoxy-D-glucose 6-phosphate + H2O = 2-Deoxy-D-glucose + Phosphate', u'2-Deoxyglucose 6-phosphate + H2O = 2-Deoxyglucose + Phosphate', u"3'-Phosphoadenylyl sulfate + H2O = Phosphate + Adenylylsulfate", u'4-Pyridoxic acid 5-phosphate + H2O = 4-Pyridoxate + Phosphate', u"5'-Ribonucleotide + H2O = Nucleoside + Phosphate", u'ADP + H2O = Diphosphate + Adenosine', u'ADP + H2O = Phosphate + AMP', u'ATP + H2O = Adenosine + Triphosphate', u'ATP + H2O = Phosphate + ADP', u"Adenosine 3',5'-bisphosphate + H2O = Phosphate + AMP", u'Arabinose 5-phosphate + H2O = Arabinose + Phosphate', u'CDP + H2O = CMP + Phosphate', u"Cytidine 3'-phosphate + H2O = Cytidine + Phosphate", u'D-Fructose 1,6-bisphosphate + H2O = D-Fructose monophosphate + Phosphate', u'D-Glucose 1-phosphate + H2O = D-Glucose + Phosphate', u'D-Glucose 6-phosphate + H2O = D-Glucose + Phosphate', u'D-Mannitol 1-phosphate + H2O = Phosphate + D-Mannitol', u'D-myo-Inositol 1,2,4,5,6-pentakisphosphate + H2O = Inositol 1,2,4,6-tetrakisphosphate + Phosphate', u'Diphosphate + H2O = Phosphate', u'Fructose 1,6-bisphosphate + H2O = Unknown + Phosphate', u'Fructose 1-phosphate + H2O = Fructose + Phosphate', u'Fructose 6-phosphate + H2O = Fructose + Phosphate', u'GDP + H2O = Phosphate + GMP', u'Galactose 1-phosphate + H2O = Galactose + Phosphate', u'Glucose + 3-Phospho-D-glyceroyl phosphate = Glucose 6-phosphate + 3-Phospho-D-glycerate', u'Glucose 1-phosphate + H2O = Glucose + Phosphate', u'Glycerate 3-phosphate + H2O = Glycerate + Phosphate', u'Glycerate 3-phosphate + H2O = Phosphate + D-Glycerate', u'Glycerol 3-phosphate + H2O = Glycerol + Phosphate', u'H2O + 1-Phosphatidyl-1D-myo-inositol 3,5-bisphosphate = Phosphate + 1-Phosphatidyl-1D-myo-inositol 5-phosphate', u'H2O + 1-Phosphatidyl-1D-myo-inositol 3-phosphate = Phosphate + 1-Phosphatidyl-D-myo-inositol', u'H2O + 1-Phospho-D-glycerate = D-Glycerate + Phosphate', u'H2O + 1D-myo-Inositol 1,3,4,5-tetrakisphosphate = Phosphate + 1D-myo-Inositol 1,3,4-trisphosphate', u'H2O + 1D-myo-Inositol 1,3,4-trisphosphate = Phosphate + D-myo-Inositol 3,4-bisphosphate', u'H2O + 1D-myo-Inositol 1,4-bisphosphate = Phosphate + myo-Inositol 4-phosphate', u'H2O + 1D-myo-Inositol 1-phosphate = Phosphate + myo-Inositol', u'H2O + 1D-myo-Inositol 3-phosphate = Phosphate + myo-Inositol', u'H2O + 2-Chloro-4-nitrophenyl phosphate = Phosphate + 2-Chloro-4-nitrophenol', u'H2O + 2-Phospho-D-glycerate = Phosphate + D-Glycerate', u'H2O + 2-Phosphoglycolate = Phosphate + Glycolate', u'H2O + 3-O-Methylfluorescein phosphate = Phosphate + 3-O-Methylfluorescein', u'H2O + 4-Chlorophenyl phosphate = Phosphate + 4-Chlorophenol', u'H2O + 4-Cyanophenyl phosphate = Phosphate + 4-Cyanophenol', u'H2O + 4-Methylumbelliferyl phosphate = Phosphate + 4-Methylumbelliferone', u'H2O + 4-Nitrophenyl phenyl phosphonate = p-Nitrophenol + Phenyl phosphonate', u'H2O + 4-Nitrophenyl phosphate = Phosphate + p-Nitrophenol', u'H2O + 4-Trifluoromethylphenyl phosphate = Phosphate + 4-Trifluoromethylphenol', u"H2O + 5'-Phosphopolynucleotide = Phosphate + Polynucleotide", u'H2O + 5-Fluoro-4-methylumbelliferyl phosphate = Phosphate + 5-Fluoro-4-methylumbelliferone', u'H2O + 6,8-Difluoro-4-methylumbelliferyl phosphate = Phosphate + 6,8-Difluoro-4-methylumbelliferone', u'H2O + 6-Fluoro-4-methylumbelliferyl phosphate = Phosphate + 6-Fluoro-4-methylumbelliferone', u'H2O + 6-Phosphogluconate = Phosphate + Gluconate', u"H2O + 7-Methylguanosine 5'-phosphate = Phosphate + 7-Methylguanosine", u'H2O + 8-Fluoro-4-methylumbelliferyl phosphate = Phosphate + 8-Fluoro-4-methylumbelliferone', u'H2O + AMP = Phosphate + Adenosine', u"H2O + Adenosine 2'-phosphate = Adenosine + Phosphate", u"H2O + Adenosine 3'-phosphate = Phosphate + Adenosine", u'H2O + Bis-4-nitrophenyl phosphate = p-Nitrophenol + 4-Nitrophenyl phosphate', u'H2O + CMP = Cytidine + Phosphate', u'H2O + CTP = CDP + Phosphate', u'H2O + Casein kinase I epsilon phosphorylated = Phosphate + Casein kinase I epsilon', u'H2O + Ceramide 1-phosphate = Phosphate + N-Acylsphingosine', u'H2O + Choline phosphate = Phosphate + Choline', u'H2O + D-Fructose 1,6-bisphosphate = Phosphate + D-Fructose 6-phosphate', u'H2O + D-Fructose 1,6-bisphosphate = Phosphate + beta-D-Fructose 6-phosphate', u'H2O + D-Fructose 2,6-bisphosphate = Phosphate + D-Fructose 6-phosphate', u'H2O + D-Fructose 2,6-bisphosphate = Phosphate + beta-D-Fructose 6-phosphate', u'H2O + D-Galactose 1-phosphate = Phosphate + D-Galactose', u'H2O + D-Mannose 6-phosphate = D-Mannose + Phosphate', u'H2O + D-O-Phosphoserine = Phosphate + D-Serine', u'H2O + D-myo-Inositol 1,3-bisphosphate = Phosphate + 1D-myo-Inositol 1-phosphate', u'H2O + D-myo-Inositol 1,4,5-trisphosphate = Phosphate + 1D-myo-Inositol 1,4-bisphosphate', u"H2O + Deoxythymidine 3',5'-diphosphate = Phosphate + dTMP", u"H2O + Deoxythymidine 3'-phosphate = Phosphate + Thymidine", u'H2O + Ethanolamine phosphate = Phosphate + Ethanolamine', u'H2O + Fructose 1,6-bisphosphate = Phosphate + Fructose 6-phosphate', u'H2O + GMP = Phosphate + Guanosine', u'H2O + GTP = GDP + Phosphate', u'H2O + Glucosamine 6-phosphate = Phosphate + Glucosamine', u'H2O + Glucose 6-phosphate = Glucose + Phosphate', u'H2O + Glycerate 2,3-bisphosphate = Glycerate 3-phosphate + Phosphate', u'H2O + Glycerate 2,3-bisphosphate = Phosphate + 3-Phospho-D-glycerate', u'H2O + Glycerol 2-phosphate = Phosphate + Glycerol', u"H2O + Guanosine 3'-phosphate = Phosphate + Guanosine", u'H2O + IMP = Phosphate + Inosine', u'H2O + Inositol 1,2,3,4,5-pentakisphosphate = Phosphate + Inositol 1,2,3,4-tetrakisphosphate', u'H2O + Inositol 1,2,4,5-tetrakisphosphate = Inositol 1,2,4-trisphosphate + Phosphate', u'H2O + Inositol 4,5-bisphosphate = Phosphate + myo-Inositol 4-phosphate', u'H2O + L-Galactose 1-phosphate = Phosphate + L-Galactose', u'H2O + L-Phosphotyrosine = Phosphate + L-Tyrosine', u'H2O + L-Threonine O-3-phosphate = L-Threonine + Phosphate', u'H2O + N-(5-Phospho-4-pyridoxyl)glycine = 4-Pyridoxylglycine + Phosphate', u'H2O + N-Acetylneuraminate 9-phosphate = Phosphate + N-Acetylneuraminate', u'H2O + N-Acylneuraminate 9-phosphate = Phosphate + N-Acylneuraminate', u'H2O + NADP+ = Phosphate + NAD+', u'H2O + NADPH = NADH + Phosphate', u'H2O + O-Phospho-L-serine = L-Serine + Phosphate', u'H2O + O-Phospho-L-serine = Serine + Phosphate', u'H2O + O-Phospho-tau-protein = Phosphate + tau-Protein', u'H2O + Phenolic phosphate = Phosphate + Phenol', u'H2O + Phenolphthalein diphosphate = Phosphate + Phenolphthalein', u'H2O + Phosphoenolpyruvate = Phosphate + Pyruvate', u'H2O + Phosphotyrosine = Phosphate + Tyrosine', u'H2O + Propan-1-ol 2-phosphate = Phosphate + 1-Propanol', u'H2O + Pyridoxal phosphate = Phosphate + Pyridoxal', u'H2O + Ribulose 5-phosphate = Ribulose + Phosphate', u'H2O + Sorbitol 6-phosphate = Phosphate + Sorbitol', u'H2O + Sphinganine 1-phosphate = Phosphate + Sphinganine', u'H2O + SpoIIAA-phosphorylated = Phosphate + SpoIIAA', u'H2O + Sucrose 6-phosphate = Sucrose + Phosphate', u'H2O + TDP = Phosphate + Thiamine monophosphate', u'H2O + Tetrapolyphosphate = Phosphate + Triphosphate', u'H2O + UDP = UMP + Phosphate', u"H2O + Uridine 3'-phosphate = Uridine + Phosphate", u'H2O + XMP = Phosphate + Xanthosine', u'H2O + beta-D-Fructose 2,6-bisphosphate = Phosphate + D-Fructose 6-phosphate', u'H2O + beta-Naphthyl phosphate = Phosphate + beta-Naphthol', u'H2O + dGDP = Phosphate + dGMP', u'H2O + dIMP = Phosphate + Deoxyinosine', u'H2O + dTMP = Phosphate + Thymidine', u'H2O + dTTP = Phosphate + TDP', u'H2O + myo-Inositol 4-phosphate = Phosphate + myo-Inositol', u'H2O + myo-Inositol hexakisphosphate = Phosphate + D-myo-Inositol 1,2,4,5,6-pentakisphosphate', u'H2O + myo-Inositol phosphate = myo-Inositol + Phosphate', u'H2O + o-Carboxyphenyl phosphate = Phosphate + Salicylate', u'H2O + p-Nitrophenylthymidine phosphate = p-Nitrophenol + dTMP', u'Inositol 2,4,5,6-tetrakisphosphate + H2O = Phosphate + Inositol 2,4,6-trisphosphate', u'Inositol 2,4,5-trisphosphate + H2O = Phosphate + Inositol 2,4-bisphosphate', u'Inositol 4,5,6-trisphosphate + H2O = Inositol 4,6-bisphosphate + Phosphate', u'Mannitol 1-phosphate + H2O = Mannitol + Phosphate', u'Mannose 6-phosphate + H2O = D-Mannose + Phosphate', u'N-(5-Phospho-4-pyridoxyl)benzylamine + H2O = 4-Pyridoxylbenzylamine + Phosphate', u'N-(5-Phospho-4-pyridoxyl)ethanolamine + H2O = 4-Pyridoxylethanolamine + Phosphate', u'N-(5-Phospho-4-pyridoxyl)phenylalanine + H2O = 4-Pyridoxylphenylalanine + Phosphate', u'N-Acetyl-D-mannosamine 6-phosphate + H2O = N-Acetylmannosamine + Phosphate', u'N-Acetylglucosamine 6-phosphate + H2O = N-Acetylglucosamine + Phosphate', u'Phosphate + Pyridoxamine = H2O + Pyridoxamine phosphate', u'Phosphate + Pyridoxine = H2O + Pyridoxine phosphate', u'Phosphorylase a + H2O = Phosphorylase b + Phosphate', u'Proteine tyrosine phosphate + H2O = Phosphate + Protein tyrosine', u'Riboflavin-5-phosphate + H2O = Phosphate + Riboflavin', u'Ribose 5-phosphate + H2O = Phosphate + beta-D-Ribopyranose', u'Sedoheptulose 1,7-bisphosphate + H2O = Phosphate + Sedoheptulose 7-phosphate', u'Thymolphthalein monophosphate + H2O = Thymolphthalein + Phosphate', u'UMP + H2O = Uridine + Phosphate', u'UTP + H2O = Phosphate + UDP', u'alpha-Naphthyl phosphate + H2O = alpha-Naphthol + Phosphate', u'beta-D-Thiogalactopyranoside 6-phosphate + H2O = beta-D-Thiogalactopyranoside + Phosphate', u'dAMP + H2O = Phosphate + Deoxyadenosine', u"dCMP + H2O = Phosphate + 2'-Deoxycytidine", u'dGMP + H2O = Phosphate + Deoxyguanosine', u"dUMP + H2O = Phosphate + 2'-Deoxyuridine", u'sn-Glycerol 1-phosphate + H2O = Phosphate + Glycerol', u'sn-Glycerol 3-phosphate + H2O = Phosphate + Glycerol'])
				#Each kinetic info class contains Entry objects from the SabioInterface module. 
				#The following tests will tests ensure the KineticInfo fields contains the right Entry objects
				
				closestEntryIDs = vmaxData.closestEntryIDs
				self.assertTrue(len(closestEntryIDs)>6 and len(closestEntryIDs)<10)
				closestEntries = vmaxData.closestEntries #this is a list of the most closely related expiremental entries
				self.assertTrue(len(closestEntries)>6 and len(closestEntries) < 10) #when I made this test, there were 45 entries. If there are more than 50 entries, something probably went wrong 

				#from the set of closestEntries, the entry with median km value is the medianEntry
				#the medianEntry is an Entry object defined in SabioInterface
				medianEntry = vmaxData.medianEntry
				self.assertEqual(medianEntry.ECNumber, '3.1.3.22')
				self.assertEqual(medianEntry.entryID, '22920') #this is the specific Entry Number sabio assigns to each expiremental entry
				self.assertEqual(medianEntry.km, '9.0E-6') #this is the km value for this expiremental entry
				self.assertEqual(medianEntry.vmax, '4.33333333E-5') #this is the vmax for this expiremental entry
				self.assertEqual(medianEntry.numParticipants, [2,2]) #this is a list to record the number of substrates and products
				self.assertEqual(medianEntry.species, 'Streptococcus bovis') #this is the species the expirement was done in
				self.assertEqual(medianEntry.proximity, 6) #this is closely related the expiremental species is to the modeler's species
				self.assertEqual(medianEntry.reactionID, '581') #this is the Sabio assigned ID for this reaction in general
				
				#check the entryID of the min and the max within kmData
				self.assertEqual(vmaxData.minEntry.entryID, '22921')
				self.assertEqual(vmaxData.maxEntry.entryID, '16039')
		self.assertTrue(found)


		#this Entry is an example of where less than 3 relevant entries are found
		#If two relevant entries are found the min and max are filled in but the median is blank
		#if one relevant entry is found, the median is filled in, but the in and max are blank
		found = False
		for formattedData in FormattedDataList:
			if formattedData.id == '[c]: FMN + H + NADH ==> FMNH2 + NAD':
				found = True

				#test the FormattedData fields

				self.assertEqual(formattedData.reactionIDs, ['5301'])
				#both km and vmax should have KineticInfo objects in their fields
				self.assertFalse(formattedData.KmData==None)
				self.assertFalse(formattedData.VmaxData==None)

				#Now the KineticInfo fields will be tested
				#First is some summaries of the expiremental data. 
				kmData = formattedData.KmData
				
				self.assertEqual(kmData.liftInfo, 'Lift Not Used')

				#Even though all the reactions are the same, the reaction has many EC numbers (because EC numbers are an imperfect system)
				self.assertEqual(kmData.ECNumbers, [u'1.14.13', u'1.14.13.7', u'1.5.1', u'1.5.1.30', u'1.5.1.39'])
				#only one reaction ID should be present
				self.assertEqual(kmData.SabioReactionIDs, ['5301'])
				self.assertEqual(kmData.reactionList, [u'NAD+ + Reduced FMN = NADH + H+ + Riboflavin-5-phosphate'])
				#Each kinetic info class contains Entry objects from the SabioInterface module. 
				#The following tests will tests ensure the KineticInfo fields contains the right Entry objects
				
				closestEntries = kmData.closestEntries
				self.assertEqual(len(closestEntries), 2)
				#make sure that median is blank but min and max are filled in
				self.assertTrue(kmData.medianEntry==None)
				self.assertFalse(kmData.minEntry==None or kmData.maxEntry==None)
				
				vmaxData = formattedData.VmaxData
				
				self.assertEqual(vmaxData.liftInfo, 'Lift Not Used')

				#Even though all the reactions are the same, the reaction has many EC numbers (because EC numbers are an imperfect system)
				self.assertEqual(vmaxData.ECNumbers, [u'1.14.13', u'1.14.13.7', u'1.5.1', u'1.5.1.30', u'1.5.1.39'])
				#only one reaction ID should be present
				self.assertEqual(vmaxData.SabioReactionIDs, ['5301'])
				self.assertEqual(vmaxData.reactionList, [u'NAD+ + Reduced FMN = NADH + H+ + Riboflavin-5-phosphate'])
				#Each kinetic info class contains Entry objects from the SabioInterface module. 
				#The following tests will tests ensure the KineticInfo fields contains the right Entry objects
				
				closestEntries = vmaxData.closestEntries
				self.assertEqual(len(closestEntries), 1)
				#make sure that median is blank but min and max are filled in
				self.assertFalse(vmaxData.medianEntry==None)
				self.assertTrue(vmaxData.minEntry==None or vmaxData.maxEntry==None)
				
		self.assertTrue(found)




		#the next entry is a case where Sabio did find the queried reaction, however it did not find any vmax infromation
		#it only found km





		#medianKmEntry is an object of the Entry class found in the SabioInterface module
		#the following are a series of tests that make sure that the entry fields are
		#filled in properly
		"""
		medianKmEnry = formattedData.KmData.medianKmEnry
		self.assertFalse(medianKmEnry==None)
		medianKm = medianKmEnry.km
		self.assertEqual(medianKmEnry.entryID, 42062)
		self.assertEqual(medianKmEnry.vmax, "")
		self.assertEqual(medianKmEnry.proximity, 6)
		"""



			#print formattedData.__dict__
				
		print FormattedDataList
		self.assertEqual(1,1)



