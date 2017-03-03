import sys
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



	#there is a bug here, I need to fix this

if __name__ == '__main__':
	#main()
	#unittest.main()
	#test_getSabioData()
	TestProgram().test_getSubstrateProductQueryString()








