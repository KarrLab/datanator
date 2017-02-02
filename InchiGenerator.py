import requests
import openbabel
from jxmlease import parse

#first the openbabel part

#!!!!!!!!!!!!! I need to improve the inchi formatter
#!!!!!!!!!!!!!!!!
#NEEED TO RUN MORE TESTS!!!!

def generateGenericInchi(smilesOrInchi):
	if smilesOrInchi[:5].lower() == "inchi":
		return trimInchi(smilesOrInchi)

	else:
		smiles = smilesOrInchi

		obConversion = openbabel.OBConversion()
		obConversion.SetInAndOutFormats("smi", "inchi")
		mol = openbabel.OBMol()
		obConversion.ReadString(mol, "{}".format(smiles))
		outMDL = obConversion.WriteString(mol)
		#remove the hydrogen part of the inchi string
		outMDL = trimInchi(outMDL)
		"""
		if "/h" in outMDL:
			end = outMDL.index("/h")
			outMDL = outMDL[:end]
		"""

		return outMDL

def trimInchi(inchi):
	if "/h" in inchi:
		end = inchi.index("/h")
		inchi = inchi[:end]
	return inchi

	


def createInchiString(inchiArray):
	searchString = ""
	for inchi in inchiArray:
		subString = "("
		for entry in inchi:
			subString = subString + "InChI:" + '"' + inchi[0] + '"' + " OR "
		subString = subString[:-4] + ")"
		searchString = searchString + subString + " AND "
	searchString = searchString[:-5]
	return searchString

def getSabioNameToInchiDict():
	inchiCompoundTranslator = {}
	with open("InchiToCompound.txt") as f:
		for line in f:
			values = line.split(" - ")
			#print values
			if values[1] != "No Inchi Found":
				#inchiCompoundTranslator[values[1][9:]] = values[2]#[9:]
				inchiCompoundTranslator[values[2][:-1]] = trimInchi(values[1][9:])
				#print values[1][9:]
			else:
				print "blue"
	return inchiCompoundTranslator


if __name__ == '__main__':
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
	answer = generateGenericInchi("InChI=1S/C10H14N5O7P/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(22-10)1-21-23(18,19)20")
	print getSabioNameToInchiDict()


