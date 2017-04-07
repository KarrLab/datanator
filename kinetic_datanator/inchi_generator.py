""" 
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-06
:Copyright: 2017, Karr Lab
:License: MIT
"""

from pkg_resources import resource_filename
import openbabel
import os
import requests

#first the openbabel part

#!!!!!!!!!!!!! I need to improve the inchi formatter
#!!!!!!!!!!!!!!!!
#NEEED TO RUN MORE TESTS!!!!

def generateGenericInchi(smiles_or_inchi):
	if smiles_or_inchi[:5].lower() == "inchi":
		return trimInchi(smiles_or_inchi)

	else:
		smiles = smiles_or_inchi

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
	filename = resource_filename('kinetic_datanator', 'data/InchiToCompound.txt')
	with open(filename) as f:
		for line in f:
			values = line.split(" - ")
			#print(values)
			if values[1] != "No Inchi Found":
				#inchiCompoundTranslator[values[1][9:]] = values[2]#[9:]
				inchiCompoundTranslator[values[2][:-1]] = trimInchi(values[1][9:])
				#print(values[1][9:])
			else:
				pass
	return inchiCompoundTranslator
