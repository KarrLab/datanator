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
	
	def test_Datanator(self):
		#Find Excel sheet with reaction data
		inputFileName = os.path.join("/home/yosef/Desktop/KineticDatanator/tests/", "ExcelToTestDatanator.xlsx")
		#turn Excel sheet into openpyxl workbook
		outputFilename = "doesntmatter.xlsx"
		species = 'mycoplasma pneumoniae'

		#formatted data list 
		FormattedDataList = Datanator.getKineticData(inputFileName, outputFilename, species)

		for formattedData in FormattedDataList:
			if formattedData.id == "A Reacion Name #1 (ATP + UMP <==> UDP + ADP)":

				#test the formatted data fields
				self.assertEqual(formattedData.reacionIDs, [77])


				#medianKmEntry is an object of the Entry class found in the SabioInterface module
				#the following are a series of tests that make sure that the entry fields are
				#filled in properly
				medianKmEnry = formattedData.KmData.medianKmEnry
				self.assertFalse(medianKmEnry==None)
				medianKm = medianKmEnry.km
				self.assertEqual(medianKmEnry.entryID, 42062)
				self.assertEqual(medianKmEnry.vmax, "")
				self.assertEqual(medianKmEnry.proximity, 6)



			print formattedData.__dict__
				
		print FormattedDataList
		self.assertEqual(1,1)




if __name__ == '__main__':
    unittest.main()
