from datanator.data_source import uniprot_ncbi
import datetime
import dateutil
import os
import shutil
import tempfile
import unittest
import pytest
import sqlalchemy
from sqlalchemy import create_engine
from pathlib import Path
from datanator.util import warning_util

warning_util.disable_warnings()


class TestDownloader(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		src = uniprot_ncbi
		cls.root_folder = tempfile.mkdtemp()
		cls.file_name = ["oma-uniprot", "oma-ncbi"]
		cls.file_type = ".txt.gz"
		cls.db_type = '.sqlite'

	@classmethod
	def tearDownClass(cls):
		shutil.rmtree(cls.root_folder)

	def setUp(self):
		self.src = uniprot_ncbi

	# def tearDown(self):
	# 	src = self.src
	# 	src.close()
	# 	src.dispose()

	def test_download_data(self):
		src = self.src
		src.download_data(self.file_name[0], self.file_type, self.root_folder)
		src.download_data(self.file_name[1], self.file_type, self.root_folder)
		
		storage_location1 = self.root_folder + self.file_name[0] + self.file_type
		storage_location2 = self.root_folder + self.file_name[1] + self.file_type
		p = Path(storage_location1)
		q = Path(storage_location2)
		self.assertTrue(p.exists())
		self.assertTrue(q.exists())

	def test_load_data(self):
		src = self.src
		print('file_name[0]: ' + self.file_name[0] + '\n')
		print('file_type: ' + self.file_type + '\n')
		df = src.load_data(self.file_name[0],self.file_type,self.root_folder)
		temp = df[0:6]
		col0 = self.file_name[0].split('-')[0].upper()
		col1 = self.file_name[0].split('-')[1].upper()
		self.assertEqual(temp[col0].tolist(), ['HEIAB00001', 'HEIAB00002', 'HEIAB00003','HEIAB00004','HEIAB00005','HEIAB00006'])
		self.assertEqual(temp[col1].tolist(), ['A0A2U3CKJ0', 'A0A2U3CKK2', 'A0A2U3CKK9','A0A2U3CKI6','A0A2U3CKJ5','A0A2U3CKL8'])
	
	def test_sql_data(self):
		src = self.src
		df1 = src.load_data(self.file_name[0],self.file_type,self.root_folder)
		df2 = src.load_data(self.file_name[1],self.file_type,self.root_folder)
		temp1 = df1[1:7]
		temp2 = df2[1:7]
		src.sql_data(temp1, temp2, self.file_name,self.root_folder,self.db_type)
		name1 = self.file_name[0].split('-')[1]
		name2 = self.file_name[1].split('-')[1]
		db_name = name1 + '_' + name2
		engine = sqlalchemy.create_engine('sqlite:///' + self.root_folder + db_name + self.db_type)
		connection = engine.connect()
		command1 = 'select ' + name1 + ' from ' + db_name
		command2 = 'select ' + name2 + ' from ' + db_name
		result1 = engine.execute(command1)
		result2 = engine.execute(command2)
		list1 = []
		list2 = []
		for row in result1:
			list1.append(row[name1.upper()])
		for row in result2:
			list2.append(row[name2.upper()])
		self.assertEqual(list1, ['A0A2U3CKK2','A0A2U3CKK9', 'A0A2U3CKI6','A0A2U3CKJ5','A0A2U3CKL8','A0A2U3CKJ8'])
		self.assertEqual(list2, ["PWI49548.1",
								"PWI49549.1",
								"PWI49550.1",
								"PWI49551.1",
								"PWI49552.1",
								"PWI49553.1"])
