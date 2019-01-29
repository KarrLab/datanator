from datanator.data_source import uniprot_ncbi
import datetime
import dateutil
import os
import shutil
import tempfile
import unittest
from pathlib import Path

warning_util.disable_warnings()

class TestDownloader(unittest.TestCase):

    def setUp(self):
    	src = uniprot_ncbi
        self.root_folder = tempfile.mkdtemp()
        self.file_name = ["oma-uniprot", "oma-ncbi"]
        self.file_type = ".txt.gz"
        self.db_type = '.sqlite'
		src.download_data(self.file_name[0],self.file_type,self.root_folder)
		src.download_data(self.file_name[1],self.file_type,self.root_folder)

    def tearDown(self):
        shutil.rmtree(self.root_folder)

    def test_download_data(self):
    	src = self.src
    	storage_location1 = self.root_folder + self.file_name[0] + self.file_type
    	storage_location2 = self.root_folder + self.file_name[1] + self.file_type
    	p = Path(storage_location1)
    	q = Path(storage_location2)
    	self.assertTrue(p.exists())
    	self.assertTrue(q.exists())

    def test_load_data(self):
    	src = self.src
    	df = src.load_data(self.file_name[0],self.file_type,self.root_folder)
    	temp = df[4:10]
    	col0 = file_name[0].split('-')[0]
    	col1 = file_name[0].split('-')[1]
    	self.assertEqual(temp[col0].tolist(), ['HEIAB00001', 'HEIAB00002', 'HEIAB00003',
    												'HEIAB00004','HEIAB00005','HEIAB00006','HEIAB00007'])
    	self.assertEqual(temp[col1].tolist(), ['A0A2U3CKJ0', 'A0A2U3CKK2','A0A2U3CKK9',
    											'A0A2U3CKI6','A0A2U3CKJ5','A0A2U3CKL8','A0A2U3CKJ8'])
    def test_sql_data(self):
    	src = self.src
    	df1 = src.load_data(self.file_name[0],self.file_type,self.root_folder)
    	df2 = src.load_data(self.file_name[1],self.file_type,self.root_folder)
    	temp = df[4:10]
    	sql_data