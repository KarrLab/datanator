from datanator.data_source import ymdb
import datetime
import dateutil
import os
import shutil
import tempfile
import unittest
	
class TestYmdbFromRemote(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		cls.cache_dirname = tempfile.mkdtemp()
		src = ymdb.ymdb(cache_dirname=cls.cache_dirname, download_backups=False, load_content=False, verbose=True, max_entries=12)
		src.load_content()
		src.session.close()
		src.engine.dispose()

	@classmethod
	def tearDownClass(cls):
		shutil.rmtree(cls.cache_dirname)

	def setUp(self):
		self.src = ecmdb.Ecmdb(cache_dirname=self.cache_dirname, download_backups=False, load_content=False, verbose=True, max_entries=12)

	def tearDown(self):
		src = self.src
		src.session.close()
		src.engine.dispose()
		

	def test_load_content(self):
		src = self.src
		session = src.session
		q = session.query(ymdb.ymdb_id)				