from datanator.data_source import sabio_rk_trial
import datanator.config.core
import unittest
import tempfile
import shutil


class TestSabioRk(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        db = 'test'
        username = datanator.config.core.get_config()[
            'datanator']['mongodb']['user']
        password = datanator.config.core.get_config(
        )['datanator']['mongodb']['password']
        MongoDB = datanator.config.core.get_config(
        )['datanator']['mongodb']['server']
        port = datanator.config.core.get_config(
        )['datanator']['mongodb']['port']
        replSet = datanator.config.core.get_config(
        )['datanator']['mongodb']['replSet']
        cls.src = sabio_rk_trial.SabioRk(cache_dirname=cls.cache_dirname,
                                         MongoDB=MongoDB,  db=db,
                                         verbose=True, max_entries=20, username=username,
                                         password=password, webservice_batch_size = 10)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.src.client.close()

    def test_load_kinetic_law_ids(self):
        ids = self.src.load_kinetic_law_ids()
        self.assertEqual(ids[0:10], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        self.assertGreater(len(ids), 55000)