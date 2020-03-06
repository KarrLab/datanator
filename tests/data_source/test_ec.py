import unittest
from datanator.data_source import ec
import datanator.config.core
import shutil
import tempfile
from pathlib import Path


class TestEC(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dir = tempfile.mkdtemp()
        db = 'test'
        username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
        cls.src = ec.EC(server=MongoDB, db=db, username=username, password=password, authSource='admin',
                        readPreference='nearest', max_entries=20, cache_dir=cls.cache_dir)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dir)
        cls.src.db.drop_collection(cls.src.collection_str)
        cls.src.client.close()

    @unittest.skip('IP')
    def test_establish_ftp(self):
        ftp = self.src.establish_ftp()
        self.assertTrue('enzyme.dat' in ftp.nlst())
        
    @unittest.skip('IP')
    def test_retrieve_content(self):
        p = Path(self.cache_dir+'/enzyme.dat')
        self.src.retrieve_content()
        self.assertTrue(p.exists())

    # @unittest.skip('for now')
    def test_parse_content(self):
        location = str(Path('~/karr_lab/datanator/docs/enzyme.dat').expanduser())
        self.src.parse_content(location)

    def test_make_doc(self):
        lines = ["ID   1.1.1.1", "DE   Alcohol dehydrogenase.", "AN   Aldehyde reductase.",
        "CA   (1) A primary alcohol + NAD(+) = an aldehyde + NADH.", "CA   (2) A secondary alcohol + NAD(+) = a ketone + NADH.",
        "CF   Zn(2+) or Fe cation."]
        result = self.src.make_doc(lines)
        self.assertEqual(result, {'ec_number': '1.1.1.1', 'ec_name': 'Alcohol dehydrogenase', 
                                  'ec_synonyms': ['Aldehyde reductase'], 
                                  'catalytic_activity': ['(1) A primary alcohol + NAD(+) = an aldehyde + NADH', '(2) A secondary alcohol + NAD(+) = a ketone + NADH'], 
                                  'cofactor': 'Zn(2+) or Fe cation'})