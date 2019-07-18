from datanator.data_source import sabio_rk
import datanator.config.core
import unittest
import tempfile
import shutil
import requests
import libsbml

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
        cls.src = sabio_rk.SabioRk(cache_dirname=cls.cache_dirname,
                                         MongoDB=MongoDB,  db=db,
                                         verbose=True, max_entries=20, username=username,
                                         password=password, webservice_batch_size = 10)
        sbml = requests.get('http://sabiork.h-its.org/sabioRestWebServices/kineticLaws', params={
                'kinlawids': ','.join('4096'),
            }).text
        reader = libsbml.SBMLReader()
        doc = reader.readSBMLFromString(sbml)
        cls.test_model = doc.getModel()
        cls.species_sbml = cls.test_model.getListOfSpecies()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.src.client.close()

    @unittest.skip('passed')
    def test_load_kinetic_law_ids(self):
        ids = self.src.load_kinetic_law_ids()
        self.assertEqual(ids[0:10], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        self.assertGreater(len(ids), 55000)

    @unittest.skip('passed')
    def test_create_cross_references_from_sbml(self):
        x_refs = self.src.create_cross_references_from_sbml(self.test_model.getListOfSpecies().get(0))
        exp = [{'namespace': 'chebi', 'id': 'CHEBI:16670'}, {'namespace': 'kegg.compound', 'id': 'C00012'}]
        self.assertEqual(exp, x_refs)

    def test_parse_enzyme_name(self):
        name, is_wildtype, variant = self.src.parse_enzyme_name(self.species_sbml.get(2).getName())
        self.assertEqual('S156E/S166D of subtilisin DSAI (N76D/N87S/S103A/V104I)', variant)
        self.assertEqual('subtilisin', name)

    def test_get_specie_from_sbml(self):
        specie, properties = self.src.get_specie_from_sbml(self.species_sbml.get(2))
        print(specie)
        print(properties)