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
        cls.sbml = requests.get('http://sabiork.h-its.org/sabioRestWebServices/kineticLaws', params={
                'kinlawids': '4096'}).text
        cls.reader = libsbml.SBMLReader()
        cls.doc = cls.reader.readSBMLFromString(cls.sbml)
        cls.test_model = cls.doc.getModel()
        cls.species_sbml = cls.test_model.getListOfSpecies()
        cls.reactions_sbml = cls.test_model.getListOfReactions()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)
        cls.src.client.close()

    @unittest.skip('passed, avoid unnecessary http requests')
    def test_load_kinetic_law_ids(self):
        ids = self.src.load_kinetic_law_ids()
        self.assertEqual(ids[0:10], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        self.assertGreater(len(ids), 55000)

    # @unittest.skip('passed')
    def test_create_cross_references_from_sbml(self):
        x_refs = self.src.create_cross_references_from_sbml(self.species_sbml.get(0))
        exp = [{'namespace': 'chebi', 'id': 'CHEBI:16810'}, {'namespace': 'chebi', 'id': 'CHEBI:30915'}, 
        {'namespace': 'kegg.compound', 'id': 'C00026'}]
        self.assertEqual(exp, x_refs)

    def test_parse_enzyme_name(self):
        name, is_wildtype, variant = self.src.parse_enzyme_name(self.species_sbml.get(5).getName())
        self.assertEqual('E211S/I50N/V80T', variant)
        self.assertEqual('4-aminobutyrate transaminase', name)

    # @unittest.skip('passed')
    def test_get_specie_from_sbml(self):
        specie, properties = self.src.get_specie_from_sbml(self.species_sbml.get(5))
        specie_exp = {'_id': 141214, 'molecular_weight': None, 'name': '4-aminobutyrate transaminase', 'subunits': [{'namespace': 'uniprot', 'id': 'P22256'}, {'namespace': 'uniprot', 'id': 'P50457'}], 
        'cross_references': []}
        properties_exp = {'is_wildtype': False, 'variant': 'E211S/I50N/V80T', 'modifier_type': 'Modifier-Catalyst'}
        self.assertEqual(specie['_id'], specie_exp['_id'])
        self.assertEqual(properties_exp['variant'], properties['variant'])

    # @unittest.skip('passed')
    def test_get_specie_reference_from_sbml(self):
        species = []
        for i_specie in range(self.species_sbml.size()):
            specie_sbml = self.species_sbml.get(i_specie)
            specie, properties = self.src.get_specie_from_sbml(specie_sbml)
            species.append(specie)
        specie, compartment = self.src.get_specie_reference_from_sbml('ENZ_141214_Cell', species)
        self.assertEqual(compartment, None)
        self.assertEqual(specie[0]['subunits'], [{'namespace': 'uniprot', 'id': 'P22256'}, 
            {'namespace': 'uniprot', 'id': 'P50457'}])

    def test_create_kinetic_law_from_sbml(self):
        species = []
        specie_properties = {}
        for i_specie in range(self.species_sbml.size()):
            specie_sbml = self.species_sbml.get(i_specie)
            specie, properties = self.src.get_specie_from_sbml(specie_sbml)
            species.append(specie)
            specie_properties[specie_sbml.getId()] = properties
        units = {}
        units_sbml = self.test_model.getListOfUnitDefinitions()
        for i_unit in range(units_sbml.size()):
            unit_sbml = units_sbml.get(i_unit)
            units[unit_sbml.getId()] = unit_sbml.getName()

        functions = {}
        functions_sbml = self.test_model.getListOfFunctionDefinitions()
        for i_function in range(functions_sbml.size()):
            function_sbml = functions_sbml.get(i_function)
            math_sbml = function_sbml.getMath()
            if math_sbml.isLambda() and math_sbml.getNumChildren():
                eq = libsbml.formulaToL3String(math_sbml.getChild(math_sbml.getNumChildren() - 1))
            else:
                eq = None
            if eq in ('', 'NaN'):
                eq = None
            functions[function_sbml.getId()] = eq

        result = self.src.create_kinetic_law_from_sbml(4096, self.reactions_sbml.get(0), species, 
                                                        specie_properties, functions, units)
        test_1 = 1922
        self.assertEqual(result['reactants'][0]['compound'][0]['_id'], test_1)

    def test_create_kinetic_laws_from_sbml(self):
        ids = [4096]
        self.src.create_kinetic_laws_from_sbml(ids, self.sbml)
        doc = self.src.collection.find_one({'kinlaw_id':ids[0]})
        test_1 = doc['compartments'][0]
        self.assertEqual(test_1, None)
        test_2 = doc['species'][0]['_id']
        self.assertEqual(test_2, 1922)