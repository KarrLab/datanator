from datanator.data_source import sabio_rk
from datanator.util import file_util
import datanator.config.core
import unittest
import tempfile
import shutil
import requests
import libsbml
import bs4

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
        cls.file_manager = file_util.FileUtil()

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
        exp = [{'chebi': 'CHEBI:16810'}, {'chebi': 'CHEBI:30915'}, 
        {'kegg.compound': 'C00026'}]
        self.assertEqual(exp, x_refs)

    # @unittest.skip('passed')
    def test_parse_enzyme_name(self):
        name, is_wildtype, variant = self.src.parse_enzyme_name(self.species_sbml.get(5).getName())
        self.assertEqual('E211S/I50N/V80T', variant)
        self.assertEqual('4-aminobutyrate transaminase', name)

    # @unittest.skip('passed')
    def test_get_specie_from_sbml(self):
        specie, properties = self.src.get_specie_from_sbml(self.species_sbml.get(5))
        specie_exp = {'_id': 141214, 'molecular_weight': None, 'name': '4-aminobutyrate transaminase', 'subunits': [{'uniprot': 'P22256'}, {'uniprot': 'P50457'}], 
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
        self.assertEqual(specie[0]['subunits'], [{'uniprot': 'P22256'}, 
            {'uniprot': 'P50457'}])

    # @unittest.skip('passed')
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
        self.assertEqual(result['reactants'][0]['_id'], test_1)

    # @unittest.skip('passed')
    def test_create_kinetic_laws_from_sbml(self):
        ids = [4096]
        self.src.create_kinetic_laws_from_sbml(ids, self.sbml)
        doc = self.src.collection.find_one({'kinlaw_id':ids[0]})
        test_1 = doc.get('enzyme', None)
        self.assertEqual(test_1[0]['_id'], 141214)
        test_2 = doc.get('media', None)
        self.assertEqual(test_2, '25 mM Tris or MOPS, 10 mM MgCl2, 0.04 % sodium azide')

    # @unittest.skip('passed')
    def test_load_compounds(self):
        compound_1 = {
            "_id" : 1922,
            "name" : "2-Oxoglutarate",
            "cross_references" : [
                {
                    "namespace" : "chebi",
                    "id" : "CHEBI:16810"
                },
                {
                    "namespace" : "chebi",
                    "id" : "CHEBI:30915"
                },
                {
                    "namespace" : "kegg.compound",
                    "id" : "C00026"
                }
            ]
        }

        compound_2 = {
            "_id" : 21128,
            "name" : "2-Methylaspartic acid",
            "cross_references" : []
        }

        self.src.load_compounds(compounds = [compound_1, compound_2])
        test_1 = self.src.collection_compound.find_one({'_id': compound_1['_id']})
        test_2 = self.src.collection_compound.find_one({'_id': compound_2['_id']})
        self.assertTrue('synonyms' in test_1)
        self.assertTrue(isinstance(test_2['structures'], list))

    # @unittest.skip('passed')
    def test_get_parameter_by_properties(self):
        kinetic_law_mock = {'kinlaw_id': 4096, 'mechanism': 'mock_mechanism',
                            'tissue': 'mock_tissue', 'enzyme_type': 'mock_et',
                            'parameters': [{'observed_type': ['mock_ot', 'ssss'], 'compound': None,
                                        'observed_value': ['mock_ov', 'some_1']}]}
        parameter_properties_mock = {'type_code': ['mock_ot'], 'associatedSpecies': None,
                                'startValue': ['mock_ov', 'some_2'], 'type': 'some_type'}
        result = self.src.get_parameter_by_properties(kinetic_law_mock, parameter_properties_mock)
        exp = {'observed_type': ['mock_ot', 'ssss'], 'compound': None, 'observed_value': ['mock_ov', 'some_1']}
        self.assertEqual(result, exp)

    # @unittest.skip('passed')
    def test_load_missing_kinetic_law_information_from_tsv_helper(self):
        url = 'http://sabiork.h-its.org/entry/exportToExcelCustomizable'
        response = requests.get(url, params={
            'entryIDs[]': [4096],
            'fields[]': [
                'EntryID',
                'KineticMechanismType',
                'Tissue',
                'Parameter',
            ],
            'preview': False,
            'format': 'tsv',
            'distinctRows': 'false',
        })
        tsv = response.text
        self.src.load_missing_kinetic_law_information_from_tsv_helper(tsv)
        result = self.src.collection.find_one({'kinlaw_id': 4096})
        self.assertEqual(result.get('mechanism', 'no mechanism filed'), 'no mechanism filed')

    # @unittest.skip('passed')
    def test_infer_compound_structures_from_names(self):
        compound_1 = {
            "_id" : 73,
            "name" : "L-Glutamate",
            "cross_references" : [
                {
                    "namespace" : "chebi",
                    "id" : "CHEBI:16015"
                },
                {
                    "namespace" : "chebi",
                    "id" : "CHEBI:29972"
                },
                {
                    "namespace" : "chebi",
                    "id" : "CHEBI:29985"
                },
                {
                    "namespace" : "chebi",
                    "id" : "CHEBI:29988"
                },
                {
                    "namespace" : "kegg.compound",
                    "id" : "C00025"
                }
            ]
        }
        compound_2 = {
            "_id" : 1922,
            "name" : "2-Oxoglutarate",
            "cross_references" : [
                {
                    "namespace" : "kegg.compound",
                    "id" : "C00026"
                },
                {
                    "namespace" : "pubchem.substance",
                    "id" : "3328"
                },
                {
                    "namespace" : "chebi",
                    "id" : "CHEBI:16810"
                },
                {
                    "namespace" : "chebi",
                    "id" : "CHEBI:30915"
                },
                {
                    "namespace" : "reactome",
                    "id" : "113671"
                },
                {
                    "namespace" : "biocyc",
                    "id" : "2-KETOGLUTARATE"
                },
                {
                    "namespace" : "metanetx.chemical",
                    "id" : "MNXM20"
                },
                {
                    "namespace" : "BioModels",
                    "id" : "16810"
                },
                {
                    "namespace" : "BioModels",
                    "id" : "30915"
                }
            ],
            "structures" : [
                {
                    "inchi" : "InChI=1S/C5H6O5/c6-3(5(9)10)1-2-4(7)8/h1-2H2,(H,7,8)(H,9,10)"
                },
                {
                    "smiles" : "OC(=O)CCC(=O)C(O)=O"
                }
            ]
        }
        result = self.src.infer_compound_structures_from_names([compound_1, compound_2])
        self.assertEqual(result[1], compound_2)
        self.assertTrue('structures' in result[0])

    # @unittest.skip('passed')
    def test_calc_inchi_formula_connectivity(self):
        s = {'smiles': '[H]O[H]'}
        test_1 = self.src.calc_inchi_formula_connectivity(s)
        self.assertEqual(test_1['_value_inchi'], 'InChI=1S/H2O/h1H2')
        self.assertEqual(test_1['_value_inchi_formula_connectivity'], 'H2O')

        s = {'inchi': 'InChI=1S/H2O/h1H2'}
        test_2 = self.src.calc_inchi_formula_connectivity(s)
        self.assertEqual(test_2['_value_inchi'], 'InChI=1S/H2O/h1H2')
        self.assertEqual(test_2['_value_inchi_formula_connectivity'], 'H2O')

        s = {'inchi': 'InChI=1S/C9H10O3/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8,10H,6H2,(H,11,12)/t8-/m1/s1'}
        test_3 = self.src.calc_inchi_formula_connectivity(s)
        self.assertEqual(test_3['_value_inchi'], 'InChI=1S/C9H10O3/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8,10H,6H2,(H,11,12)/t8-/m1/s1')
        self.assertEqual(test_3['_value_inchi_formula_connectivity'], 'C9H10O3/c10-8(9(11)12)6-7-4-2-1-3-5-7')

    # @unittest.skip('passed')
    def test_parse_complex_subunit_structure(self):
        response = requests.get('http://sabiork.h-its.org/kindatadirectiframe.jsp', params={
                                'kinlawid': 4096, 'newinterface': 'true'})
        doc = bs4.BeautifulSoup(response.text, 'html.parser')
        td = doc.find('td', text='Modifier-Catalyst')
        tr = td.parent
        td = tr.find_all('td')[-1]
        inner_html = td.decode_contents(formatter='html').strip() + ' '
        test_1 = self.src.parse_complex_subunit_structure(inner_html)
        self.assertEqual({'P22256': 2, 'P50457': 1}, test_1)

    # @unittest.skip('passed')
    def test_load_missing_enzyme_information_from_html(self):
        ids = [4096]
        self.src.load_missing_enzyme_information_from_html(ids)
        projection = {'enzyme':1}
        test_doc = self.src.collection.find_one(filter={'kinlaw_id': { '$in': ids }}, projection=projection)
        l = self.file_manager.search_dict_list(test_doc['enzyme'], 'coeffcient')
        self.assertFalse(len(l)>0)

    # @unittest.skip('passed')
    def test_calc_enzyme_molecular_weights(self):
        null = None
        enzyme = [{
            "_id" : 141214,
            "molecular_weight" : null,
            "name" : "4-aminobutyrate transaminase",
            "subunits" : [
                {
                    "uniprot" : "P22256",
                    "coefficient" : 2
                },
                {
                    "uniprot" : "P50457",
                    "coefficient" : 1
                }
            ],
            "cross_references" : [ ]
        }]
        results = self.src.calc_enzyme_molecular_weights(enzyme, len(enzyme))
        self.assertTrue(results[0]['molecular_weight'] != None)

    # @unittest.skip('takes too long')
    def test_load_content(self):
        self.src.load_content()

    def test_add_inchi_hash(self):
        result = self.src.add_inchi_hash()
        query = {}
        projection = {'products': 1, 'kinlaw_id': 1, 'reactants': 1, '_id': 0}
        col = self.src.collection.find(filter=query, projection=projection)
        for doc in col:
            self.assertTrue('InChI_Key' in doc['products'][0]['structures'][0] )