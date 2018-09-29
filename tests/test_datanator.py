""" Tests of kinetic_datanator

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-06
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator import datanator
from kinetic_datanator.util import warning_util
from os import path
import os
import unittest

warning_util.disable_warnings()


class TestDatanator(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.fixtures_dir = path.join(path.dirname(__file__), "fixtures")

        out_dir = cls.output_dir = path.join(path.dirname(__file__), "output")
        if not path.isdir(out_dir):
            os.makedirs(out_dir)

    def test_init(self):
        d = datanator.Datanator(max_taxon_dist=10, include_variants=True,
                                min_temp=34, max_temp=40, min_ph=6, max_ph=8)
        self.assertEqual(d.max_taxon_dist, 10)
        self.assertEqual(d.include_variants, True)
        self.assertEqual(d.min_temp, 34)
        self.assertEqual(d.max_temp, 40)
        self.assertEqual(d.min_ph, 6)
        self.assertEqual(d.max_ph, 8)

    def test_read_model(self):
        i_file = path.join(self.fixtures_dir, "ump_kinase.xlsx")
        d = datanator.Datanator()
        model = d.read_model(i_file)
        taxon, compartments, molecules, reactions = model
        self.assertEqual(taxon.taxon, 'Mycoplasma pneumoniae M129')

    @unittest.skip('skipping/unclear how relevant')
    def test_annotate_model(self):
        i_file = path.join(self.fixtures_dir, "ump_kinase.xlsx")
        d = datanator.Datanator()
        model = d.read_model(i_file)
        d.annotate_model(model)

        taxon, compartments, molecules, reactions = model

        mol = next(mol for mol in molecules if mol.id == 'UMP')
        ids = set(xr.id for xr in mol.cross_references if xr.namespace == 'sabio-id')
        names = set(xr.id for xr in mol.cross_references if xr.namespace == 'sabio-name')
        print(ids)
        self.assertEqual(ids, set([7458, 1300, 25061]))
        self.assertEqual(names, set(["2,4-Dioxotetrahydropyrimidine D-ribonucleotide", "UMP", "Uridine 5'-phosphate"]))

        rxn = next(rxn for rxn in reactions if rxn.id == 'ump_kinase')
        self.assertEqual(rxn.get_ec_number(), '2.7.4')  # true EC is 2.7.4.22

    @unittest.skip('skipping/unclear how relevant')
    def test_annotate_molecules(self):
        i_file = path.join(self.fixtures_dir, "ump_kinase.xlsx")
        d = datanator.Datanator()
        model = d.read_model(i_file)
        d.annotate_molecules(model)

        taxon, compartments, molecules, reactions = model

        mol = next(mol for mol in molecules if mol.id == 'UMP')
        ids = set(xr.id for xr in mol.cross_references if xr.namespace == 'sabio-id')
        names = set(xr.id for xr in mol.cross_references if xr.namespace == 'sabio-name')
        self.assertEqual(ids, set([7458, 1300, 25061]))
        self.assertEqual(names, set(["2,4-Dioxotetrahydropyrimidine D-ribonucleotide", "UMP", "Uridine 5'-phosphate"]))

    @unittest.skip('skipping/unclear how relevant')
    def test_annotate_reactions(self):
        i_file = path.join(self.fixtures_dir, "ump_kinase.xlsx")
        d = datanator.Datanator()
        model = d.read_model(i_file)
        d.annotate_reactions(model)

        taxon, compartments, molecules, reactions = model
        rxn = next(rxn for rxn in reactions if rxn.id == 'ump_kinase')
        self.assertEqual(rxn.get_ec_number(), '2.7.4')  # true EC is 2.7.4.22

    @unittest.skip('skipping/unclear how relevant')
    def test_datanator(self):
        # Test getting kinetic data about a reaction that is relevant for a species

        i_file = path.join(self.fixtures_dir, "ump_kinase.xlsx")
        o_file = path.join(self.output_dir, "ump_kinase.xlsx")

        # get kinetic data for reations in Excel sheet
        rxns = datanator.get_sabio_data(i_file, o_file)

        # test correct number of output reactions
        self.assertEqual(len(rxns), 1)

        # test outputted reaction is UMP kinases
        rxn = rxns[0]
        self.assertEqual(rxn.id, 'ump_kinase')

        # test the formatted data fields
        self.assertEqual(rxn.reaction_ids, ['201'])

        # median_km_entry is an object of the Entry class found in the sabio_rk module
        # the following are a series of tests that make sure that the entry fields are
        #filled in properly
        median_entry = rxn.km_data.median_entry
        self.assertNotEqual(median_entry, None)
        self.assertEqual(median_entry.entry_id, '17927')
        self.assertEqual(median_entry.vmax, '0.00665')
        self.assertEqual(median_entry.proximity, 6)


@unittest.skip('skipping/unclear how relevant')
class TestProgram(unittest.TestCase):

    def test_datanator(self):
        # Find Excel sheet with reaction data
        inputFileName = path.join(path.dirname(__file__), "fixtures", "five_reactions.xlsx")
        # turn Excel sheet into openpyxl workbook
        if not path.isdir(path.join(path.dirname(__file__), "output")):
            os.makedirs(path.join(path.dirname(__file__), "output"))
        outputFilename = path.join(path.dirname(__file__), "output", "five_reactions.xlsx")
        species = 'mycoplasma pneumoniae'

        # formatted data list
        # The formatted data list is a list of SabioResult objects
        # SabioResult objects store and organize the final relevant expiremental data
        # SabioResult is the central Class that should contain all the data needed to respond to a query
        # Ideally, a SabioResult objective should be comprehensive enough to pass to a user, and the user should be able
        # to derive all the information he wants from that object

        rxns = datanator.get_kinetic_data(inputFileName, outputFilename, species)
        self.assertEqual(set([rxn.id for rxn in rxns]), set([
            'ump_kinase',
            'gmp_kinase',
            'oligonucleotidase_dcmp_dtmp',
            'fmn_reductase',
            'nucleotidase_7gmp',
        ]))

        ##################################################
        ##################################################
        """ UMP kinase """
        rxn = next((rxn for rxn in rxns if rxn.id == 'ump_kinase'), None)
        self.assertNotEqual(rxn, None)

        # reaction_ids is there to classify how many distinct sabio-RK Ids reference the desired reaction
        # there are two reasons a reaction ID list may contain more than a single entry. Either if there is an error
        # or if there are two reactions that differ only in stereochemistry (the inchi strings and smiles cannot detect
        # stereochemical distinctions).
        # the reaction IDs for this reaction should just be 201
        self.assertEqual(rxn.reaction_ids, ['201'])
        # both km and vmax should have KineticInfo objects in their fields
        self.assertNotEqual(rxn.km_data, None)
        self.assertNotEqual(rxn.vmax_data, None)

        # Now the KineticInfo fields will be tested
        # First is some summaries of the expiremental data.
        km_data = rxn.km_data
        # it has two EC numbers for the same reaction. EC numbers are not a perfect system
        self.assertEqual(km_data.ec_numbers, ['2.7.4.14', '2.7.4.22'])
        self.assertEqual(km_data.sabio_reaction_ids, ['201'])
        self.assertEqual(km_data.lift_info, "Lift Not Used")
        self.assertEqual(km_data.reaction_list, ['ATP + UMP = UDP + ADP'])

        # Each kinetic info class contains Entry objects from the sabio_rk module.
        # The following tests will tests ensure the KineticInfo fields contains the right Entry objects

        closest_entry_ids = km_data.closest_entry_ids
        self.assertTrue(len(closest_entry_ids) > 12 and len(closest_entry_ids) < 25)
        closest_entries = km_data.closest_entries  # this is a list of the most closely related expiremental entries
        # when I made this test, there were 13 entries. If there are more than 25 entries, something probably went wrong
        self.assertTrue(len(closest_entries) > 12 and len(closest_entries) < 25)

        # from the set of closest_entries, the entry with median km value is the median_entry
        # the median_entry is an Entry object defined in sabio_rk
        median_entry = km_data.median_entry
        self.assertEqual(median_entry.ec_number, '2.7.4.22')
        self.assertEqual(median_entry.entry_id, '17927')  # this is the specific Entry Number sabio assigns to each expiremental entry
        self.assertEqual(median_entry.km, '1.0E-4')  # this is the km value for this expiremental entry
        self.assertEqual(median_entry.vmax, '0.00665')  # this is the vmax for this expiremental entry
        self.assertEqual(median_entry.num_participants, [2, 2])  # this is a list to record the number of substrates and products
        self.assertEqual(median_entry.species, 'Streptococcus pneumoniae')  # this is the species the expirement was done in
        self.assertEqual(median_entry.proximity, 6)  # this is closely related the expiremental species is to the modeler's species
        self.assertEqual(median_entry.reaction_id, '201')  # this is the Sabio assigned ID for this reaction in general

        # check the entry_id of the min and the max within km_data
        self.assertEqual(km_data.min_entry.entry_id, '51904')
        self.assertEqual(km_data.max_entry.entry_id, '51903')

        # now the vmax infor fields will  be tested
        vmaxData = rxn.vmax_data
        # it has two EC numbers for the same reaction. EC numbers are not a perfect system
        self.assertEqual(vmaxData.ec_numbers, ['2.7.4.14', '2.7.4.22'])
        self.assertEqual(vmaxData.sabio_reaction_ids, ['201'])
        self.assertEqual(vmaxData.lift_info, "Lift Not Used")
        self.assertEqual(vmaxData.reaction_list, ['ATP + UMP = UDP + ADP'])

        # Each kinetic info class contains Entry objects from the sabio_rk module.
        # The following tests will tests ensure the KineticInfo fields contains the right Entry objects

        closest_entry_ids = vmaxData.closest_entry_ids
        self.assertTrue(len(closest_entry_ids) > 44 and len(closest_entry_ids) < 50)
        closest_entries = vmaxData.closest_entries  # this is a list of the most closely related expiremental entries
        # when I made this test, there were 45 entries. If there are more than 50 entries, something probably went wrong
        self.assertTrue(len(closest_entries) > 44 and len(closest_entries) < 50)

        # from the set of closest_entries, the entry with median km value is the median_entry
        # the median_entry is an Entry object defined in sabio_rk
        median_entry = vmaxData.median_entry
        self.assertEqual(median_entry.ec_number, '2.7.4.22')
        self.assertEqual(median_entry.entry_id, '17907')  # this is the specific Entry Number sabio assigns to each expiremental entry
        self.assertEqual(median_entry.km, '')  # this is the km value for this expiremental entry
        self.assertEqual(median_entry.vmax, '0.00456666667')  # this is the vmax for this expiremental entry
        self.assertEqual(median_entry.num_participants, [2, 2])  # this is a list to record the number of substrates and products
        self.assertEqual(median_entry.species, 'Streptococcus pneumoniae')  # this is the species the expirement was done in
        self.assertEqual(median_entry.proximity, 6)  # this is closely related the expiremental species is to the modeler's species
        self.assertEqual(median_entry.reaction_id, '201')  # this is the Sabio assigned ID for this reaction in general

        # check the entry_id of the min and the max within km_data
        self.assertEqual(vmaxData.min_entry.entry_id, '40434')
        self.assertEqual(vmaxData.max_entry.entry_id, '17934')

        ##################################################
        ##################################################
        rxn = next((rxn for rxn in rxns if rxn.id == 'oligonucleotidase_dcmp_dtmp'), None)
        self.assertNotEqual(rxn, None)

        # the reaction IDs for this reaction should just be empty becuase no exact mathches were found in the Sabio Database
        self.assertEqual(rxn.reaction_ids, [])
        # both km and vmax should have KineticInfo objects in their fields
        self.assertNotEqual(rxn.km_data, None)
        self.assertNotEqual(rxn.vmax_data, None)

        # Now the KineticInfo fields will be tested
        # First is some summaries of the expiremental data.
        km_data = rxn.km_data
        # No Sabio reaction should match the query, and therefore Sabio will gather data from any reaction in
        # the same EC classification subclass (the first three digits of the four digit EC number)
        self.assertEqual(km_data.lift_info, 'Lifted From 3.1.3')
        # The query gathers data from many different reactions, and therefore many the reactions will have different EC numbers.
        self.assertEqual(km_data.ec_numbers, ['3.1.3.1', '3.1.3.11', '3.1.3.13', '3.1.3.16', '3.1.3.17', '3.1.3.18', '3.1.3.2', '3.1.3.22', '3.1.3.23', '3.1.3.24', '3.1.3.25', '3.1.3.29', '3.1.3.3', '3.1.3.31', '3.1.3.33', '3.1.3.34', '3.1.3.36', '3.1.3.37', '3.1.3.4', '3.1.3.43',
                                              '3.1.3.46', '3.1.3.48', '3.1.3.5', '3.1.3.56', '3.1.3.57', '3.1.3.6', '3.1.3.64', '3.1.3.7', '3.1.3.74', '3.1.3.76', '3.1.3.8', '3.1.3.9', '3.1.3.91', '3.1.3.93', '3.1.3.95'])  # it has two EC numbers for the same reaction. EC numbers are not a perfect system
        # Similarly, the reactions gatheres will have different Sabio reaction IDs
        self.assertEqual(km_data.sabio_reaction_ids, ['10261', '10263', '10269', '10270', '10271', '10330', '10377', '10378', '10388', '10394', '10399', '10400', '10401', '10631', '10700', '10707', '10711', '10726', '10734', '1117', '1118', '11265', '11266', '11267', '11288', '11318', '1150', '1182', '12091', '12166', '124', '12404', '12405', '129', '13054', '13055', '1325', '13296', '13453', '13635', '13655', '13970', '14184', '14198', '14219', '1422', '1423', '1424', '1683', '1686', '1851', '186', '1891', '1923', '1965', '200', '207', '2084', '209', '210', '211', '2219', '2424', '246', '2493', '2503', '27', '2717', '2730', '282', '2874', '295', '304', '305', '309', '3140', '3170',
                                                      '3177', '3197', '3198', '3224', '3227', '339', '4095', '433', '444', '479', '480', '487', '501', '508', '581', '6096', '6098', '6445', '6446', '664', '6734', '6769', '697', '7185', '7186', '7187', '7188', '7189', '7190', '7191', '7192', '7193', '7194', '7195', '7196', '7278', '728', '736', '742', '75', '76', '7712', '7713', '7952', '7953', '7954', '796', '797', '7970', '7971', '7972', '7973', '7974', '7975', '7976', '7977', '8072', '8082', '816', '8203', '8204', '8205', '8206', '8749', '8902', '8953', '9216', '937', '9556', '9571', '9583', '9586', '9599', '9603', '964', '983', '9864', '9910', '9945', '9953', '9954', '9955', '9956', '9957', '9958', '9959', '9960', '9963'])
        self.assertEqual(km_data.reaction_list, ["1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3',4',5'-trisphosphate) + H2O = 1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3',4'-bisphosphate) + Phosphate", "1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3',5'-bisphosphate) + H2O = 1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3'-phosphate) + Phosphate", "1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 4',5'-bisphosphate) + H2O = 1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 4'-phosphate) + Phosphate", '1-Phosphatidyl-D-myo-inositol 4,5-bisphosphate + H2O = Phosphate + 1-Phosphatidyl-1D-myo-inositol 4-phosphate', '1D-myo-Inositol 1,4,5,6-tetrakisphosphate + H2O = Inositol 1,4,6-trisphosphate + Phosphate', '1D-myo-Inositol 2-phosphate + H2O = myo-Inositol + Phosphate', "2'-Deoxyadenosine 3'-phosphate + H2O = Deoxyadenosine + Phosphate", "2'-Deoxycytidine 3'-phosphate + H2O = 2'-Deoxycytidine + Phosphate", "2'-Deoxyguanosine 3'-phosphate + H2O = Deoxyguanosine + Phosphate", "2'-Deoxyuridine 3'-phosphate + H2O = 2'-Deoxyuridine + Phosphate", '2,3-Diphosphoglycerate + H2O = Unknown + Phosphate', '2-Deoxy-D-glucose 6-phosphate + H2O = 2-Deoxy-D-glucose + Phosphate', '2-Deoxyglucose 6-phosphate + H2O = 2-Deoxyglucose + Phosphate', "3'-Phosphoadenylyl sulfate + H2O = Phosphate + Adenylylsulfate", '4-Pyridoxic acid 5-phosphate + H2O = 4-Pyridoxate + Phosphate', "5'-Ribonucleotide + H2O = Nucleoside + Phosphate", 'ADP + H2O = Diphosphate + Adenosine', 'ADP + H2O = Phosphate + AMP', 'ATP + H2O = Adenosine + Triphosphate', 'ATP + H2O = Phosphate + ADP', "Adenosine 3',5'-bisphosphate + H2O = Phosphate + AMP", 'Arabinose 5-phosphate + H2O = Arabinose + Phosphate', 'CDP + H2O = CMP + Phosphate', "Cytidine 3'-phosphate + H2O = Cytidine + Phosphate", 'D-Fructose 1,6-bisphosphate + H2O = D-Fructose monophosphate + Phosphate', 'D-Glucose 1-phosphate + H2O = D-Glucose + Phosphate', 'D-Glucose 6-phosphate + H2O = D-Glucose + Phosphate', 'D-Mannitol 1-phosphate + H2O = Phosphate + D-Mannitol', 'D-myo-Inositol 1,2,4,5,6-pentakisphosphate + H2O = Inositol 1,2,4,6-tetrakisphosphate + Phosphate', 'Diphosphate + H2O = Phosphate', 'Fructose 1,6-bisphosphate + H2O = Unknown + Phosphate', 'Fructose 1-phosphate + H2O = Fructose + Phosphate', 'Fructose 6-phosphate + H2O = Fructose + Phosphate', 'GDP + H2O = Phosphate + GMP', 'Galactose 1-phosphate + H2O = Galactose + Phosphate', 'Glucose + 3-Phospho-D-glyceroyl phosphate = Glucose 6-phosphate + 3-Phospho-D-glycerate', 'Glucose 1-phosphate + H2O = Glucose + Phosphate', 'Glycerate 3-phosphate + H2O = Glycerate + Phosphate', 'Glycerate 3-phosphate + H2O = Phosphate + D-Glycerate', 'Glycerol 3-phosphate + H2O = Glycerol + Phosphate', 'H2O + 1-Phosphatidyl-1D-myo-inositol 3,5-bisphosphate = Phosphate + 1-Phosphatidyl-1D-myo-inositol 5-phosphate', 'H2O + 1-Phosphatidyl-1D-myo-inositol 3-phosphate = Phosphate + 1-Phosphatidyl-D-myo-inositol', 'H2O + 1-Phospho-D-glycerate = D-Glycerate + Phosphate', 'H2O + 1D-myo-Inositol 1,3,4,5-tetrakisphosphate = Phosphate + 1D-myo-Inositol 1,3,4-trisphosphate', 'H2O + 1D-myo-Inositol 1,3,4-trisphosphate = Phosphate + D-myo-Inositol 3,4-bisphosphate', 'H2O + 1D-myo-Inositol 1,4-bisphosphate = Phosphate + myo-Inositol 4-phosphate', 'H2O + 1D-myo-Inositol 1-phosphate = Phosphate + myo-Inositol', 'H2O + 1D-myo-Inositol 3-phosphate = Phosphate + myo-Inositol', 'H2O + 2-Chloro-4-nitrophenyl phosphate = Phosphate + 2-Chloro-4-nitrophenol', 'H2O + 2-Phospho-D-glycerate = Phosphate + D-Glycerate', 'H2O + 2-Phosphoglycolate = Phosphate + Glycolate', 'H2O + 3-O-Methylfluorescein phosphate = Phosphate + 3-O-Methylfluorescein', 'H2O + 4-Chlorophenyl phosphate = Phosphate + 4-Chlorophenol', 'H2O + 4-Cyanophenyl phosphate = Phosphate + 4-Cyanophenol', 'H2O + 4-Methylumbelliferyl phosphate = Phosphate + 4-Methylumbelliferone', 'H2O + 4-Nitrophenyl phenyl phosphonate = p-Nitrophenol + Phenyl phosphonate', 'H2O + 4-Nitrophenyl phosphate = Phosphate + p-Nitrophenol', 'H2O + 4-Trifluoromethylphenyl phosphate = Phosphate + 4-Trifluoromethylphenol', "H2O + 5'-Phosphopolynucleotide = Phosphate + Polynucleotide", 'H2O + 5-Fluoro-4-methylumbelliferyl phosphate = Phosphate + 5-Fluoro-4-methylumbelliferone', 'H2O + 6,8-Difluoro-4-methylumbelliferyl phosphate = Phosphate + 6,8-Difluoro-4-methylumbelliferone', 'H2O + 6-Fluoro-4-methylumbelliferyl phosphate = Phosphate + 6-Fluoro-4-methylumbelliferone', 'H2O + 6-Phosphogluconate = Phosphate + Gluconate', "H2O + 7-Methylguanosine 5'-phosphate = Phosphate + 7-Methylguanosine", 'H2O + 8-Fluoro-4-methylumbelliferyl phosphate = Phosphate + 8-Fluoro-4-methylumbelliferone', 'H2O + AMP = Phosphate + Adenosine', "H2O + Adenosine 2'-phosphate = Adenosine + Phosphate", "H2O + Adenosine 3'-phosphate = Phosphate + Adenosine", 'H2O + Bis-4-nitrophenyl phosphate = p-Nitrophenol + 4-Nitrophenyl phosphate', 'H2O + CMP = Cytidine + Phosphate', 'H2O + CTP = CDP + Phosphate', 'H2O + Casein kinase I epsilon phosphorylated = Phosphate + Casein kinase I epsilon', 'H2O + Ceramide 1-phosphate = Phosphate + N-Acylsphingosine', 'H2O + Choline phosphate = Phosphate + Choline', 'H2O + D-Fructose 1,6-bisphosphate = Phosphate + D-Fructose 6-phosphate', 'H2O + D-Fructose 1,6-bisphosphate = Phosphate + beta-D-Fructose 6-phosphate',
                                                 'H2O + D-Fructose 2,6-bisphosphate = Phosphate + D-Fructose 6-phosphate', 'H2O + D-Fructose 2,6-bisphosphate = Phosphate + beta-D-Fructose 6-phosphate', 'H2O + D-Galactose 1-phosphate = Phosphate + D-Galactose', 'H2O + D-Mannose 6-phosphate = D-Mannose + Phosphate', 'H2O + D-O-Phosphoserine = Phosphate + D-Serine', 'H2O + D-myo-Inositol 1,3-bisphosphate = Phosphate + 1D-myo-Inositol 1-phosphate', 'H2O + D-myo-Inositol 1,4,5-trisphosphate = Phosphate + 1D-myo-Inositol 1,4-bisphosphate', "H2O + Deoxythymidine 3',5'-diphosphate = Phosphate + dTMP", "H2O + Deoxythymidine 3'-phosphate = Phosphate + Thymidine", 'H2O + Ethanolamine phosphate = Phosphate + Ethanolamine', 'H2O + Fructose 1,6-bisphosphate = Phosphate + Fructose 6-phosphate', 'H2O + GMP = Phosphate + Guanosine', 'H2O + GTP = GDP + Phosphate', 'H2O + Glucosamine 6-phosphate = Phosphate + Glucosamine', 'H2O + Glucose 6-phosphate = Glucose + Phosphate', 'H2O + Glycerate 2,3-bisphosphate = Glycerate 3-phosphate + Phosphate', 'H2O + Glycerate 2,3-bisphosphate = Phosphate + 3-Phospho-D-glycerate', 'H2O + Glycerol 2-phosphate = Phosphate + Glycerol', "H2O + Guanosine 3'-phosphate = Phosphate + Guanosine", 'H2O + IMP = Phosphate + Inosine', 'H2O + Inositol 1,2,3,4,5-pentakisphosphate = Phosphate + Inositol 1,2,3,4-tetrakisphosphate', 'H2O + Inositol 1,2,4,5-tetrakisphosphate = Inositol 1,2,4-trisphosphate + Phosphate', 'H2O + Inositol 4,5-bisphosphate = Phosphate + myo-Inositol 4-phosphate', 'H2O + L-Galactose 1-phosphate = Phosphate + L-Galactose', 'H2O + L-Phosphotyrosine = Phosphate + L-Tyrosine', 'H2O + L-Threonine O-3-phosphate = L-Threonine + Phosphate', 'H2O + N-(5-Phospho-4-pyridoxyl)glycine = 4-Pyridoxylglycine + Phosphate', 'H2O + N-Acetylneuraminate 9-phosphate = Phosphate + N-Acetylneuraminate', 'H2O + N-Acylneuraminate 9-phosphate = Phosphate + N-Acylneuraminate', 'H2O + NADP+ = Phosphate + NAD+', 'H2O + NADPH = NADH + Phosphate', 'H2O + O-Phospho-L-serine = L-Serine + Phosphate', 'H2O + O-Phospho-L-serine = Serine + Phosphate', 'H2O + O-Phospho-tau-protein = Phosphate + tau-Protein', 'H2O + Phenolic phosphate = Phosphate + Phenol', 'H2O + Phenolphthalein diphosphate = Phosphate + Phenolphthalein', 'H2O + Phosphoenolpyruvate = Phosphate + Pyruvate', 'H2O + Phosphotyrosine = Phosphate + Tyrosine', 'H2O + Propan-1-ol 2-phosphate = Phosphate + 1-Propanol', 'H2O + Pyridoxal phosphate = Phosphate + Pyridoxal', 'H2O + Ribulose 5-phosphate = Ribulose + Phosphate', 'H2O + Sorbitol 6-phosphate = Phosphate + Sorbitol', 'H2O + Sphinganine 1-phosphate = Phosphate + Sphinganine', 'H2O + SpoIIAA-phosphorylated = Phosphate + SpoIIAA', 'H2O + Sucrose 6-phosphate = Sucrose + Phosphate', 'H2O + TDP = Phosphate + Thiamine monophosphate', 'H2O + Tetrapolyphosphate = Phosphate + Triphosphate', 'H2O + UDP = UMP + Phosphate', "H2O + Uridine 3'-phosphate = Uridine + Phosphate", 'H2O + XMP = Phosphate + Xanthosine', 'H2O + beta-D-Fructose 2,6-bisphosphate = Phosphate + D-Fructose 6-phosphate', 'H2O + beta-Naphthyl phosphate = Phosphate + beta-Naphthol', 'H2O + dGDP = Phosphate + dGMP', 'H2O + dIMP = Phosphate + Deoxyinosine', 'H2O + dTMP = Phosphate + Thymidine', 'H2O + dTTP = Phosphate + TDP', 'H2O + myo-Inositol 4-phosphate = Phosphate + myo-Inositol', 'H2O + myo-Inositol hexakisphosphate = Phosphate + D-myo-Inositol 1,2,4,5,6-pentakisphosphate', 'H2O + myo-Inositol phosphate = myo-Inositol + Phosphate', 'H2O + o-Carboxyphenyl phosphate = Phosphate + Salicylate', 'H2O + p-Nitrophenylthymidine phosphate = p-Nitrophenol + dTMP', 'Inositol 2,4,5,6-tetrakisphosphate + H2O = Phosphate + Inositol 2,4,6-trisphosphate', 'Inositol 2,4,5-trisphosphate + H2O = Phosphate + Inositol 2,4-bisphosphate', 'Inositol 4,5,6-trisphosphate + H2O = Inositol 4,6-bisphosphate + Phosphate', 'Mannitol 1-phosphate + H2O = Mannitol + Phosphate', 'Mannose 6-phosphate + H2O = D-Mannose + Phosphate', 'N-(5-Phospho-4-pyridoxyl)benzylamine + H2O = 4-Pyridoxylbenzylamine + Phosphate', 'N-(5-Phospho-4-pyridoxyl)ethanolamine + H2O = 4-Pyridoxylethanolamine + Phosphate', 'N-(5-Phospho-4-pyridoxyl)phenylalanine + H2O = 4-Pyridoxylphenylalanine + Phosphate', 'N-Acetyl-D-mannosamine 6-phosphate + H2O = N-Acetylmannosamine + Phosphate', 'N-Acetylglucosamine 6-phosphate + H2O = N-Acetylglucosamine + Phosphate', 'Phosphate + Pyridoxamine = H2O + Pyridoxamine phosphate', 'Phosphate + Pyridoxine = H2O + Pyridoxine phosphate', 'Phosphorylase a + H2O = Phosphorylase b + Phosphate', 'Proteine tyrosine phosphate + H2O = Phosphate + Protein tyrosine', 'Riboflavin-5-phosphate + H2O = Phosphate + Riboflavin', 'Ribose 5-phosphate + H2O = Phosphate + beta-D-Ribopyranose', 'Sedoheptulose 1,7-bisphosphate + H2O = Phosphate + Sedoheptulose 7-phosphate', 'Thymolphthalein monophosphate + H2O = Thymolphthalein + Phosphate', 'UMP + H2O = Uridine + Phosphate', 'UTP + H2O = Phosphate + UDP', 'alpha-Naphthyl phosphate + H2O = alpha-Naphthol + Phosphate', 'beta-D-Thiogalactopyranoside 6-phosphate + H2O = beta-D-Thiogalactopyranoside + Phosphate', 'dAMP + H2O = Phosphate + Deoxyadenosine', "dCMP + H2O = Phosphate + 2'-Deoxycytidine", 'dGMP + H2O = Phosphate + Deoxyguanosine', "dUMP + H2O = Phosphate + 2'-Deoxyuridine", 'sn-Glycerol 1-phosphate + H2O = Phosphate + Glycerol', 'sn-Glycerol 3-phosphate + H2O = Phosphate + Glycerol'])

        # Each kinetic info class contains Entry objects from the sabio_rk module.
        # The following tests will tests ensure the KineticInfo fields contains the right Entry objects

        closest_entry_ids = km_data.closest_entry_ids
        self.assertTrue(len(closest_entry_ids) > 9 and len(closest_entry_ids) < 15)
        closest_entries = km_data.closest_entries  # this is a list of the most closely related expiremental entries
        # when I made this test, there were 13 entries. If there are more than 25 entries, something probably went wrong
        self.assertTrue(len(closest_entries) > 9 and len(closest_entries) < 15)

        # from the set of closest_entries, the entry with median km value is the median_entry
        # the median_entry is an Entry object defined in sabio_rk
        median_entry = km_data.median_entry
        self.assertEqual(median_entry.ec_number, '3.1.3.5')
        self.assertEqual(median_entry.entry_id, '30219')  # this is the specific Entry Number sabio assigns to each expiremental entry
        self.assertEqual(median_entry.km, '2.6E-6')  # this is the km value for this expiremental entry
        self.assertEqual(median_entry.vmax, '')  # this is the vmax for this expiremental entry
        self.assertEqual(median_entry.num_participants, [2, 2])  # this is a list to record the number of substrates and products
        self.assertEqual(median_entry.species, 'Mycoplasma fermentans')  # this is the species the expirement was done in
        self.assertEqual(median_entry.proximity, 1)  # this is closely related the expiremental species is to the modeler's species
        self.assertEqual(median_entry.reaction_id, '295')  # this is the Sabio assigned ID for this reaction in general

        # check the entry_id of the min and the max within km_data
        self.assertEqual(km_data.min_entry.entry_id, '30226')
        self.assertEqual(km_data.max_entry.entry_id, '30224')

        # now the vmax info fields will  be tested
        vmaxData = rxn.vmax_data
        self.assertEqual(vmaxData.ec_numbers, ['3.1.3.1', '3.1.3.11', '3.1.3.13', '3.1.3.16', '3.1.3.17', '3.1.3.18', '3.1.3.2', '3.1.3.22', '3.1.3.23', '3.1.3.24', '3.1.3.25', '3.1.3.29', '3.1.3.3', '3.1.3.31', '3.1.3.33', '3.1.3.34', '3.1.3.36', '3.1.3.37', '3.1.3.4', '3.1.3.43',
                                               '3.1.3.46', '3.1.3.48', '3.1.3.5', '3.1.3.56', '3.1.3.57', '3.1.3.6', '3.1.3.64', '3.1.3.7', '3.1.3.74', '3.1.3.76', '3.1.3.8', '3.1.3.9', '3.1.3.91', '3.1.3.93', '3.1.3.95'])  # it has two EC numbers for the same reaction. EC numbers are not a perfect system
        self.assertEqual(vmaxData.sabio_reaction_ids, ['10261', '10263', '10269', '10270', '10271', '10330', '10377', '10378', '10388', '10394', '10399', '10400', '10401', '10631', '10700', '10707', '10711', '10726', '10734', '1117', '1118', '11265', '11266', '11267', '11288', '11318', '1150', '1182', '12091', '12166', '124', '12404', '12405', '129', '13054', '13055', '1325', '13296', '13453', '13635', '13655', '13970', '14184', '14198', '14219', '1422', '1423', '1424', '1683', '1686', '1851', '186', '1891', '1923', '1965', '200', '207', '2084', '209', '210', '211', '2219', '2424', '246', '2493', '2503', '27', '2717', '2730', '282', '2874', '295', '304', '305', '309', '3140', '3170',
                                                       '3177', '3197', '3198', '3224', '3227', '339', '4095', '433', '444', '479', '480', '487', '501', '508', '581', '6096', '6098', '6445', '6446', '664', '6734', '6769', '697', '7185', '7186', '7187', '7188', '7189', '7190', '7191', '7192', '7193', '7194', '7195', '7196', '7278', '728', '736', '742', '75', '76', '7712', '7713', '7952', '7953', '7954', '796', '797', '7970', '7971', '7972', '7973', '7974', '7975', '7976', '7977', '8072', '8082', '816', '8203', '8204', '8205', '8206', '8749', '8902', '8953', '9216', '937', '9556', '9571', '9583', '9586', '9599', '9603', '964', '983', '9864', '9910', '9945', '9953', '9954', '9955', '9956', '9957', '9958', '9959', '9960', '9963'])
        self.assertEqual(vmaxData.lift_info, 'Lifted From 3.1.3')
        self.assertEqual(vmaxData.reaction_list, ["1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3',4',5'-trisphosphate) + H2O = 1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3',4'-bisphosphate) + Phosphate", "1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3',5'-bisphosphate) + H2O = 1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 3'-phosphate) + Phosphate", "1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 4',5'-bisphosphate) + H2O = 1,2-Dibutanoyl-sn-glycero-3-phospho-(1'-inositol 4'-phosphate) + Phosphate", '1-Phosphatidyl-D-myo-inositol 4,5-bisphosphate + H2O = Phosphate + 1-Phosphatidyl-1D-myo-inositol 4-phosphate', '1D-myo-Inositol 1,4,5,6-tetrakisphosphate + H2O = Inositol 1,4,6-trisphosphate + Phosphate', '1D-myo-Inositol 2-phosphate + H2O = myo-Inositol + Phosphate', "2'-Deoxyadenosine 3'-phosphate + H2O = Deoxyadenosine + Phosphate", "2'-Deoxycytidine 3'-phosphate + H2O = 2'-Deoxycytidine + Phosphate", "2'-Deoxyguanosine 3'-phosphate + H2O = Deoxyguanosine + Phosphate", "2'-Deoxyuridine 3'-phosphate + H2O = 2'-Deoxyuridine + Phosphate", '2,3-Diphosphoglycerate + H2O = Unknown + Phosphate', '2-Deoxy-D-glucose 6-phosphate + H2O = 2-Deoxy-D-glucose + Phosphate', '2-Deoxyglucose 6-phosphate + H2O = 2-Deoxyglucose + Phosphate', "3'-Phosphoadenylyl sulfate + H2O = Phosphate + Adenylylsulfate", '4-Pyridoxic acid 5-phosphate + H2O = 4-Pyridoxate + Phosphate', "5'-Ribonucleotide + H2O = Nucleoside + Phosphate", 'ADP + H2O = Diphosphate + Adenosine', 'ADP + H2O = Phosphate + AMP', 'ATP + H2O = Adenosine + Triphosphate', 'ATP + H2O = Phosphate + ADP', "Adenosine 3',5'-bisphosphate + H2O = Phosphate + AMP", 'Arabinose 5-phosphate + H2O = Arabinose + Phosphate', 'CDP + H2O = CMP + Phosphate', "Cytidine 3'-phosphate + H2O = Cytidine + Phosphate", 'D-Fructose 1,6-bisphosphate + H2O = D-Fructose monophosphate + Phosphate', 'D-Glucose 1-phosphate + H2O = D-Glucose + Phosphate', 'D-Glucose 6-phosphate + H2O = D-Glucose + Phosphate', 'D-Mannitol 1-phosphate + H2O = Phosphate + D-Mannitol', 'D-myo-Inositol 1,2,4,5,6-pentakisphosphate + H2O = Inositol 1,2,4,6-tetrakisphosphate + Phosphate', 'Diphosphate + H2O = Phosphate', 'Fructose 1,6-bisphosphate + H2O = Unknown + Phosphate', 'Fructose 1-phosphate + H2O = Fructose + Phosphate', 'Fructose 6-phosphate + H2O = Fructose + Phosphate', 'GDP + H2O = Phosphate + GMP', 'Galactose 1-phosphate + H2O = Galactose + Phosphate', 'Glucose + 3-Phospho-D-glyceroyl phosphate = Glucose 6-phosphate + 3-Phospho-D-glycerate', 'Glucose 1-phosphate + H2O = Glucose + Phosphate', 'Glycerate 3-phosphate + H2O = Glycerate + Phosphate', 'Glycerate 3-phosphate + H2O = Phosphate + D-Glycerate', 'Glycerol 3-phosphate + H2O = Glycerol + Phosphate', 'H2O + 1-Phosphatidyl-1D-myo-inositol 3,5-bisphosphate = Phosphate + 1-Phosphatidyl-1D-myo-inositol 5-phosphate', 'H2O + 1-Phosphatidyl-1D-myo-inositol 3-phosphate = Phosphate + 1-Phosphatidyl-D-myo-inositol', 'H2O + 1-Phospho-D-glycerate = D-Glycerate + Phosphate', 'H2O + 1D-myo-Inositol 1,3,4,5-tetrakisphosphate = Phosphate + 1D-myo-Inositol 1,3,4-trisphosphate', 'H2O + 1D-myo-Inositol 1,3,4-trisphosphate = Phosphate + D-myo-Inositol 3,4-bisphosphate', 'H2O + 1D-myo-Inositol 1,4-bisphosphate = Phosphate + myo-Inositol 4-phosphate', 'H2O + 1D-myo-Inositol 1-phosphate = Phosphate + myo-Inositol', 'H2O + 1D-myo-Inositol 3-phosphate = Phosphate + myo-Inositol', 'H2O + 2-Chloro-4-nitrophenyl phosphate = Phosphate + 2-Chloro-4-nitrophenol', 'H2O + 2-Phospho-D-glycerate = Phosphate + D-Glycerate', 'H2O + 2-Phosphoglycolate = Phosphate + Glycolate', 'H2O + 3-O-Methylfluorescein phosphate = Phosphate + 3-O-Methylfluorescein', 'H2O + 4-Chlorophenyl phosphate = Phosphate + 4-Chlorophenol', 'H2O + 4-Cyanophenyl phosphate = Phosphate + 4-Cyanophenol', 'H2O + 4-Methylumbelliferyl phosphate = Phosphate + 4-Methylumbelliferone', 'H2O + 4-Nitrophenyl phenyl phosphonate = p-Nitrophenol + Phenyl phosphonate', 'H2O + 4-Nitrophenyl phosphate = Phosphate + p-Nitrophenol', 'H2O + 4-Trifluoromethylphenyl phosphate = Phosphate + 4-Trifluoromethylphenol', "H2O + 5'-Phosphopolynucleotide = Phosphate + Polynucleotide", 'H2O + 5-Fluoro-4-methylumbelliferyl phosphate = Phosphate + 5-Fluoro-4-methylumbelliferone', 'H2O + 6,8-Difluoro-4-methylumbelliferyl phosphate = Phosphate + 6,8-Difluoro-4-methylumbelliferone', 'H2O + 6-Fluoro-4-methylumbelliferyl phosphate = Phosphate + 6-Fluoro-4-methylumbelliferone', 'H2O + 6-Phosphogluconate = Phosphate + Gluconate', "H2O + 7-Methylguanosine 5'-phosphate = Phosphate + 7-Methylguanosine", 'H2O + 8-Fluoro-4-methylumbelliferyl phosphate = Phosphate + 8-Fluoro-4-methylumbelliferone', 'H2O + AMP = Phosphate + Adenosine', "H2O + Adenosine 2'-phosphate = Adenosine + Phosphate", "H2O + Adenosine 3'-phosphate = Phosphate + Adenosine", 'H2O + Bis-4-nitrophenyl phosphate = p-Nitrophenol + 4-Nitrophenyl phosphate', 'H2O + CMP = Cytidine + Phosphate', 'H2O + CTP = CDP + Phosphate', 'H2O + Casein kinase I epsilon phosphorylated = Phosphate + Casein kinase I epsilon', 'H2O + Ceramide 1-phosphate = Phosphate + N-Acylsphingosine', 'H2O + Choline phosphate = Phosphate + Choline', 'H2O + D-Fructose 1,6-bisphosphate = Phosphate + D-Fructose 6-phosphate', 'H2O + D-Fructose 1,6-bisphosphate = Phosphate + beta-D-Fructose 6-phosphate',
                                                  'H2O + D-Fructose 2,6-bisphosphate = Phosphate + D-Fructose 6-phosphate', 'H2O + D-Fructose 2,6-bisphosphate = Phosphate + beta-D-Fructose 6-phosphate', 'H2O + D-Galactose 1-phosphate = Phosphate + D-Galactose', 'H2O + D-Mannose 6-phosphate = D-Mannose + Phosphate', 'H2O + D-O-Phosphoserine = Phosphate + D-Serine', 'H2O + D-myo-Inositol 1,3-bisphosphate = Phosphate + 1D-myo-Inositol 1-phosphate', 'H2O + D-myo-Inositol 1,4,5-trisphosphate = Phosphate + 1D-myo-Inositol 1,4-bisphosphate', "H2O + Deoxythymidine 3',5'-diphosphate = Phosphate + dTMP", "H2O + Deoxythymidine 3'-phosphate = Phosphate + Thymidine", 'H2O + Ethanolamine phosphate = Phosphate + Ethanolamine', 'H2O + Fructose 1,6-bisphosphate = Phosphate + Fructose 6-phosphate', 'H2O + GMP = Phosphate + Guanosine', 'H2O + GTP = GDP + Phosphate', 'H2O + Glucosamine 6-phosphate = Phosphate + Glucosamine', 'H2O + Glucose 6-phosphate = Glucose + Phosphate', 'H2O + Glycerate 2,3-bisphosphate = Glycerate 3-phosphate + Phosphate', 'H2O + Glycerate 2,3-bisphosphate = Phosphate + 3-Phospho-D-glycerate', 'H2O + Glycerol 2-phosphate = Phosphate + Glycerol', "H2O + Guanosine 3'-phosphate = Phosphate + Guanosine", 'H2O + IMP = Phosphate + Inosine', 'H2O + Inositol 1,2,3,4,5-pentakisphosphate = Phosphate + Inositol 1,2,3,4-tetrakisphosphate', 'H2O + Inositol 1,2,4,5-tetrakisphosphate = Inositol 1,2,4-trisphosphate + Phosphate', 'H2O + Inositol 4,5-bisphosphate = Phosphate + myo-Inositol 4-phosphate', 'H2O + L-Galactose 1-phosphate = Phosphate + L-Galactose', 'H2O + L-Phosphotyrosine = Phosphate + L-Tyrosine', 'H2O + L-Threonine O-3-phosphate = L-Threonine + Phosphate', 'H2O + N-(5-Phospho-4-pyridoxyl)glycine = 4-Pyridoxylglycine + Phosphate', 'H2O + N-Acetylneuraminate 9-phosphate = Phosphate + N-Acetylneuraminate', 'H2O + N-Acylneuraminate 9-phosphate = Phosphate + N-Acylneuraminate', 'H2O + NADP+ = Phosphate + NAD+', 'H2O + NADPH = NADH + Phosphate', 'H2O + O-Phospho-L-serine = L-Serine + Phosphate', 'H2O + O-Phospho-L-serine = Serine + Phosphate', 'H2O + O-Phospho-tau-protein = Phosphate + tau-Protein', 'H2O + Phenolic phosphate = Phosphate + Phenol', 'H2O + Phenolphthalein diphosphate = Phosphate + Phenolphthalein', 'H2O + Phosphoenolpyruvate = Phosphate + Pyruvate', 'H2O + Phosphotyrosine = Phosphate + Tyrosine', 'H2O + Propan-1-ol 2-phosphate = Phosphate + 1-Propanol', 'H2O + Pyridoxal phosphate = Phosphate + Pyridoxal', 'H2O + Ribulose 5-phosphate = Ribulose + Phosphate', 'H2O + Sorbitol 6-phosphate = Phosphate + Sorbitol', 'H2O + Sphinganine 1-phosphate = Phosphate + Sphinganine', 'H2O + SpoIIAA-phosphorylated = Phosphate + SpoIIAA', 'H2O + Sucrose 6-phosphate = Sucrose + Phosphate', 'H2O + TDP = Phosphate + Thiamine monophosphate', 'H2O + Tetrapolyphosphate = Phosphate + Triphosphate', 'H2O + UDP = UMP + Phosphate', "H2O + Uridine 3'-phosphate = Uridine + Phosphate", 'H2O + XMP = Phosphate + Xanthosine', 'H2O + beta-D-Fructose 2,6-bisphosphate = Phosphate + D-Fructose 6-phosphate', 'H2O + beta-Naphthyl phosphate = Phosphate + beta-Naphthol', 'H2O + dGDP = Phosphate + dGMP', 'H2O + dIMP = Phosphate + Deoxyinosine', 'H2O + dTMP = Phosphate + Thymidine', 'H2O + dTTP = Phosphate + TDP', 'H2O + myo-Inositol 4-phosphate = Phosphate + myo-Inositol', 'H2O + myo-Inositol hexakisphosphate = Phosphate + D-myo-Inositol 1,2,4,5,6-pentakisphosphate', 'H2O + myo-Inositol phosphate = myo-Inositol + Phosphate', 'H2O + o-Carboxyphenyl phosphate = Phosphate + Salicylate', 'H2O + p-Nitrophenylthymidine phosphate = p-Nitrophenol + dTMP', 'Inositol 2,4,5,6-tetrakisphosphate + H2O = Phosphate + Inositol 2,4,6-trisphosphate', 'Inositol 2,4,5-trisphosphate + H2O = Phosphate + Inositol 2,4-bisphosphate', 'Inositol 4,5,6-trisphosphate + H2O = Inositol 4,6-bisphosphate + Phosphate', 'Mannitol 1-phosphate + H2O = Mannitol + Phosphate', 'Mannose 6-phosphate + H2O = D-Mannose + Phosphate', 'N-(5-Phospho-4-pyridoxyl)benzylamine + H2O = 4-Pyridoxylbenzylamine + Phosphate', 'N-(5-Phospho-4-pyridoxyl)ethanolamine + H2O = 4-Pyridoxylethanolamine + Phosphate', 'N-(5-Phospho-4-pyridoxyl)phenylalanine + H2O = 4-Pyridoxylphenylalanine + Phosphate', 'N-Acetyl-D-mannosamine 6-phosphate + H2O = N-Acetylmannosamine + Phosphate', 'N-Acetylglucosamine 6-phosphate + H2O = N-Acetylglucosamine + Phosphate', 'Phosphate + Pyridoxamine = H2O + Pyridoxamine phosphate', 'Phosphate + Pyridoxine = H2O + Pyridoxine phosphate', 'Phosphorylase a + H2O = Phosphorylase b + Phosphate', 'Proteine tyrosine phosphate + H2O = Phosphate + Protein tyrosine', 'Riboflavin-5-phosphate + H2O = Phosphate + Riboflavin', 'Ribose 5-phosphate + H2O = Phosphate + beta-D-Ribopyranose', 'Sedoheptulose 1,7-bisphosphate + H2O = Phosphate + Sedoheptulose 7-phosphate', 'Thymolphthalein monophosphate + H2O = Thymolphthalein + Phosphate', 'UMP + H2O = Uridine + Phosphate', 'UTP + H2O = Phosphate + UDP', 'alpha-Naphthyl phosphate + H2O = alpha-Naphthol + Phosphate', 'beta-D-Thiogalactopyranoside 6-phosphate + H2O = beta-D-Thiogalactopyranoside + Phosphate', 'dAMP + H2O = Phosphate + Deoxyadenosine', "dCMP + H2O = Phosphate + 2'-Deoxycytidine", 'dGMP + H2O = Phosphate + Deoxyguanosine', "dUMP + H2O = Phosphate + 2'-Deoxyuridine", 'sn-Glycerol 1-phosphate + H2O = Phosphate + Glycerol', 'sn-Glycerol 3-phosphate + H2O = Phosphate + Glycerol'])
        # Each kinetic info class contains Entry objects from the sabio_interface module.

        # The following tests will tests ensure the KineticInfo fields contains the right Entry objects

        closest_entry_ids = vmaxData.closest_entry_ids
        self.assertTrue(len(closest_entry_ids) > 6 and len(closest_entry_ids) < 10)
        closest_entries = vmaxData.closest_entries  # this is a list of the most closely related expiremental entries
        # when I made this test, there were 45 entries. If there are more than 50 entries, something probably went wrong
        self.assertTrue(len(closest_entries) > 6 and len(closest_entries) < 10)

        # from the set of closest_entries, the entry with median km value is the median_entry
        # the median_entry is an Entry object defined in sabio_rk
        median_entry = vmaxData.median_entry
        self.assertEqual(median_entry.ec_number, '3.1.3.22')
        self.assertEqual(median_entry.entry_id, '22920')  # this is the specific Entry Number sabio assigns to each expiremental entry
        self.assertEqual(median_entry.km, '9.0E-6')  # this is the km value for this expiremental entry
        self.assertEqual(median_entry.vmax, '4.33333333E-5')  # this is the vmax for this expiremental entry
        self.assertEqual(median_entry.num_participants, [2, 2])  # this is a list to record the number of substrates and products
        self.assertEqual(median_entry.species, 'Streptococcus bovis')  # this is the species the expirement was done in
        self.assertEqual(median_entry.proximity, 6)  # this is closely related the expiremental species is to the modeler's species
        self.assertEqual(median_entry.reaction_id, '581')  # this is the Sabio assigned ID for this reaction in general

        # check the entry_id of the min and the max within km_data
        self.assertEqual(vmaxData.min_entry.entry_id, '22921')
        self.assertEqual(vmaxData.max_entry.entry_id, '16039')

        ##################################################
        ##################################################
        # this Entry is an example of where less than 3 relevant entries are found
        # If two relevant entries are found the min and max are filled in but the median is blank
        # if one relevant entry is found, the median is filled in, but the in and max are blank

        rxn = next((rxn for rxn in rxns if rxn.id == 'fmn_reductase'), None)
        self.assertNotEqual(rxn, None)

        self.assertEqual(rxn.reaction_ids, ['5301'])
        # both km and vmax should have KineticInfo objects in their fields
        self.assertNotEqual(rxn.km_data, None)
        self.assertNotEqual(rxn.vmax_data, None)

        # Now the KineticInfo fields will be tested
        # First is some summaries of the expiremental data.
        km_data = rxn.km_data

        self.assertEqual(km_data.lift_info, 'Lift Not Used')

        # Even though all the reactions are the same, the reaction has many EC numbers (because EC numbers are an imperfect system)
        self.assertEqual(km_data.ec_numbers, ['1.14.13', '1.14.13.7', '1.5.1', '1.5.1.30', '1.5.1.39'])
        # only one reaction ID should be present
        self.assertEqual(km_data.sabio_reaction_ids, ['5301'])
        self.assertEqual(km_data.reaction_list, ['NAD+ + Reduced FMN = NADH + H+ + Riboflavin-5-phosphate'])
        # Each kinetic info class contains Entry objects from the sabio_interface module.

        # The following tests will tests ensure the KineticInfo fields contains the right Entry objects

        closest_entries = km_data.closest_entries
        self.assertEqual(len(closest_entries), 2)
        # make sure that median is blank but min and max are filled in
        self.assertEqual(km_data.median_entry, None)
        self.assertFalse(km_data.min_entry == None or km_data.max_entry == None)

        vmaxData = rxn.vmax_data

        self.assertEqual(vmaxData.lift_info, 'Lift Not Used')

        # Even though all the reactions are the same, the reaction has many EC numbers (because EC numbers are an imperfect system)
        self.assertEqual(vmaxData.ec_numbers, ['1.14.13', '1.14.13.7', '1.5.1', '1.5.1.30', '1.5.1.39'])
        # only one reaction ID should be present
        self.assertEqual(vmaxData.sabio_reaction_ids, ['5301'])
        self.assertEqual(vmaxData.reaction_list, ['NAD+ + Reduced FMN = NADH + H+ + Riboflavin-5-phosphate'])
        # Each kinetic info class contains Entry objects from the sabio_interface module.

        # The following tests will tests ensure the KineticInfo fields contains the right Entry objects

        closest_entries = vmaxData.closest_entries
        self.assertEqual(len(closest_entries), 1)
        # make sure that median is blank but min and max are filled in
        self.assertNotEqual(vmaxData.median_entry, None)
        self.assertTrue(vmaxData.min_entry == None or vmaxData.max_entry == None)

        ##################################################
        ##################################################
        # todo the next entry is a case where Sabio did find the queried reaction, however it did not find any vmax infromation
        # it only found km
