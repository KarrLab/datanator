# -*- coding: utf-8 -*-

""" Test of common schema

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-07-31
:Copyright: 2017, Karr Lab
:License: MIT
"""
import unittest
from sqlalchemy.orm import sessionmaker
from kinetic_datanator.core import common_schema
import tempfile
import shutil


class ShortTestCommonSchema(unittest.TestCase):
    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        self._is_local = False

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_working(self):
        if self._is_local:
            cs = common_schema.CommonSchema(name = 'aggregate', clear_content = True,
                                            load_content=False, download_backup=False,
                                            max_entries = 3)
            cs._is_local = True
        else:
            cs = common_schema.CommonSchema(cache_dirname = self.cache_dirname , clear_content = True,
                load_content= False, download_backup=False, max_entries = 3)
            cs._is_local = False

        cs.load_content()
        session = cs.session

        dataset = session.query(common_schema.AbundanceDataSet).filter_by(file_name = '882/882-Desulfo_Form_Exp_SC_zhang_2006.txt').first()
        self.assertEqual(dataset.score , 2.47)

        subunit = session.query(common_schema.ProteinSubunit).filter_by(uniprot_id = 'Q72DH8').first()
        abundance = session.query(common_schema.AbundanceData).filter_by(subunit_id = subunit.subunit_id).first()
        self.assertEqual(abundance.abundance, 78.0)

        metadata = session.query(common_schema.Metadata).filter_by(name = '882/882-Desulfo_Lac_Exp_SC_zhang_2006.txt').first()
        taxon = session.query(common_schema._metadata_taxon).filter_by(_metadata_id = metadata.id).first()
        self.assertEqual(taxon.taxon_id, 882)

        subunit = session.query(common_schema.ProteinSubunit).filter_by(subunit_name = 'Histone deacetylase 5').first()
        self.assertEqual(subunit.uniprot_id, 'Q9UQL6')
        self.assertEqual(subunit.entrez_id, 10014)
        complex_  = session.query(common_schema.ProteinComplex).get(subunit.proteincomplex_id)
        self.assertEqual(complex_.complex_name, 'BCL6-HDAC5 complex')

        metadata = session.query(common_schema.Metadata).filter_by(name = 'BCL6-HDAC7 complex').first()
        resource = session.query(common_schema._metadata_resource).filter_by(_metadata_id = metadata.id).first()
        pubmed = session.query(common_schema.Resource).get(resource.resource_id)
        self.assertEqual(pubmed._id, '11929873')


        subunit = session.query(common_schema.ProteinSubunit).filter_by(subunit_name = 'TFAP2A').first()
        self.assertEqual(subunit.uniprot_id, 'P05549')
        self.assertEqual(subunit.class_name, 'Zipper-Type')
        binding = session.query(common_schema.DNABindingDataset).filter_by(subunit_id = subunit.subunit_id).first()
        data = session.query(common_schema.DNABindingData).filter_by(dataset_id = binding.dataset_id).first()
        self.assertEqual(data.position, 1)
        self.assertEqual(data.frequency_c, 94)

        metadata = session.query(common_schema.Metadata).filter_by(name = 'RUNX1').first()
        resource = session.query(common_schema._metadata_resource).filter_by(_metadata_id = metadata.id).first()
        pubmed = session.query(common_schema.Resource).get(resource.resource_id)
        self.assertEqual(pubmed._id, '8413232')


        compound = session.query(common_schema.Compound).filter_by(name = 'Deoxyuridine').first()
        self.assertEqual(compound.description, "2'-Deoxyuridine is a naturally occurring nucleoside. It is similar in chemical structure to uridine, but without the 2'-hydroxyl group.  It is considered to be an antimetabolite that is converted to deoxyuridine triphosphate during DNA synthesis.")
        structure = session.query(common_schema.Structure).get(compound.structure_id)
        self.assertEqual(structure._value_inchi , 'InChI=1S/C9H12N2O5/c12-4-6-5(13)3-8(16-6)11-2-1-7(14)10-9(11)15/h1-2,5-6,8,12-13H,3-4H2,(H,10,14,15)/t5-,6+,8+/m0/s1')

        metadata  = session.query(common_schema.Metadata).filter_by(name = 'Deoxycytidine').first()
        synonym = session.query(common_schema._metadata_synonym).filter_by(_metadata_id = metadata.id).first()
        name = session.query(common_schema.Synonym).get(synonym.synonym_id)
        self.assertEqual(name.name , '4-amino-1-[(2R,4S,5R)-4-hydroxy-5-(hydroxymethyl)oxolan-2-yl]-1,2-dihydropyrimidin-2-one')

        compound = session.query(common_schema.Compound).filter_by(compound_name = 'Peptide').first()
        self.assertEqual(compound._is_name_ambiguous , 1)
        structure = session.query(common_schema.Structure).get(compound.structure_id)
        self.assertEqual(structure._value_smiles, '[*]C(N)C(=O)NC([*])C(O)=O')

        metadata = session.query(common_schema.Metadata).filter_by(name = 'Kinetic Law 1').first()
        cell_line = session.query(common_schema._metadata_cell_line).filter_by(_metadata_id = metadata.id).first()
        name = session.query(common_schema.CellLine).get(cell_line.cell_line_id)
        self.assertEqual(name.name , 'variant DSAI (N76D/N87S/S103A/V104I)')


# class LongTestCommonSchema(unittest.TestCase):
#     def setUp(self):
#         self.cache_dirname = tempfile.mkdtemp()
#
#     def tearDown(self):
#         shutil.rmtree(self.cache_dirname)
#
#     def test_working(self):
#         # cs = common_schema.CommonSchema(cache_dirname = self.cache_dirname , clear_content = True, load_content=False, download_backup=False)
#         cs = common_schema.CommonSchema(name = 'aggregate', clear_content = True,
#                                         load_content=False, download_backup=False
#                                         max_entries = 20)
#         cs.load_content()
#         session = cs.session
