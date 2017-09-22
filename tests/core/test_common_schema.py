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
from kinetic_datanator.data_source import pax
import tempfile
import shutil
import random


class ShortTestCommonSchema(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.cs = common_schema.CommonSchema(cache_dirname = self.cache_dirname,
                                clear_content = True,
                                load_content= True, download_backup= False,
                                max_entries = 10, verbose = True)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)

    def test_fill_missing_subunit_info(self):
        session = self.cs.session
        subunit = session.query(common_schema.ProteinSubunit).filter_by(uniprot_id = 'Q01196').first()
        self.assertEqual(subunit.length, 453)
        self.assertEqual(subunit.mass, 48737)

    def test_fill_missing_ncbi_names(self):
        session = self.cs.session
        taxon = session.query(common_schema.Taxon).filter_by(ncbi_id = 882).first()
        self.assertEqual(taxon.name, 'Desulfovibrio vulgaris str. Hildenborough')

    def test_pax_added(self):
        session = self.cs.session
        pax_compare = pax.Pax(cache_dirname = self.cache_dirname, download_backup = True, load_content = False)
        pax_session = pax_compare.session

        dataset = session.query(common_schema.AbundanceDataSet).first()
        comparison = pax_session.query(pax.Dataset).filter_by(file_name = dataset.file_name).first()
        self.assertEqual(dataset.score, comparison.score)
        self.assertEqual(dataset.weight, comparison.weight)
        self.assertEqual(dataset.coverage, comparison.coverage)

        metadata = session.query(common_schema.Metadata).filter_by(name = dataset.file_name).first()
        taxon = session.query(common_schema._metadata_taxon).filter_by(_metadata_id = metadata.id).first()
        self.assertEqual(taxon.taxon_id, comparison.taxon_ncbi_id)

    def test_corum_added(self):
        session = self.cs.session
        subunit = session.query(common_schema.ProteinSubunit).filter_by(subunit_name = 'Histone deacetylase 5 (HD5) (EC 3.5.1.98) (Antigen NY-CO-9)').first()
        self.assertEqual(subunit.uniprot_id, 'Q9UQL6')
        self.assertEqual(subunit.entrez_id, 10014)
        complex_  = session.query(common_schema.ProteinComplex).get(subunit.proteincomplex_id)
        self.assertEqual(complex_.complex_name, 'BCL6-HDAC5 complex')

        metadata = session.query(common_schema.Metadata).filter_by(name = 'BCL6-HDAC7 complex').first()
        resource = session.query(common_schema._metadata_resource).filter_by(_metadata_id = metadata.id).first()
        pubmed = session.query(common_schema.Resource).get(resource.resource_id)
        self.assertEqual(pubmed._id, '11929873')

    def test_jaspar_added(self):
        session = self.cs.session
        subunit = session.query(common_schema.ProteinSubunit).filter_by(subunit_name = 'Transcription factor AP-2-alpha (AP2-alpha) (AP-2 transcription factor) (Activating enhancer-binding protein 2-alpha) (Activator protein 2) (AP-2)').first()
        self.assertEqual(subunit.uniprot_id, 'P05549')
        self.assertEqual(subunit.class_name, 'Zipper-Type')
        binding = session.query(common_schema.DNABindingDataset).filter_by(subunit_id = subunit.subunit_id).first()
        data = session.query(common_schema.DNABindingData).filter_by(dataset_id = binding.dataset_id).first()
        self.assertEqual(data.position, 1)
        self.assertEqual(data.frequency_g, 185)

        metadata = session.query(common_schema.Metadata).filter_by(name = 'RUNX1').first()
        resource = session.query(common_schema._metadata_resource).filter_by(_metadata_id = metadata.id).first()
        pubmed = session.query(common_schema.Resource).get(resource.resource_id)
        self.assertEqual(pubmed._id, '8413232')

    def test_ecmdb_added(self):
        session = self.cs.session
        compound = session.query(common_schema.Compound).filter_by(name = 'Deoxyuridine').first()
        self.assertEqual(compound.description, "2'-Deoxyuridine is a naturally occurring nucleoside. It is similar in chemical structure to uridine, but without the 2'-hydroxyl group.  It is considered to be an antimetabolite that is converted to deoxyuridine triphosphate during DNA synthesis.")
        structure = session.query(common_schema.Structure).get(compound.structure_id)
        self.assertEqual(structure._value_inchi , 'InChI=1S/C9H12N2O5/c12-4-6-5(13)3-8(16-6)11-2-1-7(14)10-9(11)15/h1-2,5-6,8,12-13H,3-4H2,(H,10,14,15)/t5-,6+,8+/m0/s1')

        metadata  = session.query(common_schema.Metadata).filter_by(name = 'Deoxycytidine').first()
        synonym = session.query(common_schema._metadata_synonym).filter_by(_metadata_id = metadata.id).first()
        name = session.query(common_schema.Synonym).get(synonym.synonym_id)
        self.assertEqual(name.name , '4-amino-1-[(2R,4S,5R)-4-hydroxy-5-(hydroxymethyl)oxolan-2-yl]-1,2-dihydropyrimidin-2-one')

    def test_sabio_added(self):
        session = self.cs.session
        compound = session.query(common_schema.Compound).filter_by(compound_name = 'Peptide').first()
        self.assertEqual(compound._is_name_ambiguous , 1)
        structure = session.query(common_schema.Structure).get(compound.structure_id)
        self.assertEqual(structure._value_smiles, '[*]C(N)C(=O)NC([*])C(O)=O')

        metadata = session.query(common_schema.Metadata).filter_by(name = 'Kinetic Law 1').first()
        cell_line = session.query(common_schema._metadata_cell_line).filter_by(_metadata_id = metadata.id).first()
        name = session.query(common_schema.CellLine).get(cell_line.cell_line_id)
        self.assertEqual(name.name , 'variant DSAI (N76D/N87S/S103A/V104I)')


# class LongTestCommonSchema(unittest.TestCase):
#     @classmethod
#     def setUpClass(self):
#         self.cache_dirname = tempfile.mkdtemp()
#         self.cs = common_schema.CommonSchema(cache_dirname = self.cache_dirname,
#                                 clear_content = True,
#                                 load_content= True, download_backup= False,
#                                 max_entries = 20, verbose = True)
#
#     @classmethod
#     def tearDownClass(self):
#         shutil.rmtree(self.cache_dirname)
#
#
#     def test_pax_added(self):
#         session = self.cs.session
#         pax_compare = pax.Pax(cache_dirname = self.cache_dirname, download_backup = True, load_content = False, clear_content = False)
#         pax_session = pax_compare.session
#
#         dataset_id = random.choice([1,2,3,4])
#         pax_dataset = pax_session.query(pax.Dataset).get(dataset_id)
#
#         abundance_dataset = session.query(common_schema.AbundanceDataSet).filter_by(file_name = pax_dataset.file_name).first()
#         if abundance_dataset:
#             self.assertEqual(abundance_dataset.score, pax_dataset.score)
#             self.assertEqual(abundance_dataset.weight, pax_dataset.weight)
#             self.assertEqual(abundance_dataset.coverage, pax_dataset.coverage)
#
#         metadata = session.query(common_schema.Metadata).filter_by(name = abundance_dataset.file_name).first()
#         taxon = session.query(common_schema._metadata_taxon).filter_by(_metadata_id = metadata.id).first()
#         self.assertEqual(taxon.taxon_id, pax_dataset.taxon_ncbi_id)
#
#     def test_corum_added(self):
#         session = self.cs.session
#         complex_ = session.query(common_schema.ProteinComplex).filter_by(complex_name = 'Condensin I complex').first()
#         self.assertEqual(complex_.funcat_id, '10.03.01.01.11;10.03.04.03;10.03.04.05;42.10.03;70.10.03')
#
#         subunit = session.query(common_schema.ProteinSubunit).filter_by(subunit_name = 'Hermansky-Pudlak syndrome 1 protein').first()
#         self.assertEqual(subunit.entrez_id, 3257)
#         cmplx = session.query(common_schema.ProteinComplex).get(subunit.proteincomplex_id)
#         self.assertEqual(cmplx.complex_name, 'BLOC-3 (biogenesis of lysosome-related organelles complex 3)')
#
#         metadata_subq = session.query(common_schema.Metadata).filter_by(name = 'Condensin I complex').subquery()
#         resource = session.query(common_schema.Resource).join((metadata_subq, common_schema.Resource._metadata)).first()
#         self.assertEqual(resource._id, '11136719')
#
#     def test_jaspar_added(self):
#         session = self.cs.session
#         subunit = session.query(common_schema.ProteinSubunit).filter_by(uniprot_id = 'Q01295').first()
#         self.assertEqual(subunit.entrez_id, 44505)
#         self.assertEqual(subunit.gene_name , 'br')
#
#         binding = session.query(common_schema.DNABindingDataset).filter_by(subunit_id = subunit.subunit_id).first()
#         data = session.query(common_schema.DNABindingData).filter_by(dataset_id = binding.dataset_id).filter_by(position = 1).first()
#         self.assertEqual(data.frequency_g, 4)
#
#         metadata_subq = session.query(common_schema.Metadata).filter_by(name = 'PAX5').subquery()
#         taxon = session.query(common_schema.Taxon).join((metadata_subq, common_schema.Taxon._metadata)).first()
#         self.assertEqual(taxon.name, 'Mus musculus')
#
#     def test_ecmdb_added(self):
#         session = self.cs.session
#         compound = session.query(common_schema.Compound).filter_by(compound_name = 'Acetic acid').first()
#         structure = session.query(common_schema.Structure).get(compound.structure_id)
#         self.assertEqual(structure._value_inchi, 'InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)')
#
#
#         concentration = session.query(common_schema.Concentration).filter_by(compound_id = compound.compound_id)
#         self.assertEqual(concentration.count(),1)
#         self.assertEqual(concentration.first().value, 658.0)
#
#         metadata_subq = session.query(common_schema.Metadata).filter_by(name = 'Adenosine monophosphate').subquery()
#         cell_line = session.query(common_schema.CellLine).join((metadata_subq, common_schema.CellLine._metadata)).all()
#         self.assertEqual(len(cell_line), 3)
#
#     def test_sabio_added(self):
#         session = self.cs.session
#         compound = session.query(common_schema.Compound).filter_by(compound_name = 'Reduced FMN')
#         self.assertEqual(compound.count(),1)
