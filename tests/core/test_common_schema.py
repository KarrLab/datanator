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



class LoadingTestCommonSchema(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.cs = common_schema.CommonSchema(cache_dirname = self.cache_dirname,
                                clear_content = True, load_entire_small_DBs = False,
                                download_backup= False, load_content = True, max_entries = 10,
                                verbose = True)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)


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
        subunit = session.query(common_schema.ProteinSubunit).filter_by(subunit_name = 'Histone deacetylase 5').first()
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
        subunit = session.query(common_schema.ProteinSubunit).filter_by(gene_name = 'TFAP2A').first()
        self.assertEqual(subunit.uniprot_id, 'P05549')
        self.assertEqual(subunit.class_name, 'Basic helix-span-helix factors (bHSH)')
        binding = session.query(common_schema.DNABindingDataset).filter_by(subunit_id = subunit.subunit_id).first()
        data = session.query(common_schema.DNABindingData).filter_by(dataset_id = binding.dataset_id).first()
        self.assertEqual(data.position, 1)
        self.assertEqual(data.frequency_g, 185)

        metadata = session.query(common_schema.Metadata).filter_by(name = 'RUNX1 Binding Motif').first()
        self.assertEqual(metadata.resource[0]._id, '8413232')

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


        q = session.query(common_schema.KineticLaw) \
            .join((common_schema.Metadata, common_schema.KineticLaw._metadata)).join((common_schema.Resource, common_schema.Metadata.resource))\
            .filter(common_schema.Resource.namespace == 'ec-code').filter(common_schema.Resource._id.in_(['3.4.21.62']))

        compare = session.query(common_schema.KineticLaw).filter_by(enzyme_id = q.first().enzyme_id).all()

        self.assertEqual(set([n.kineticlaw_id for n in q.all()]),
            set([c.kineticlaw_id for c in compare]))



        resource = session.query(common_schema.Resource).filter_by(namespace = 'ec-code').filter_by(_id = '3.4.21.62').all()
        self.assertEqual(len(resource), 1)
        print resource[0]._metadata

    def test_intact_added(self):
        session = self.cs.session

        interact = session.query(common_schema.ProteinInteractions).filter_by(participant_a = 'uniprotkb:P49418').all()
        self.assertEqual(set([c.site_b for c in interact]), set(['binding-associated region:1063-1070(MINT-376288)', '-']))

        interact = session.query(common_schema.ProteinInteractions).filter_by(participant_a = 'intact:EBI-7121765').first()
        self.assertEqual(interact._metadata.resource[0]._id, '10542231|mint')

    def test_uniprot_added(self):
        session = self.cs.session

        uni = session.query(common_schema.ProteinSubunit).filter_by(uniprot_id = 'Q72DQ8').first()
        self.assertEqual(uni.subunit_name, 'PYRH_DESVH')
        self.assertEqual(uni.entrez_id, 2795170)
        self.assertEqual(uni.length, 238)
