# -*- coding: utf-8 -*-

""" Test of common schema

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-07-31
:Copyright: 2017, Karr Lab
:License: MIT
"""
import unittest
from kinetic_datanator.core import flask_common_schema, models
from kinetic_datanator.data_source import pax
import flask_whooshalchemy
import flask
import tempfile
import shutil
import random
import os
from six.moves import reload_module

class DownloadTestFlaskCommonSchema(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.flk = flask_common_schema.FlaskCommonSchema(cache_dirname=self.cache_dirname)

        for item in self.flk.text_indicies:
            flask_whooshalchemy.whoosh_index(self.flk.app, item)

    @classmethod
    def tearDownClass(self):
        models.db.session.remove()
        models.db.drop_all()
        shutil.rmtree(self.cache_dirname)

    def test_serialization(self):
        test_json = models.Compound.query.filter_by(compound_name='2-Hydroxyisocaproate').first().serialize()
        self.assertEqual(test_json['compound_name'], '2-Hydroxyisocaproate')
        self.assertEqual(test_json['_is_name_ambiguous'], False)

        test_json = models.AbundanceDataSet.query.filter_by(file_name='882/882-Desulfo_Form_Exp_SC_zhang_2006.txt').first().serialize()
        self.assertEqual(test_json['score'], 2.47)
        self.assertEqual(test_json['weight'], 100)

        test_json = models.Compound.query.filter_by(compound_name='Adenine').first().serialize(metadata=True, relationships=True)
        self.assertEqual(test_json['compound_name'], 'Adenine')
        self.assertEqual(test_json['relationships']['structure']['_value_inchi'], 'InChI=1S/C5H5N5/c6-4-3-5(9-1-7-3)10-2-8-4/h1-2H,(H3,6,7,8,9,10)')
        self.assertEqual(test_json['metadata']['cell_compartment'][0]['id'], 21)


    def test_data_loaded(self):
        session = self.flk.session
        taxon = session.query(models.Taxon).filter_by(ncbi_id = 882).first()
        self.assertEqual(taxon.name, 'Desulfovibrio vulgaris str. Hildenborough')

        subunit = session.query(models.ProteinSubunit).filter_by(gene_name = 'TFAP2A').first()
        self.assertEqual(subunit.uniprot_id, 'P05549')
        self.assertEqual(subunit.class_name, 'Basic helix-span-helix factors (bHSH)')
        binding = session.query(models.DNABindingDataset).filter_by(subunit_id = subunit.subunit_id).first()
        data = session.query(models.DNABindingData).filter_by(dataset_id = binding.dataset_id).first()
        self.assertEqual(data.position, 1)
        self.assertEqual(data.frequency_g, 185)

    def test_whoosh_indexing(self):
        for c in models.Compound.query.whoosh_search('2-Oxooctanoate').all():
            self.assertEqual(c.name, '2-Oxooctanoate')
            break

    def test_size(self):
        session = self.flk.session

        subunits = session.query(models.ProteinSubunit).all()
        self.assertGreater(len(subunits), 20000)


class LoadingTestFlaskCommonSchema(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.cs = flask_common_schema.FlaskCommonSchema(cache_dirname=self.cache_dirname,
                                clear_content = True, download_backups= False,
                                load_content = True, max_entries = 10,
                                verbose = True, test=True)

        for item in self.cs.text_indicies:
            flask_whooshalchemy.whoosh_index(self.cs.app, item)

    @classmethod
    def tearDownClass(self):
        models.db.session.remove()
        models.db.drop_all()
        shutil.rmtree(self.cache_dirname)



    def test_ncbi(self):
        session = self.cs.session
        tax = session.query(models.Taxon).filter_by(name = 'Homo sapiens').all()
        self.assertGreaterEqual(len(tax), 1)
        self.assertEqual(tax[0].ncbi_id, 9606)


        taxon = session.query(models.Taxon).filter_by(ncbi_id = 882).first()
        self.assertEqual(taxon.name, 'Desulfovibrio vulgaris str. Hildenborough')

    def test_pax(self):
        session = self.cs.session

        pax_compare = pax.Pax(cache_dirname = self.cache_dirname, download_backups = True, load_content = False)
        pax_session = pax_compare.session

        dataset = session.query(models.AbundanceDataSet).first()
        comparison = pax_session.query(pax.Dataset).filter_by(file_name = dataset.file_name).first()
        self.assertEqual(dataset.score, comparison.score)
        self.assertEqual(dataset.weight, comparison.weight)
        self.assertEqual(dataset.coverage, comparison.coverage)

        metadata = session.query(models.Metadata).filter_by(name = dataset.file_name).first()
        taxon = session.query(models._metadata_taxon).filter_by(_metadata_id = metadata.id).first()
        self.assertEqual(taxon.taxon_id, comparison.taxon_ncbi_id)

    def test_jaspar(self):
        session = self.cs.session
        subunit = session.query(models.ProteinSubunit).filter_by(gene_name = 'TFAP2A').first()
        self.assertEqual(subunit.uniprot_id, 'P05549')
        self.assertEqual(subunit.class_name, 'Basic helix-span-helix factors (bHSH)')
        binding = session.query(models.DNABindingDataset).filter_by(subunit_id = subunit.subunit_id).first()
        data = session.query(models.DNABindingData).filter_by(dataset_id = binding.dataset_id).first()
        self.assertEqual(data.position, 1)
        self.assertEqual(data.frequency_g, 185)

        metadata = session.query(models.Metadata).filter_by(name = 'RUNX1 Binding Motif').first()
        self.assertEqual(metadata.resource[0]._id, '8413232')

    def test_ecmdb(self):
        session = self.cs.session
        compound = session.query(models.Compound).filter_by(name = 'Deoxyuridine').first()
        self.assertEqual(compound.description, "2'-Deoxyuridine is a naturally occurring nucleoside. It is similar in chemical structure to uridine, but without the 2'-hydroxyl group.  It is considered to be an antimetabolite that is converted to deoxyuridine triphosphate during DNA synthesis.")
        structure = session.query(models.Structure).get(compound.structure_id)
        self.assertEqual(structure._value_inchi , 'InChI=1S/C9H12N2O5/c12-4-6-5(13)3-8(16-6)11-2-1-7(14)10-9(11)15/h1-2,5-6,8,12-13H,3-4H2,(H,10,14,15)/t5-,6+,8+/m0/s1')

    def test_arrayexpress(self):
        session = self.cs.session
        sample = session.query(models.RNASeqDataSet).filter_by(experiment_accession_number = 'E-MTAB-6272', sample_name='Sample 13').first()
        self.assertEqual(sample.assay, 'Sample 13')
        self.assertEqual(sample.ensembl_organism_strain, "mus_musculus")
        self.assertEqual(sample.read_type, "single")
        self.assertEqual(sample.full_strain_specificity, True)
        self.assertEqual(sample.reference_genome[0].download_url, "ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz")
        #print("here:{}".format(sample._metadata_id))
        self.assertEqual(len(sample._metadata.characteristic), 7)
        #print("here:{}".format(sample._metadata_id))
        self.assertEqual(len(sample._metadata.variable), 2)


        experiment = session.query(models.RNASeqExperiment).filter_by(accession_number = 'E-MTAB-6454').first()
        self.assertEqual(len(experiment.samples), 36)
        self.assertEqual(experiment.accession_number, 'E-MTAB-6454')
        self.assertEqual(experiment.exp_name, "Gene expression profiling across ontogenetic stages in wood white (Leptidea sinapis) reveals pathways linked to butterfly diapause regulation")
        self.assertEqual(experiment._experimentmetadata.description[:12], "In temperate")
        #self.assertEqual(experiment.has_fastq_files

    def test_sabio(self):
        session = self.cs.session
        compound = session.query(models.Compound).filter_by(compound_name = 'Peptide').first()
        self.assertEqual(compound._is_name_ambiguous , 1)
        structure = session.query(models.Structure).get(compound.structure_id)
        self.assertEqual(structure._value_smiles, '[*]C(N)C(=O)NC([*])C(O)=O')

        metadata = session.query(models.Metadata).filter_by(name = 'Kinetic Law 1').first()
        cell_line = session.query(models._metadata_cell_line).filter_by(_metadata_id = metadata.id).first()
        name = session.query(models.CellLine).get(cell_line.cell_line_id)
        self.assertEqual(name.name , 'variant DSAI (N76D/N87S/S103A/V104I)')


        q = session.query(models.KineticLaw) \
            .join((models.Metadata, models.KineticLaw._metadata)).join((models.Resource, models.Metadata.resource))\
            .filter(models.Resource.namespace == 'ec-code').filter(models.Resource._id.in_(['3.4.21.62']))

        compare = session.query(models.KineticLaw).filter_by(enzyme_id = q.first().enzyme_id).all()

        self.assertEqual(set([n.kinetic_law_id for n in q.all()]),
            set([c.kinetic_law_id for c in compare]))


        resource = session.query(models.Resource).filter_by(namespace = 'ec-code').filter_by(_id = '3.4.21.62').all()
        self.assertEqual(len(resource), 1)


    def test_corum(self):
        session = self.cs.session
        subunit = session.query(models.ProteinSubunit).filter_by(subunit_name = 'Histone deacetylase 5').first()
        self.assertEqual(subunit.uniprot_id, 'Q9UQL6')
        self.assertEqual(subunit.entrez_id, 10014)
        complex_  = session.query(models.ProteinComplex).get(subunit.proteincomplex_id)
        self.assertEqual(complex_.complex_name, 'BCL6-HDAC5 complex')

        metadata = session.query(models.Metadata).filter_by(name = 'BCL6-HDAC7 complex').first()
        resource = session.query(models._metadata_resource).filter_by(_metadata_id = metadata.id).first()
        pubmed = session.query(models.Resource).get(resource.resource_id)
        self.assertEqual(pubmed._id, '11929873')


    def test_intact_interactions(self):
        session = self.cs.session

        interact = session.query(models.ProteinInteractions).filter_by(participant_a = 'uniprotkb:P49418').all()
        self.assertEqual(set([c.site_b for c in interact]), set(['binding-associated region:1063-1070(MINT-376288)', '-']))

        for testcase in [interaction.protein_subunit for interaction in interact\
            if interaction.participant_b != 'uniprotkb:O43426']:
            self.assertEqual(len(testcase), 1)
            self.assertEqual(testcase[0].uniprot_id, 'P49418')

        for items in interact:
            self.assertTrue(items._metadata)

        interact = session.query(models.ProteinInteractions).filter_by(participant_a = 'intact:EBI-7121765').first()
        self.assertEqual(interact._metadata.resource[0]._id, '10542231|mint')

    def test_intact_complex_added(self):
        session = self.cs.session

        plex = session.query(models.ProteinComplex).filter_by(complex_name = 'CPAP-STIL complex').all()
        self.assertEqual(len(plex), 1)
        self.assertEqual(plex[0].go_id, '1905832|0110028|0046601|0005737|0008017|1900087|0120099')
        self.assertEqual(plex[0].su_cmt, 'Q7ZVT3(0)|Q8JGS1(0)|E7FCY1(0)')
        self.assertEqual(len(plex[0].protein_subunit), 3)

    def test_uniprot_added(self):
        session = self.cs.session

        uni = session.query(models.ProteinSubunit).filter_by(uniprot_id = 'Q72DQ8').first()
        self.assertEqual(uni.subunit_name, 'PYRH_DESVH')
        self.assertEqual(uni.entrez_id, 2795170)
        self.assertEqual(uni.length, 238)

    def test_whoosh(self):
        self.assertEqual(set([c.name for c in models.Compound.query.whoosh_search('adenine').all()]),
            set(['Adenosine', 'Adenosine monophosphate', 'Cyclic AMP', "Adenosine 3',5'-diphosphate", 'Adenine']))


# class DownloadLoadingTestFlaskCommonSchema(unittest.TestCase):
#
#     @classmethod
#     def setUpClass(self):
#         self.cache_dirname = tempfile.mkdtemp()
#
#         self.cs = flask_common_schema.FlaskCommonSchema(cache_dirname=self.cache_dirname,
#                                 download_backups= True, load_content = True, max_entries = 10,
#                                 verbose = True, test=True)
#
#         for item in self.cs.text_indicies:
#             flask_whooshalchemy.whoosh_index(self.cs.app, item)
#
#         self.session = self.cs.session
#
#
#     @classmethod
#     def tearDownClass(self):
#         models.db.session.remove()
#         models.db.drop_all()
#
#         shutil.rmtree(self.cache_dirname)
#
#     def test_progress_check(self):
#         for prog in self.session.query(models.Progress).all():
#             self.assertGreaterEqual(prog.amount_loaded, self.prev_load[prog.database_name])
