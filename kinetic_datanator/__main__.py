""" Command line utilities

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-04-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

# from kinetic_datanator import datanator #todo
from __future__ import print_function
from cement.core import controller
from cement.core import foundation
from kinetic_datanator import io
from kinetic_datanator.core import data_model, common_schema, upload_data, data_query#, json_schema
from kinetic_datanator.core.render_form import render_html_from_schema
from kinetic_datanator.data_source import *
from kinetic_datanator.data_source import refseq
from kinetic_datanator.util import molecule_util
from kinetic_datanator.util import taxonomy_util
from pkg_resources import resource_filename
import bioservices
import kinetic_datanator
import pubchempy
import re
import shutil
import sys
import os
from Bio import SeqIO
from kinetic_datanator.api.query import reaction_kinetics
#from kinetic_datanator.core import data_query


CACHE_DIRNAME = os.path.join(os.path.dirname(__file__), 'data_source', 'cache')


class BaseController(controller.CementBaseController):

    class Meta:
        label = 'base'
        description = 'Utilities for aggregating data for biochemical models'
        arguments = [
            (['-v', '--version'], dict(action='version', version=kinetic_datanator.__version__)),
        ]

    @controller.expose(hide=True)
    def default(self):
        self.app.args.print_help()


class UploadDataController(controller.CementBaseController):

    class Meta:
        label = 'upload'
        description = "Upload reference genome into Datanator. The reference genome must be in Genbank or Ensembl format"
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = []

    @controller.expose(hide=True)
    def default(self):
        self.app.args.print_help()


class UploadReferenceGenome(controller.CementBaseController):

    class Meta:
        label = 'reference-genome'
        description = 'Upload a ref-seq genome'
        stacked_on = 'upload'
        stacked_type = 'nested'
        arguments = [
            (['path_to_annotation_file'], dict(type=str, help="path to reference genome file", default=CACHE_DIRNAME)),
            (['--path_to_database'], dict(type=str, help="path to build the database", default=CACHE_DIRNAME)),
        ]

    @controller.expose(hide=True)
    def default(self):
        pargs = self.app.pargs
        print(pargs.__dict__)
        #bio_seqio_object = SeqIO.parse(pargs.path_to_annotation_file, "genbank")
        #list_of_bio_seqio_objects = [bio_seqio_object]
        upload_data.UploadData(cache_dirname=pargs.path_to_database).upload_reference_genome(pargs.path_to_annotation_file)
        # refseq.Refseq(cache_dirname=pargs.path_to_database).load_content(list_of_bio_seqio_objects)


class UploadRNASeqExperiment(controller.CementBaseController):

    class Meta:
        label = 'rna-seq-experiment'
        description = 'Upload an RNA-seq experiment'
        stacked_on = 'upload'
        stacked_type = 'nested'
        arguments = [
            (['path_to_folder with excel files'], dict(type=str, help="path to reference genome file", default=CACHE_DIRNAME)),
            (['--path_to_database'], dict(type=str, help="path to build the database", default=CACHE_DIRNAME)),
        ]

    @controller.expose(hide=True)
    def default(self):
        pargs = self.app.pargs
        print(pargs.__dict__)
        bio_seqio_object = SeqIO.parse(pargs.path_to_annotation_file, "genbank")
        list_of_bio_seqio_objects = [bio_seqio_object]
        refseq.Refseq(cache_dirname=pargs.path_to_database).load_content(list_of_bio_seqio_objects)

class UploadData(controller.CementBaseController):

    class Meta:
        label = 'general-data'
        description = 'Upload data from overall schema'
        stacked_on = 'upload'
        stacked_type = 'nested'
        arguments = [
            (['data_type'], dict(type=str, help="type of data to be uploaded")),
        ]

    @controller.expose(hide=True)
    def default(self):
        pargs = self.app.pargs
        data_type = pargs.data_type
        
        
        print(data_type)
        a_json_schema = json_schema.get_json_schema(data_type)

        the_json = render_html_from_schema.RenderForms().render(a_json_schema)
        

        #the_json = """{"_experimentmetadata":[],"accession_number":"Test-Accesion","exp_name":"test_experiment","has_fastq_files":false,"samples":[{"_metadata":[],"assay":"blue","ensembl_organism_strain":"green","experiment_accession_number":"test","full_strain_specificity":false,"reference_genome":[]}]}"""
        print(the_json)


class BuildController(controller.CementBaseController):

    class Meta:
        label = 'build'
        description = "Build aggregated database"
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['--max-entries'], dict(type=int, help="number of normalized entries to add per database. Default: Full Database", default=float('inf'))),
            (['--path'], dict(type=str, help="path to build the database", default=CACHE_DIRNAME)),
            (['--clear-existing-content'], dict(type=bool, help="clears existing content of the db if exists", default=False)),
            (['--verbose'], dict(type=bool, help="verbosity", default=CACHE_DIRNAME))
        ]

    @controller.expose(help='Builds Corum Complex DB from source')
    def corum(self):
        pargs = self.app.pargs
        corum.Corum(cache_dirname=pargs.path, load_content=True, download_backups=False,
                    max_entries=pargs.max_entries, verbose=pargs.verbose)

    @controller.expose(help='Builds IntAct Interactions and Complex DB from source')
    def intact(self):
        pargs = self.app.pargs
        intact.IntAct(cache_dirname=pargs.path, load_content=True, download_backups=False,
                      max_entries=pargs.max_entries, verbose=pargs.verbose)

    @controller.expose(help='Builds Sabio Reaction Kinetics DB from source')
    def sabio(self):
        pargs = self.app.pargs
        sabio_rk.SabioRk(cache_dirname=pargs.path, load_content=True, download_backups=False,
                         max_entries=pargs.max_entries, verbose=pargs.verbose)

    @controller.expose(help='Builds Pax Protein Abundance DB from source')
    def pax(self):
        pargs = self.app.pargs
        pax.Pax(cache_dirname=pargs.path, load_content=True, download_backups=False, max_entries=pargs.max_entries, verbose=pargs.verbose)

    @controller.expose(help='Builds Array Express RNA Seq DB from source')
    def array_express(self):
        pargs = self.app.pargs
        array_express.ArrayExpress(cache_dirname=pargs.path, load_content=True, download_backups=False,
                                   max_entries=pargs.max_entries, verbose=pargs.verbose)

    @controller.expose(help='Builds Jaspar DNA protein interaction DB from source')
    def jaspar(self):
        pargs = self.app.pargs
        jaspar.Jaspar(cache_dirname=pargs.path, load_content=True, download_backups=False,
                      max_entries=pargs.max_entries, verbose=pargs.verbose)

    @controller.expose(help='Builds Uniprot Protein DB from source')
    def uniprot(self):
        pargs = self.app.pargs
        uniprot.Uniprot(cache_dirname=pargs.path, load_content=True, download_backups=False,
                        max_entries=pargs.max_entries, verbose=pargs.verbose)

    @controller.expose(help='Builds ECMDB metabolite DB from source')
    def ecmdb(self):
        pargs = self.app.pargs
        ecmdb.Ecmdb(cache_dirname=pargs.path, load_content=True, download_backups=False,
                    max_entries=pargs.max_entries, verbose=pargs.verbose)

    @controller.expose(hide=True)
    def default(self):
        self.app.args.print_help()


class AggregateBuildController(controller.CementBaseController):

    class Meta:
        label = 'aggregate'
        description = "Builds Aggregated Database"
        stacked_on = 'build'
        stacked_type = 'nested'
        arguments = [
            (['--path'], dict(type=str, help="path to build the database", default=CACHE_DIRNAME)),
            (['--verbose'], dict(type=bool, help="verbosity", default=CACHE_DIRNAME)),
            (['--max-entries'], dict(type=int, help="number of normalized entries to add per database. Default: Full Database", default=float('inf'))),
            (['--build-on-existing'], dict(type=bool, help="load from existing database on karr lab server", default=False)),
            (['--load-full-small-dbs'], dict(type=bool, help="loads entire small database modules", default=True))
        ]

    @controller.expose(help='Controller that controls aggregated')
    def default(self):
        pargs = self.app.pargs
        common_schema.CommonSchema(cache_dirname=pargs.path, load_content=True,
                                              download_backups=False, max_entries=pargs.max_entries, verbose=pargs.verbose)


class DownloadController(controller.CementBaseController):

    class Meta:
        label = 'download'
        description = "Download existing databases from Karr Lab Server"
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['--path'], dict(type=str, help="path to download the modules", default=CACHE_DIRNAME))
        ]

    @controller.expose(help='Loads Corum Complex DB from Karr Lab Server')
    def corum(self):
        pargs = self.app.pargs
        corum.Corum(cache_dirname=pargs.path, download_backups=True)

    @controller.expose(help='Loads IntAct Interactions and Complex DB from Karr Lab Server')
    def intact(self):
        pargs = self.app.pargs
        intact.IntAct(cache_dirname=pargs.path, download_backups=True)

    @controller.expose(help='Loads Sabio Reaction Kinetics DB from Karr Lab Server')
    def sabio(self):
        pargs = self.app.pargs
        sabio_rk.SabioRk(cache_dirname=pargs.path, download_backups=True)

    @controller.expose(help='Loads Pax Protein Abundance DB from Karr Lab Server')
    def pax(self):
        pargs = self.app.pargs
        pax.Pax(cache_dirname=pargs.path, download_backups=True)

    @controller.expose(help='Loads Array Express RNA Seq DB from Karr Lab Server')
    def array_express(self):
        pargs = self.app.pargs
        array_express.ArrayExpress(cache_dirname=pargs.path, download_backups=True)

    @controller.expose(help='Loads Jaspar DNA protein interaction DB from Karr Lab Server')
    def jaspar(self):
        pargs = self.app.pargs
        jaspar.Jaspar(cache_dirname=pargs.path, download_backups=True)

    @controller.expose(help='Loads Uniprot Protein DB from Karr Lab Server')
    def uniprot(self):
        pargs = self.app.pargs
        uniprot.Uniprot(cache_dirname=pargs.path, download_backups=True)

    @controller.expose(help='Loads ECMDB metabolite DB from Karr Lab Server')
    def ecmdb(self):
        pargs = self.app.pargs
        ecmdb.Ecmdb(cache_dirname=pargs.path, download_backups=True)

    @controller.expose(help='Loads Aggregated DB from Karr Lab Server')
    def aggregate(self):
        pargs = self.app.pargs
        common_schema.CommonSchema(cache_dirname=pargs.path, download_backups=True)

    @controller.expose(hide=True)
    def default(self):
        self.app.args.print_help()


class GetDataController(controller.CementBaseController):

    class Meta:
        label = 'get-data'
        description = "Get relevant data for a model of a taxon"
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['input_file'], dict(type=str, help="path to the input data spreadsheet (.xlsx)")),
            (['output_file'], dict(type=str, help="path to the output data spreadsheet (.xlsx)")),
            (['--max-taxon-dist'], dict(help="Maximum acceptable taxonomic distance",
                                        type=int, default=None)),
            (['--taxon-dist-scale'], dict(help="Exponential constant for scoring taxonomic distance",
                                          type=float, default=float('nan'))),
            (['--include-variants'], dict(help="If set, include data from genetic variants",
                                          action='store_true', default=False)),
            (['--temperature'], dict(help="Target temperature", type=float, default=37.)),
            (['--temperature-std'], dict(help="Standard deviation for scoring temperatures",
                                         type=float, default=1.)),
            (['--ph'], dict(help="Target pH", type=float, default=7.5)),
            (['--ph-std'], dict(help="Standard deviation for scoring pHs",
                                type=float, default=0.3)),
        ]

    @controller.expose(hide=True)
    def default(self):
        pargs = self.app.pargs
        genetics, compartments, species, reactions = io.InputReader().run(pargs.input_file)
        #print(reactions)
        #for thing in reactions[0].get_reactants():
        #    print(thing.specie.to_inchi())
        #    print(thing.specie.to_smiles())
        cache_dirname = '{}/data_source/cache'.format(os.path.dirname(__file__))
        #print(cache_dirname)
        observed_values = reaction_kinetics.ReactionKineticsQuery(cache_dirname=cache_dirname).get_observed_result(reactions[0])

        class ConcreteDataQueryGenerator(data_query.DataQueryGenerator):

            def get_observed_result(self):
                pass
        print("the length is initially {}".format(len(observed_values)))
        gen = ConcreteDataQueryGenerator()
        gen.filters = [
            data_query.TemperatureRangeFilter(min=36., max=38.),
            data_query.TaxonomicDistanceFilter(taxon='Mycoplasma genitalium')#, max=7)
        ]
        filtered = gen.filter_observed_results(None, observed_values)
        print(filtered.observed_value_indices)
        print(data_query.ConsensusGenerator().calc_average(filtered.observed_value_indices, method='median'))
            #print(thing.)
        
        #con_gen =  data_query.ConsensusGenerator().run(filtered, "median")
        #print(con_gen)
        
        #print(gen.get_consensus(None, filtered))
        #print("the length after filtering is {}".format(len(filtered.observed_results)))
        #consensus = data_query.ConsensusGenerator().run(filtered, 'median')
        #print(consensus)


        # todo: implement
        return

        data = datanator.get_kinetic_data(
            taxon=genetics.taxon, max_taxon_dist=pargs.max_taxon_dist, taxon_dist_scale=pargs.taxon_dist_scale,
            include_variants=pargs.include_variants,
            temperature=pargs.temperature, temperature_std=pargs.temperature_std,
            ph=pargs.ph, ph_std=pargs.ph_std)
        io.OutputWriter().run(pargs.pargs.output_file, data)


class GenerateTemplateController(controller.CementBaseController):

    class Meta:
        label = 'generate-template'
        description = "Generate an Excel template for specifying which reactions to aggregate kinetic data about"
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['filename'], dict(type=str, help="path to save the Excel template")),
        ]

    @controller.expose(hide=True)
    def default(self):
        # todo: generate template
        template_filename = resource_filename('kinetic_datanator', 'data/InputTemplate.xlsx')
        shutil.copyfile(template_filename, self.app.pargs.filename)


class GenerateRNASeqTemplate(controller.CementBaseController):

    class Meta:
        label = 'generate-rna-seq-template'
        description = "Generate a folder with excel tables to upload rna-seq experiments"
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['directory_name'], dict(type=str, help="path to save the directory")),
        ]

    @controller.expose(hide=True)
    def default(self):
        # todo: generate template
        template_directory = resource_filename('kinetic_datanator', 'data/RNA-Seq_Experiment_Template')
        shutil.copytree(template_directory, "{}/RNA-Seq_Experiment_Template".format(self.app.pargs.directory_name))


class TaxonomyController(controller.CementBaseController):

    class Meta:
        label = 'taxonomy'
        description = 'Taxonomy utilities'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = []

    @controller.expose(hide=True)
    def default(self):
        self.app.args.print_help()


class TaxonomyGetRankController(controller.CementBaseController):

    class Meta:
        label = 'get-rank'
        description = 'Get the rank of a taxon'
        stacked_on = 'taxonomy'
        stacked_type = 'nested'
        arguments = [
            (['taxon_id_or_name'], dict(type=str, help="Taxon id or name (examples: 9606, 'Homo sapiens')"))
        ]

    @controller.expose(hide=True)
    def default(self):
        taxon = create_taxon(self.app.pargs.taxon_id_or_name)
        if taxon.distance_from_nearest_ncbi_taxon != 0:
            raise ValueError('The NCBI taxonomy database does not contain a taxon with id or name {}'
                             .format(self.app.pargs.taxon_id_or_name))
        print(taxon.get_rank())


class TaxonomyGetParentsController(controller.CementBaseController):

    class Meta:
        label = 'get-parents'
        description = 'Get the parents of a taxon'
        stacked_on = 'taxonomy'
        stacked_type = 'nested'
        arguments = [
            (['taxon_id_or_name'], dict(type=str, help="Taxon id or name (examples: 9606, 'Homo sapiens')"))
        ]

    @controller.expose(hide=True)
    def default(self):
        taxon = create_taxon(self.app.pargs.taxon_id_or_name)
        if taxon.id_of_nearest_ncbi_taxon is None:
            raise ValueError('The NCBI taxonomy database does not contain a taxon with id or name {}'
                             .format(self.app.pargs.taxon_id_or_name))

        parents = taxon.get_parent_taxa()
        for i_parent, parent in enumerate(parents):
            print(parent.name)


class TaxonomyGetCommonAncestorController(controller.CementBaseController):

    class Meta:
        label = 'get-common-ancestor'
        description = "Get the latest common ancestor between two taxa"
        stacked_on = 'taxonomy'
        stacked_type = 'nested'
        arguments = [
            (['taxon_id_or_name_1'], dict(type=str, help="Taxon id or name (examples: 9606, 'Homo sapiens')")),
            (['taxon_id_or_name_2'], dict(type=str, help="Taxon id or name (examples: 9606, 'Homo sapiens')")),
        ]

    @controller.expose(hide=True)
    def default(self):
        taxon_1 = create_taxon(self.app.pargs.taxon_id_or_name_1)
        taxon_2 = create_taxon(self.app.pargs.taxon_id_or_name_2)

        if taxon_1.id_of_nearest_ncbi_taxon is None:
            raise ValueError('The NCBI taxonomy database does not contain a taxon with id or name {}'.format(
                self.app.pargs.taxon_id_or_name_1))
        if taxon_2.id_of_nearest_ncbi_taxon is None:
            raise ValueError('The NCBI taxonomy database does not contain a taxon with id or name {}'.format(
                self.app.pargs.taxon_id_or_name_2))

        print(taxon_1.get_common_ancestor(taxon_2).name)


class TaxonomyGetDistanceToCommonAncestorController(controller.CementBaseController):

    class Meta:
        label = 'get-distance-to-common-ancestor'
        description = "Get the distance to the latest common ancestor between two taxa"
        stacked_on = 'taxonomy'
        stacked_type = 'nested'
        arguments = [
            (['taxon_id_or_name_1'], dict(type=str, help="Taxon id or name (examples: 9606, 'Homo sapiens')")),
            (['taxon_id_or_name_2'], dict(type=str, help="Taxon id or name (examples: 9606, 'Homo sapiens')")),
        ]

    @controller.expose(hide=True)
    def default(self):
        taxon_1 = create_taxon(self.app.pargs.taxon_id_or_name_1)
        taxon_2 = create_taxon(self.app.pargs.taxon_id_or_name_2)

        if taxon_1.id_of_nearest_ncbi_taxon is None:
            raise ValueError('The NCBI taxonomy database does not contain a taxon with id or name {}'.format(
                self.app.pargs.taxon_id_or_name_1))
        if taxon_2.id_of_nearest_ncbi_taxon is None:
            raise ValueError('The NCBI taxonomy database does not contain a taxon with id or name {}'.format(
                self.app.pargs.taxon_id_or_name_2))

        print(taxon_1.get_distance_to_common_ancestor(taxon_2))


class TaxonomyGetDistanceToRoot(controller.CementBaseController):

    class Meta:
        label = 'get-distance-to-root'
        description = "Get the distance to from a taxon to the root of the taxonomic tree"
        stacked_on = 'taxonomy'
        stacked_type = 'nested'
        arguments = [
            (['taxon_id_or_name'], dict(type=str, help="Taxon id or name (examples: 9606, 'Homo sapiens')")),
        ]

    @controller.expose(hide=True)
    def default(self):
        taxon = create_taxon(self.app.pargs.taxon_id_or_name)

        if taxon.id_of_nearest_ncbi_taxon is None:
            raise ValueError('The NCBI taxonomy database does not contain a taxon with id or name {}'
                             .format(self.app.pargs.taxon_id_or_name))

        print(taxon.get_distance_to_root())


class MoleculeController(controller.CementBaseController):

    class Meta:
        label = 'molecule'
        description = 'Molecule utilities'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = []

    @controller.expose(hide=True)
    def default(self):
        self.app.args.print_help()


class MoleculeGetStructureController(controller.CementBaseController):

    class Meta:
        label = 'get-structure'
        description = 'Get the structure of a molecule by its name or id'
        stacked_on = 'molecule'
        stacked_type = 'nested'
        arguments = [
            (['--by-name'], dict(dest='by_name', action='store_true', default=True,
                                 help="If set, lookup structure by name")),
            (['--by-id'], dict(dest='by_name', action='store_false', default=True,
                               help="If set, lookup structure by id")),
            (['--namespace'], dict(type=str, help="Namespace of id")),
            (['name_or_id'], dict(type=str, help="Name or id of the molecule")),
        ]

    @controller.expose(hide=True)
    def default(self):
        if self.app.pargs.by_name:
            compounds = pubchempy.get_compounds(self.app.pargs.name_or_id, 'name')
            results = [[compound.synonyms[0], 'pubchem.compound', compound.cid, compound.inchi] for compound in compounds]
        else:
            unichem = bioservices.UniChem()
            structure = unichem.get_structure(int(float(self.app.pargs.name_or_id)), self.app.pargs.namespace)
            if structure:
                results = [['', self.app.pargs.namespace, self.app.pargs.name_or_id, structure['standardinchi']]]
            else:
                results = []

        if not results:
            print('Unable to find structure', file=sys.stderr)
            return

        lens = [
            max(4, max(len(str(r[0])) for r in results)),
            max(9, max(len(str(r[1])) for r in results)),
            max(2, max(len(str(r[2])) for r in results)),
            max(9, max(len(str(r[3])) for r in results)),
        ]
        format = '{{:<{}}}  {{:<{}}}  {{:<{}}}  {{:<{}}}'.format(*lens)
        print(format.format('Name', 'Namespace', 'ID', 'Structure'))
        print(format.format('=' * lens[0], '=' * lens[1], '=' * lens[2], '=' * lens[3]))
        for result in results:
            print(format.format(*result))


class MoleculeConvertStructureController(controller.CementBaseController):

    class Meta:
        label = 'convert-structure'
        description = 'Convert molecule structure'
        stacked_on = 'molecule'
        stacked_type = 'nested'
        arguments = [
            (['structure'], dict(type=str, help="Structure in InChI, MOL, or canonical SMILES format")),
            (['format'], dict(type=str, help="Output format: inchi (InChI), mol (MOL), or can (canonical SMILES)")),
        ]

    @controller.expose(hide=True)
    def default(self):
        structure = self.app.pargs.structure
        format = self.app.pargs.format
        print(molecule_util.Molecule(structure=structure).to_format(format))


class ReactionController(controller.CementBaseController):

    class Meta:
        label = 'reaction'
        description = 'Reaction utilities'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = []

    @controller.expose(hide=True)
    def default(self):
        self.app.args.print_help()


class ReactionGetEcNumberController(controller.CementBaseController):

    class Meta:
        label = 'get-ec-number'
        description = 'Use Ezyme to predict the EC number of a reaction'
        stacked_on = 'reaction'
        stacked_type = 'nested'
        arguments = [
            (['reaction'], dict(
                type=str, help=("The reaction (e.g. {ATP(4-) + H2O --> ADP(3-) + PI(2-) + proton}, "
                                "{inchi-atp + inchi-h2o --> inchi-adp + inchi-pi + inchi-h}), or "
                                "{smiles-atp + smiles-h2o --> smiles-adp + smiles-pi + smiles-h})"
                                ))),
        ]

    # todo: add find_ec
    @controller.expose(hide=True)
    def default(self):
        # parse input
        def parse_participants(side, coefficient, reaction, errors):
            for participant in side.split(' + '):
                participant = participant.strip()

                pubchem_compounds = pubchempy.get_compounds(participant, 'name')
                if len(pubchem_compounds) == 1:
                    structure = pubchem_compounds[0].inchi
                elif molecule_util.Molecule(structure=participant).get_format():
                    structure = participant
                else:
                    structure = ''
                    errors.append(participant)

                reaction.participants.append(data_model.ReactionParticipant(
                    specie=data_model.Specie(structure=structure),
                    coefficient=coefficient,
                ))

        match = re.match('^(.*?)<*[\-=]{1,2}>(.*?)$', self.app.pargs.reaction)
        if not match:
            print('The reaction is ill-formed', file=sys.stderr)
            return

        reaction = data_model.Reaction()
        errors = []
        parse_participants(match.group(1), -1, reaction, errors)
        parse_participants(match.group(2),  1, reaction, errors)
        if errors:
            print('Unable to interpret participants:\n  '.format('\n  '.join(errors)), file=sys.stderr)
            return

        # predict EC number
        results = ezyme.Ezyme().run(reaction)

        # print results
        if results:
            print('EC number  Weight')
            print('=========  ======')
            for result in results:
                print('{:<9s} {:>6.2f}'.format(result.ec_number, result.score))
        else:
            print('There is no appropriate EC number for the reaction')


class App(foundation.CementApp):

    class Meta:
        label = "kinetic_datanator"
        base_controller = "base"
        handlers = [
            BaseController,

            UploadDataController,
            UploadReferenceGenome,
            UploadRNASeqExperiment,
            UploadData,

            DownloadController,
            BuildController,
            AggregateBuildController,

            GetDataController,
            GenerateTemplateController,
            GenerateRNASeqTemplate,

            # TaxonomyController2,
            # TaxonomyGetRankController2,

            TaxonomyController,
            TaxonomyGetRankController,
            TaxonomyGetParentsController,
            TaxonomyGetCommonAncestorController,
            TaxonomyGetDistanceToCommonAncestorController,
            TaxonomyGetDistanceToRoot,

            MoleculeController,
            MoleculeGetStructureController,
            MoleculeConvertStructureController,

            ReactionController,
            ReactionGetEcNumberController,
        ]


def create_taxon(id_or_name):
    """ Create a taxon with NCBI id=:obj:`id_or_name` or name=:obj:`id_or_name`

    Args:
        id_or_name (:obj:`str`): NCBI id or name

    Returns:
        :obj:`taxonomy_util.Taxon`: taxon
    """
    ncbi_id = None
    name = None
    try:
        ncbi_id = float(id_or_name)
    except ValueError:
        name = id_or_name

    return taxonomy_util.Taxon(ncbi_id=ncbi_id, name=name)


def main():
    with App() as app:
        app.run()


if __name__ == "__main__":
    main()
