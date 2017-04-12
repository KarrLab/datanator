""" Command line utilities

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

from . import datanator
from .util import taxonomy_util
from cement.core.foundation import CementApp
from cement.core.controller import CementBaseController, expose
from pkg_resources import resource_filename
import openpyxl
import os.path


class BaseController(CementBaseController):

    class Meta:
        label = 'base'
        description = 'Utilities for aggregating data for biochemical models'


class GetKineticsController(CementBaseController):

    class Meta:
        label = 'get-kinetics'
        description = "Get kinetic data"
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            # we can make these arguments "positional" (required) by not adding "--arg" and just writing "arg"
            # the action 'store' will store the value passed for the option in self.app.pargs
            (['input_data_file'], dict(metavar='input-data-file',
                                       type=str, help="path to the input data spreadsheet (xlsx)")),
            (['output_data_file'], dict(metavar='output-data-file',
                                        type=str, help="path to the output data spreadsheet (xlsx)")),
            (['species'], dict(metavar='input-data-file',
                               type=str, help="name of the species you are searching")),
            # the following arguments are optional arguments, we connote this by adding "--" to the argument
            (['--min-temp'], dict(action='store',
                                  metavar='FLOAT', help="minimum temperature", default=15)),
            (['--max-temp'], dict(action='store',
                                  metavar='FLOAT', help="maximum temperature", default=40)),
            (['--min-ph'], dict(action='store',
                                metavar='FLOAT', help="minimum ph", default=5)),
            (['--max-ph'], dict(action='store',
                                metavar='FLOAT', help="maximum ph", default=9)),
            # the action 'store_true' turns this option into a boolean that will be store in self.app.pargs
            (['--include-mutants'], dict(action='store_true',
                                         help="include mutants")),
            (['--proxim-limit'], dict(action='store',
                                      metavar='FLOAT', help="the maximum acceptable taxonomic distance", default=1000))
        ]

        # print(ec_number_finder.format_ec_for_sabio("a"))

    @expose(hide=True)
    def default(self):
        self.app.log.info("In get_kinetics")
        self.app.log.info("Input file: '{}'".format(self.app.pargs.input_data_file))
        self.app.log.info("Output file: '{}'".format(self.app.pargs.output_data_file))
        self.app.log.info("Species: '{}'".format(self.app.pargs.species))
        self.app.log.info("Minimum temperature: '{}'".format(self.app.pargs.min_temp))
        self.app.log.info("Maximum temperature: '{}'".format(self.app.pargs.max_temp))
        self.app.log.info("Minimum ph: '{}'".format(self.app.pargs.min_ph))
        self.app.log.info("Maximum ph: '{}'".format(self.app.pargs.max_ph))
        self.app.log.info("Include mutants: '{}'".format(self.app.pargs.include_mutants))
        self.app.log.info("Proximity Limit: '{}'".format(self.app.pargs.proxim_limit))

        # temp_range = [30, 40], enzyme_type = "wildtype", ph_range = [5,9], proxim_limit=1000):

        enzyme_type = "wildtype"
        if self.app.pargs.include_mutants:
            enzyme_type = "(wildtype OR mutant)"
        datanator.get_kinetic_data(self.app.pargs.input_data_file,
                                   self.app.pargs.output_data_file,
                                   self.app.pargs.species,
                                   temp_range=[float(self.app.pargs.min_temp), float(self.app.pargs.max_temp)],
                                   ph_range=[float(self.app.pargs.min_ph), float(self.app.pargs.max_ph)],
                                   enzyme_type=enzyme_type,
                                   proxim_limit=float(self.app.pargs.proxim_limit)
                                   )


class GenerateTemplateController(CementBaseController):

    class Meta:
        label = 'generate-template'
        description = "Generate an Excel template for specifying which reactions to aggregate kinetic data about"
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = []  # todo: add argument for template output file

    @expose(hide=True)
    def default(self):
        # todo: generate template
        template_filename = resource_filename('kinetic_datanator', 'data/InputTemplate.xlsx')
        output_filename = 'InputTemplate.xlsx'
        wb = openpyxl.load_workbook(filename)
        wb.save(output_filename)


class TaxonomyController(CementBaseController):

    class Meta:
        label = 'taxonomy'
        description = 'Get taxonomic information'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = []


class TaxonomyGetRankController(CementBaseController):

    class Meta:
        label = 'get-rank'
        description = 'Get the rank of a taxon'
        stacked_on = 'taxonomy'
        stacked_type = 'nested'
        arguments = [
            (['taxon_id_or_name'], dict(metavar='taxon_id_or_name',
                                        type=str, help="Taxon id or name (examples: 9606, 'Homo sapiens')"))
        ]

    @expose(hide=True)
    def default(self):
        taxon_id_or_name = self.app.pargs.taxon_id_or_name
        taxon = taxonomy_util.Taxon(taxon_id_or_name)
        if taxon.distance_from_nearest_ncbi_taxon != 0:
            raise ValueError('The NCBI taxonomy database does not contain a taxon with id or name {}'.format(taxon_id_or_name))
        print(taxon.get_rank())


class TaxonomyGetParentsController(CementBaseController):

    class Meta:
        label = 'get-parents'
        description = 'Get the parents of a taxon'
        stacked_on = 'taxonomy'
        stacked_type = 'nested'
        arguments = [
            (['taxon_id_or_name'], dict(metavar='taxon_id_or_name',
                                        type=str, help="Taxon id or name (examples: 9606, 'Homo sapiens')"))
        ]

    @expose(hide=True)
    def default(self):
        taxon_id_or_name = self.app.pargs.taxon_id_or_name
        taxon = taxonomy_util.Taxon(taxon_id_or_name)
        if taxon.id_of_nearest_ncbi_taxon is None:
            raise ValueError('The NCBI taxonomy database does not contain a taxon with id or name {}'.format(taxon_id_or_name))

        parents = taxon.get_parent_taxa()
        for i_parent, parent in enumerate(parents):
            print(parent.name)


class TaxonomyGetCommonAncestorController(CementBaseController):

    class Meta:
        label = 'get-common-ancestor'
        description = "Get the latest common ancestor between two taxa"
        stacked_on = 'taxonomy'
        stacked_type = 'nested'
        arguments = [
            (['taxon_id_or_name_1'], dict(metavar='taxon_id_or_name_1',
                                          type=str, help="Taxon id or name (examples: 9606, 'Homo sapiens')")),
            (['taxon_id_or_name_2'], dict(metavar='taxon_id_or_name_2',
                                          type=str, help="Taxon id or name (examples: 9606, 'Homo sapiens')")),
        ]

    @expose(hide=True)
    def default(self):
        taxon_id_or_name_1 = self.app.pargs.taxon_id_or_name_1
        taxon_id_or_name_2 = self.app.pargs.taxon_id_or_name_2
        taxon_1 = taxonomy_util.Taxon(taxon_id_or_name_1)
        taxon_2 = taxonomy_util.Taxon(taxon_id_or_name_2)

        if taxon_1.id_of_nearest_ncbi_taxon is None:
            raise ValueError('The NCBI taxonomy database does not contain a taxon with id or name {}'.format(taxon_id_or_name_1))
        if taxon_2.id_of_nearest_ncbi_taxon is None:
            raise ValueError('The NCBI taxonomy database does not contain a taxon with id or name {}'.format(taxon_id_or_name_2))

        print(taxon_1.get_common_ancestor(taxon_2).name)


class TaxonomyGetDistanceToCommonAncestorController(CementBaseController):

    class Meta:
        label = 'get-distance-to-common-ancestor'
        description = "Get the distance to the latest common ancestor between two taxa"
        stacked_on = 'taxonomy'
        stacked_type = 'nested'
        arguments = [
            (['taxon_id_or_name_1'], dict(metavar='taxon_id_or_name_1',
                                          type=str, help="Taxon id or name (examples: 9606, 'Homo sapiens')")),
            (['taxon_id_or_name_2'], dict(metavar='taxon_id_or_name_2',
                                          type=str, help="Taxon id or name (examples: 9606, 'Homo sapiens')")),
        ]

    @expose(hide=True)
    def default(self):
        taxon_id_or_name_1 = self.app.pargs.taxon_id_or_name_1
        taxon_id_or_name_2 = self.app.pargs.taxon_id_or_name_2
        taxon_1 = taxonomy_util.Taxon(taxon_id_or_name_1)
        taxon_2 = taxonomy_util.Taxon(taxon_id_or_name_2)

        if taxon_1.id_of_nearest_ncbi_taxon is None:
            raise ValueError('The NCBI taxonomy database does not contain a taxon with id or name {}'.format(taxon_id_or_name_1))
        if taxon_2.id_of_nearest_ncbi_taxon is None:
            raise ValueError('The NCBI taxonomy database does not contain a taxon with id or name {}'.format(taxon_id_or_name_2))

        print(taxon_1.get_distance_to_common_ancestor(taxon_2))


class TaxonomyGetDistanceToRoot(CementBaseController):

    class Meta:
        label = 'get-distance-to-root'
        description = "Get the distance to from a taxon to the root of the taxonomic tree"
        stacked_on = 'taxonomy'
        stacked_type = 'nested'
        arguments = [
            (['taxon_id_or_name'], dict(metavar='taxon_id_or_name',
                                        type=str, help="Taxon id or name (examples: 9606, 'Homo sapiens')")),
        ]

    @expose(hide=True)
    def default(self):
        taxon_id_or_name = self.app.pargs.taxon_id_or_name
        taxon = taxonomy_util.Taxon(taxon_id_or_name)

        if taxon.id_of_nearest_ncbi_taxon is None:
            raise ValueError('The NCBI taxonomy database does not contain a taxon with id or name {}'.format(taxon_id_or_name))

        print(taxon.get_distance_to_root())


class App(CementApp):

    class Meta:
        label = "kinetic_datanator"
        base_controller = "base"
        handlers = [
            BaseController,

            GetKineticsController,
            GenerateTemplateController,

            TaxonomyController,
            TaxonomyGetRankController,
            TaxonomyGetParentsController,
            TaxonomyGetCommonAncestorController,
            TaxonomyGetDistanceToCommonAncestorController,
            TaxonomyGetDistanceToRoot,
        ]


def main():
    with App() as app:
        app.run()

if __name__ == "__main__":
    main()
