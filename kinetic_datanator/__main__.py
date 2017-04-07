""" 
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-06
:Copyright: 2017, Karr Lab
:License: MIT
"""

from . import datanator
from . import taxon_finder
from cement.core.foundation import CementApp
from cement.core.controller import CementBaseController, expose
from pkg_resources import resource_filename
import openpyxl
import os.path


class BaseController(CementBaseController):
	class Meta:
		label = 'base'
		description = 'This app find information about reactions'


class GetKineticsController(CementBaseController):
	class Meta:
		label = 'get-kinetics'
		description = "This command gets kinetics"
		stacked_on = 'base'
		stacked_type = 'nested'
		arguments = [
		#we can make these arguments "positional" (required) by not adding "--arg" and just writing "arg"
		#the action 'store' will store the value passed for the option in self.app.pargs
		(['input_data_file'],dict(metavar='input-data-file', 
			type=str, help="path to the input data spreadsheet (xlsx)")),
		(['output_data_file'],dict(metavar='output-data-file', 
			type=str, help="path to the output data spreadsheet (xlsx)")),
		(['species'],dict(metavar='input-data-file', 
			type=str, help="name of the species you are searching")),
		#the following arguments are optional arguments, we connote this by adding "--" to the argument
		(['--min-temp'], dict(action='store', 
			metavar='FLOAT', help = "minimum temperature", default=15)),
		(['--max-temp'], dict(action='store', 
			metavar='FLOAT', help = "maximum temperature", default=40)),
		(['--min-ph'], dict(action='store', 
			metavar='FLOAT', help = "minimum ph",default=5)),
		(['--max-ph'], dict(action='store', 
			metavar='FLOAT', help = "maximum ph", default=9)),
		# the action 'store_true' turns this option into a boolean that will be store in self.app.pargs
		(['--include-mutants'], dict(action='store_true', 
			help = "include mutants")),
		(['--proxim-limit'], dict(action='store', 
			metavar='FLOAT', help = "the maximum acceptable taxonomic distance", default=1000))
		]



		#print(ec_number_finder.format_ec_for_sabio("a"))
		

	@expose(help="This command gets kinetics", hide=True)
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

#temp_range = [30, 40], enzyme_type = "wildtype", ph_range = [5,9], proxim_limit=1000):

		enzyme_type = "wildtype"
		if self.app.pargs.include_mutants:
			enzyme_type = "(wildtype OR mutant)"
		datanator.get_kinetic_data(self.app.pargs.input_data_file, self.app.pargs.output_data_file,
		 self.app.pargs.species, temp_range=[float(self.app.pargs.min_temp), float(self.app.pargs.max_temp)], ph_range = [float(self.app.pargs.min_ph), float(self.app.pargs.max_ph)],
		 enzyme_type = enzyme_type, proxim_limit = float(self.app.pargs.proxim_limit))



class GenerateTemplateController(CementBaseController):
	class Meta:
		label = 'generate-template'
		description = "This command generates a template"
		stacked_on = 'base'
		stacked_type = 'nested'
		arguments = [] #todo: add argument for template output file		


	@expose(help="This command generates a template", hide=True)
	def default(self):
		# todo: generate template
		template_filename = resource_filename('kinetic_datanator', 'data/InputTemplate.xlsx')
		output_filename = 'InputTemplate.xlsx'
		wb = openpyxl.load_workbook(filename)
		wb.save(output_filename)


class GetTaxonomicLineage(CementBaseController):
	class Meta:
		label = 'get-taxonomic-lineage'
		description = "This command prints out the taxonomic lineage of a species by number, corresponing to increasing nodes on the NCBI taxonomic tree"
		stacked_on = 'base'
		stacked_type = 'nested'
		arguments = [
		(['species_name'],dict(metavar='species_name', 
			type=str, help="name of the species whose lineage you are searching (example: 'homo sapiens')"))
		]


	@expose(help="This command prints out the taxonomic lineage of a species by number, corresponing to increasing nodes on the NCBI taxonomic tree", hide=True)
	def default(self):

		lineage = (taxon_finder.get_taxonomic_lineage(self.app.pargs.species_name))
		i = 1
		for node in lineage:
			print("{}: {}".format(i, node))
			i = i+1


class FindKineticsApp(CementApp):
	class Meta:
		label = "find-kinetics"
		base_controller = "base"
		handlers = [BaseController, GetKineticsController, GenerateTemplateController, GetTaxonomicLineage]


def main():
    with FindKineticsApp() as app:
        app.run()

if __name__ == "__main__":
    main()
