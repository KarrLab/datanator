from cement.core.foundation import CementApp
from cement.core.controller import CementBaseController, expose
import os.path
import openpyxl
from KineticDatanator import Datanator
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
			metavar='FLOAT', help = "minimum temperature", default=30)),
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



		#print ECNumberFinder.formatECForSabio("a")
	#__main__.main(self.app.pargs.input_data_file, self.app.pargs.output_data_file, self.app.pargs.species)
		

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

#tempRange = [30, 40], enzymeType = "wildtype", phRange = [5,9], proximLimit=1000):

		enzymeType = "wildtype"
		if self.app.pargs.include_mutants:
			enzymeType = "(wildtype OR mutant)"
		Datanator.getKineticData(self.app.pargs.input_data_file, self.app.pargs.output_data_file,
		 self.app.pargs.species, tempRange=[self.app.pargs.min_temp, self.app.pargs.max_temp], phRange = [self.app.pargs.min_ph, self.app.pargs.max_ph],
		 enzymeType = enzymeType, proximLimit = self.app.pargs.proxim_limit)



class GenerateTemplateController(CementBaseController):
	class Meta:
		label = 'generate-template'
		description = "This command generates a template"
		stacked_on = 'base'
		stacked_type = 'nested'
		arguments = []
		#arguments = [(['input_data_file'],dict(metavar='input-data-file', 
		#	type=str, help="path to the input data spreadsheet (xlsx)"))]

		#path = "/home/yosef/Desktop/messingWithCement/blue"
		#filename = "TemplateForDatanator2.xlsx"
		#if not os.path.exists(path):
		#	os.makedirs(path)
		#if not os.path.isfile('SmilesStuff23.xlsx'):
		#wb = openpyxl.load_workbook(filename='./KineticDatanator/TemplateDocument.xlsx') 
		#wb.save(".TheTemplateDocument.xlsx")#, as_template=False)


	@expose(help="This command generates a template", hide=True)
	def default(self):
		self.app.log.info("In generate_template")
	#	self.app.log.info("Input file: '{}'".format(self.app.pargs.input_data_file))
		print(os.path.dirname(os.path.realpath('inchiGenerator.py')))
		print os.path.realpath('inchiGenerator.py')
		print os.getcwd()
		wb = openpyxl.load_workbook(os.path.join(os.getcwd(), 'TemplateDocument.xlsx'))
		#wb = openpyxl.load_workbook(filename = 'TemplateDocument.xlsx')
		wb.save(".TheTemplateDocument.xlsx")#, as_template=False)

class FindKineticsApp(CementApp):
	class Meta:
		label = "find-kinetics"
		base_controller = "base"
		handlers = [BaseController, GetKineticsController, GenerateTemplateController]


with FindKineticsApp() as app:
	app.run()
