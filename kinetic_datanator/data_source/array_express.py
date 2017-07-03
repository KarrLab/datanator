import datetime
import dateutil.parser
import io
import json
import jxmlease
import requests.exceptions
import sqlalchemy
import sqlalchemy.ext.declarative
import sqlalchemy.orm
import warnings
import zipfile
from kinetic_datanator.core import data_source



Base = sqlalchemy.ext.declarative.declarative_base()
# :obj:`Base`: base model for local sqlite database

sample_characteristic = sqlalchemy.Table(
	'sample_characteristic', Base.metadata,
	sqlalchemy.Column('sample__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('sample._id'), index=True),
	sqlalchemy.Column('characteristic_id', sqlalchemy.Integer, sqlalchemy.ForeignKey('characteristic._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Sample:Characteristic many-to-many association table

sample_variable = sqlalchemy.Table(
	'sample_variable', Base.metadata,
	sqlalchemy.Column('sample__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('sample._id'), index=True),
	sqlalchemy.Column('variable_id', sqlalchemy.Integer, sqlalchemy.ForeignKey('variable._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Sample:Variable many-to-many association table


class Characteristic(Base):
	""" Represents an expiremental characteristic
	Attributes:
		samples
		name (:obj:`str`): name of the characteristic (e.g. organism)
		value (:obj:`str`): value of characteristic (e.g Mus musculus)
	"""
	_id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
	name = sqlalchemy.Column(sqlalchemy.String())
	value = sqlalchemy.Column(sqlalchemy.String())

	#sqlalchemy.schema.UniqueConstraint(name, value)

	__tablename__ = 'characteristic'


class Variable(Base):
	""" Represents an expiremental variable
	Attributes:
		samples
		name (:obj:`str`): name of the variable (e.g. genotype)
		value (:obj:`str`): value of variable (e.g control)
	"""
	_id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
	name = sqlalchemy.Column(sqlalchemy.String())
	value = sqlalchemy.Column(sqlalchemy.String())
	unit = sqlalchemy.Column(sqlalchemy.String())

	#sqlalchemy.schema.UniqueConstraint(name, value)

	__tablename__ = 'variable'


class Sample(Base):
	""" Represents an observed concentration
	Attributes:
		experiment
		index: name of the extract of the sample
		name (:obj:`str`): name of the source of the sample
		characteristics
		variables
	"""
	_id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
	experiment_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('experiment._id'), index=True)
	experiment = sqlalchemy.orm.relationship('Experiment', backref=sqlalchemy.orm.backref('samples'), foreign_keys=[experiment_id])
	#index = sqlalchemy.Column(sqlalchemy.Integer())
	name = sqlalchemy.Column(sqlalchemy.String())
	assay = sqlalchemy.Column(sqlalchemy.String())
	extract = sqlalchemy.Column(sqlalchemy.String())
	characteristics = sqlalchemy.orm.relationship('Characteristic',
												  secondary=sample_characteristic, backref=sqlalchemy.orm.backref('samples'))
	variables = sqlalchemy.orm.relationship('Variable',
											secondary=sample_variable, backref=sqlalchemy.orm.backref('samples'))

	#sqlalchemy.schema.UniqueConstraint(experiment_id, index)
	#sqlalchemy.schema.UniqueConstraint(experiment_id, name)
	#sqlalchemy.schema.UniqueConstraint(experiment_id, extract, index, name)

	__tablename__ = 'sample'


class Experiment(Base):
	"""
	Attributes:
		id
		samples
	"""
	_id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
	id = sqlalchemy.Column(sqlalchemy.String(), unique=True)

	organism = sqlalchemy.Column(sqlalchemy.String())
	name = sqlalchemy.Column(sqlalchemy.String())
	description = sqlalchemy.Column(sqlalchemy.String())
	experiment_type = sqlalchemy.Column(sqlalchemy.String())
	# samples = sqlalchemy.orm.relationship('Sample',
	#                                             secondary=experiment_sample, backref=sqlalchemy.orm.backref('experiment'))

	__tablename__ = 'experiment'


class ExperimentDesign(Base):
	_id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
	experiment_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('experiment._id'))
	experiment = sqlalchemy.orm.relationship('Experiment', backref=sqlalchemy.orm.backref(
		'experiment_designs'), foreign_keys=[experiment_id])
	name = sqlalchemy.Column(sqlalchemy.String())

	__tablename__ = 'experiment_design'

class ExperimentType(Base):
	_id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
	experiment_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('experiment._id'))
	experiment = sqlalchemy.orm.relationship('Experiment', backref=sqlalchemy.orm.backref(
		'experiment_types'), foreign_keys=[experiment_id])
	name = sqlalchemy.Column(sqlalchemy.String())

	__tablename__ = 'experiment_type'

class Error(Base):
	_id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
	exp_samp_id = sqlalchemy.Column(sqlalchemy.String())
	error_message = sqlalchemy.Column(sqlalchemy.String())

	__tablename__ = 'error'




class ArrayExpress(data_source.HttpDataSource):
	""" A local sqlite copy of the ECMDB database
	Attributes:
		DOWNLOAD_INDEX_URL (:obj:`str`): URL to download an index of array express experiments
		DOWNLOAD_COMPOUND_URL (:obj:`str`): URL pattern to download samples of a single experiment
	"""

	base_model = Base
	ENDPOINT_DOMAIN = 'https://www.ebi.ac.uk/arrayexpress/json/v3/experiments'
	DOWNLOAD_SAMPLE_URL = ENDPOINT_DOMAIN + '/{}/samples'
	DOWNLOAD_COMPLETE_SAMPLE_URL = 'https://www.ebi.ac.uk/arrayexpress/xml/v3/experiments/samples'
	SAMPLE_EXCLUDE = ['E-GEOD-7812', 'E-GEOD-7264']

	def load_samples(self, experiments):
		""" """

		def load_single_sample(sample_jxmlease_object):
			sample = Sample()
			sample.experiment = experiment
			print("Sample: {}".format(experiment.id))
			#sample.index = self.get_node_text(entry_details['experiment']['sample'][num]['extract']['name'])
			print sample_jxmlease_object['source']['name']
			if 'source' in sample_jxmlease_object:
				sample.name = sample_jxmlease_object['source']['name']
			if 'extract' in self.get_node_text(sample_jxmlease_object):
				sample.extract = self.get_node_text(sample_jxmlease_object['extract']['name'])
			if 'assay' in self.get_node_text(sample_jxmlease_object):
				sample.assay = self.get_node_text(sample_jxmlease_object['assay']['name'])

			# create a characteristic object for each characteristic and append that to the sample's characteristic field
			if 'characteristic' in sample_jxmlease_object:
				characteristics = self.get_node_text(sample_jxmlease_object['characteristic'])
				for entry in characteristics:
					new_charachteristic = Characteristic()
					new_charachteristic.name = entry['category']
					new_charachteristic.value = entry['value']
					sample.characteristics.append(new_charachteristic)


			# create a variable object for each variable and append that to the sample's variable field
			#print sample.experiment.id
			if 'variable' in sample_jxmlease_object:
				variables = sample_jxmlease_object['variable']
				for entry in variables:
					new_variable = Variable()
					new_variable.name = entry['name']
					new_variable.value = entry['value']
					if "unit" in entry:
						new_variable.unit = entry['unit']
					sample.variables.append(new_variable)
			db_session.add(sample)


		req_session = self.requests_session
		db_session = self.session

		#xml_parser = jxmlease.Parser()
		for experiment in experiments:
			try:
				if experiment.id not in self.SAMPLE_EXCLUDE:
					response = req_session.get(self.DOWNLOAD_SAMPLE_URL.format(experiment.id))
					response.raise_for_status()
					entry_details = json.loads(response.text)#xml_parser(response.text)

					if isinstance(entry_details['experiment']['sample'], list):
						for num, entry in enumerate(entry_details['experiment']['sample']):
							load_single_sample(entry_details['experiment']['sample'][num])

					if isinstance(entry_details['experiment']['sample'], dict): #in case there is only one sample
						load_single_sample(entry_details['experiment']['sample'])
			
			except (requests.exceptions.HTTPError, KeyError, TypeError) as e:# TypeError or KeyError or requests.exceptions.HTTPError as e:
				#print e
				error = Error()
				error.exp_samp_id = experiment.id
				error.error_message = "{}".format(e)
				db_session.add(error)
			




	def load_experiments(self, experiment_ids=None, test_url=""):
		db_session = self.session
		req_session = self.requests_session

		if test_url:
			url = test_url
		elif experiment_ids is None:
			url = self.ENDPOINT_DOMAIN
		else:
			url = self.ENDPOINT_DOMAIN + '/' + ','.join([str(id) for id in experiment_ids])
		response = req_session.get(url)
		response.raise_for_status()
		xml_parser = jxmlease.Parser()


		entry_details = json.loads(response.text)#xml_parser(response.text)
		#print xml_parser(response.text)
		#print entry_details


		for single_entry in entry_details['experiments']['experiment']:
			try:
				experiment = Experiment()
				#experiment.id = self.get_node_text(single_entry['accession'])
				experiment.id = single_entry['accession']
				#experiment.name = self.get_node_text(single_entry['name'])
				experiment.name = single_entry['name']

				#experiment.experiment_type = self.get_node_text(single_entry['experimenttype'])
				if 'organism' in single_entry:
					if len(single_entry['organism'])>1:
						pass
					else:
						experiment.organism = single_entry['organism'][0]
				if 'description' in single_entry:
					experiment.description = single_entry['description'][0]['text']
				
				if 'experimenttype' in single_entry:
					#if isinstance(single_entry['experimenttype'], list):
					entries = single_entry['experimenttype']

					#else:
					#	entries = [single_entry['experimenttype']]
					for entry in entries:
						experiment.experiment_types.append(ExperimentType(name=entry))

				if 'experimentdesign' in single_entry:
					#if isinstance(single_entry['experimentdesign'], list):
					entries = single_entry['experimentdesign']
					
					#else:
					#	entries = [single_entry['experimentdesign']]
					for entry in entries:
						experiment.experiment_designs.append(ExperimentDesign(name=entry))
					

				#print experiment.__dict__

				db_session.add(experiment)
			
			except TypeError or KeyError or requests.exceptions.HTTPError, e:
				#print e
				error = Error()
				error.exp_samp_id = experiment.id
				error.error_message = "{}".format(e)
				db_session.add(error)
			


	def load_content(self, experiment_ids=None, test_url=""):
		db_session = self.session


		# retrieve list of all experiments
		self.load_experiments(experiment_ids, test_url)
		#db_session.commit()


		# retrieve the samples for all of the experiments
		experiments = db_session.query(Experiment).all()
		self.load_samples(experiments)

		# save the changes to the database file
		db_session.commit()


	def load_content_in_chunks(self):
		"""
		Iterate through the dates and use load_content() to download everything to database in chunks
		"""


		self.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments?date=[2001-06-31+2001-12-30]")

		i = 2002
		current_year = datetime.datetime.now().year
		while i <= datetime.datetime.now().year:
			winter_spring = "[{}-01-01+{}-06-31]".format(i,i)
			summer_fall = "[{}-07-01+{}-12-30]".format(i,i)
			self.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments?date={}".format(winter_spring))
			self.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/json/v3/experiments?date={}".format(summer_fall))
			i = i + 1





	def get_node_text(self, node):
		""" Get the next of a XML node
		Args:
			node (:obj:`jxmlease.cdatanode.XMLCDATANode` or :obj:`str`): XML node or its text
		Returns:
			:obj:`str`: text of the node
		"""
		if isinstance(node, jxmlease.cdatanode.XMLCDATANode):
			return node.get_cdata()
		return node


if __name__ == '__main__':
	a = ArrayExpress()
	a.load_content_in_chunks()
	#a.load_content(test_url="https://www.ebi.ac.uk/arrayexpress/xml/v3/experiments/E-GEOD-7812,E-MTAB-1234")
	#becuase of E-GEOD-7812, add in http requests exceptions







