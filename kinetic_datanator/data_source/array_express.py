""" Downloads and parses the ArrayExpress database

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-08-16
:Copyright: 2017, Karr Lab
:License: MIT
"""


import datetime
import dateutil.parser
import demjson
import io
import jxmlease
import os
import pkg_resources
import requests.exceptions
import sqlalchemy
import sqlalchemy.ext.declarative
import sqlalchemy.orm
import sys
import warnings
import zipfile
from kinetic_datanator.core import data_source
from kinetic_datanator.data_source import download_ax


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

experiment_organism = sqlalchemy.Table(
	'experiment_organism', Base.metadata,
	sqlalchemy.Column('experiment__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('experiment._id'), index=True),
	sqlalchemy.Column('organism__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('organism_._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Experiment:Organism many-to-many association table

experiment_experiment_type= sqlalchemy.Table(
	'experiment_experiment_type', Base.metadata,
	sqlalchemy.Column('experiment__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('experiment._id'), index=True),
	sqlalchemy.Column('experiment_type_id', sqlalchemy.Integer, sqlalchemy.ForeignKey('experiment_type._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Experiment:ExperimentType many-to-many association table

experiment_experiment_design= sqlalchemy.Table(
	'experiment_experiment_design', Base.metadata,
	sqlalchemy.Column('experiment__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('experiment._id'), index=True),
	sqlalchemy.Column('experiment_design_id', sqlalchemy.Integer, sqlalchemy.ForeignKey('experiment_design._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Experiment:ExperimentDesign many-to-many association table

experiment_data_format = sqlalchemy.Table(
	'experiment_data_format', Base.metadata,
	sqlalchemy.Column('experiment__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('experiment._id'), index=True),
	sqlalchemy.Column('data_format_id', sqlalchemy.Integer, sqlalchemy.ForeignKey('data_format._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Experiment:DataFormat many-to-many association table


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

	sqlalchemy.schema.UniqueConstraint(name, value)

	__tablename__ = 'characteristic'


class Variable(Base):
	""" Represents an expiremental variable
	Attributes:
		samples
		name (:obj:`str`): name of the variable (e.g. genotype)
		value (:obj:`str`): value of variable (e.g control)
		unit (:obj:`str`): units of value (e.g control). This field is not always filled. 
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
		experiment_id: (:obj:`str`): the accesion number of the experiment
		index(:obj:`str`): name of the extract of the sample
		name (:obj:`str`): name of the source of the sample

		characteristics (:obj:`list` of :obj:`Characteristic'):
		variables (:obj:`list` of :obj:`Variable'):
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
	name = sqlalchemy.Column(sqlalchemy.String())
	name2 = sqlalchemy.Column(sqlalchemy.String())
	description = sqlalchemy.Column(sqlalchemy.String())
	organism_ = sqlalchemy.orm.relationship('Organism', secondary=experiment_organism, backref=sqlalchemy.orm.backref('experiments'))
	experiment_types = sqlalchemy.orm.relationship('ExperimentType', secondary=experiment_experiment_type, backref=sqlalchemy.orm.backref('experiments'))
	experiment_designs = sqlalchemy.orm.relationship('ExperimentDesign', secondary=experiment_experiment_design, backref=sqlalchemy.orm.backref('experiments'))
	submission_date = sqlalchemy.Column(sqlalchemy.String())
	release_date = sqlalchemy.Column(sqlalchemy.String())
	data_formats = sqlalchemy.orm.relationship('DataFormat', secondary=experiment_data_format, backref=sqlalchemy.orm.backref('experiments'))

	__tablename__ = 'experiment'


class ExperimentDesign(Base):
	_id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
	name = sqlalchemy.Column(sqlalchemy.String())
	sqlalchemy.schema.UniqueConstraint(name)

	__tablename__ = 'experiment_design'


class ExperimentType(Base):
	_id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
	name = sqlalchemy.Column(sqlalchemy.String())
	sqlalchemy.schema.UniqueConstraint(name)

	__tablename__ = 'experiment_type'


class DataFormat(Base):
	_id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
	experiment_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('experiment._id'))
	#experiment = sqlalchemy.orm.relationship('Experiment', backref=sqlalchemy.orm.backref(
	#	'data_format'), foreign_keys=[experiment_id])
	name = sqlalchemy.Column(sqlalchemy.String())

	__tablename__ = 'data_format'


class Organism(Base):
	_id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
	name = sqlalchemy.Column(sqlalchemy.String())
	sqlalchemy.schema.UniqueConstraint(name)

	__tablename__ = 'organism_'


class Extract(Base):
	_id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
	sample_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('sample._id'))
	sample = sqlalchemy.orm.relationship('Sample', backref=sqlalchemy.orm.backref(
		'extracts'), foreign_keys=[sample_id])
	name = sqlalchemy.Column(sqlalchemy.String())

	__tablename__ = 'sample_extract'


class Error(Base):
	_id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
	exp_samp_id = sqlalchemy.Column(sqlalchemy.String())
	error_message = sqlalchemy.Column(sqlalchemy.String())

	__tablename__ = 'error'


class ArrayExpress(data_source.HttpDataSource):
	""" A local sqlite copy of the ArrayExpress database

	Attributes:
		EXCLUDED_DATASET_IDS (:obj:`list` of :obj:`str`): list of IDs of datasets to exclude
	"""

	base_model = Base

	def __init__(self, name=None, cache_dirname=None, clear_content=False, load_content=False, max_entries=float('inf'),
                 commit_intermediate_results=False, download_backup=True, verbose=False,
                 clear_requests_cache=False, download_request_backup=False):
        """
        Args:
            name (:obj:`str`, optional): name
            cache_dirname (:obj:`str`, optional): directory to store the local copy of the data source and the HTTP requests cache
            clear_content (:obj:`bool`, optional): if :obj:`True`, clear the content of the sqlite local copy of the data source
            load_content (:obj:`bool`, optional): if :obj:`True`, load the content of the local sqlite database from the external source
            max_entries (:obj:`float`, optional): maximum number of entries to save locally
            commit_intermediate_results (:obj:`bool`, optional): if :obj:`True`, commit the changes throughout the loading
                process. This is particularly helpful for restarting this method when webservices go offline.
            download_backup (:obj:`bool`, optional): if :obj:`True`, load the local copy of the data source from the Karr Lab server
            verbose (:obj:`bool`, optional): if :obj:`True`, print status information to the standard output
            clear_requests_cache (:obj:`bool`, optional): if :obj:`True`, clear the HTTP requests cache
            download_request_backup (:obj:`bool`, optional): if :obj:`True`, download the request backup
        """
		super(ArrayExpress, self).__init__(name=name, cache_dirname=cache_dirname, clear_content=clear_content,
	                                       load_content=load_content, max_entries=max_entries,
	                                       commit_intermediate_results=commit_intermediate_results,
	                                       download_backup=download_backup, verbose=verbose,
	                                       clear_requests_cache=clear_requests_cache, download_request_backup=download_request_backup)

		with open(pkg_resources.resource_filename('kinetic_datanator', 'data_source/array_express_excluded_dataset_ids.txt'), 'r') as file:
			self.EXCLUDED_DATASET_IDS = [line.rstrip() for line in file]

	def load_experiments(self, json_object=None):
		db_session = self.session
		req_session = self.requests_session

		entry_details=json_object

		for single_entry in entry_details['experiments']['experiment']:
			experiment = Experiment()
			experiment.id = single_entry['accession']
			if type(single_entry['name']) is list:
				experiment.name  = single_entry['name'][0]
				experiment.name2 = single_entry['name'][1]
			else:		
				experiment.name = single_entry['name']

			if 'organism' in single_entry:
				for org in single_entry['organism']:
					experiment.organism_.append(self.return_relevant_object(Organism, {"name":org}))
			if 'description' in single_entry:
				experiment.description = single_entry['description'][0]['text']
			
			if 'experimenttype' in single_entry:
				entries = single_entry['experimenttype']
				for entry in entries:
					experiment.experiment_types.append(self.return_relevant_object(ExperimentType, {'name':entry}))

			if 'experimentdesign' in single_entry:
				entries = single_entry['experimentdesign']
				for entry in entries:
					experiment.experiment_designs.append(self.return_relevant_object(ExperimentDesign, {'name':entry}))

			if 'bioassaydatagroup' in single_entry:
				for entry in single_entry['bioassaydatagroup']:
					experiment.data_formats.append(self.return_relevant_object(DataFormat, {'name': entry['dataformat']}))

			if 'submissiondate' in single_entry:
				experiment.submission_date = single_entry['submissiondate']

			if 'releasedate' in single_entry:
				experiment.release_date = single_entry['releasedate']

			db_session.add(experiment)


	def load_content(self, experiment_ids=None, test_url=""):
		pass

	def load_experiments_from_text(self):
		"""
		Iterate through the dates and use load_content() to download everything to database in chunks
		"""
		db_session = self.session
		for filename in os.listdir('AllExperiments'):
			self.load_experiments(json_object=demjson.decode_file(os.path.join('AllExperiments', filename), encoding='utf-8'))
		db_session.commit()

	def load_single_sample(self, sample_jxmlease_object, ax_num):
		db_session = self.session
		sample = Sample()
		sample.experiment_id = ax_num
		sample.experiment = db_session.query(Experiment).filter_by(id = ax_num).first()
		if 'source' in sample_jxmlease_object:
			if ax_num not in ['E-MTAB-2067', 'E-MTAB-2617', 'E-MTAB-2325']:
				sample.name = sample_jxmlease_object['source']['name']
			else: #'E-MTAB-2067', 'E-MTAB-2617', 'E-MTAB-2325' has an error where the source is a list of two duplicates
				sample.name = sample_jxmlease_object['source'][0]['name']
		
		if 'extract' in sample_jxmlease_object:
			entries = sample_jxmlease_object['extract']
			if type(entries) is not list:
				entries = [entries]
			for entry in entries:
				sample.extracts.append(Extract(name=entry['name']))

		if 'assay' in sample_jxmlease_object:
			sample.assay = sample_jxmlease_object['assay']['name']

		# create a characteristic object for each characteristic and append that to the sample's characteristic field
		if 'characteristic' in sample_jxmlease_object:
			characteristics = sample_jxmlease_object['characteristic']
			for entry in characteristics:
				sample.characteristics.append(self.return_relevant_object(Characteristic, {'name': entry['category'], 'value':entry['value']}))


		# create a variable object for each variable and append that to the sample's variable field
		#print sample.experiment.id
		if 'variable' in sample_jxmlease_object:
			variables = sample_jxmlease_object['variable']
			for entry in variables:
				unit = None
				if "unit" in entry:
					unit = entry['unit']
				sample.variables.append(self.return_relevant_object(Variable, {'name': entry['name'], "value": entry['value'], 'unit':unit}))
		db_session.add(sample)

	def return_relevant_object(self, db_class, filters):
		"""
		For classes that have unique constraints. This method checks to see if an object
		already exists. If it does, then it return that object. If not, it creates a new object
		and return it. This method can only be used if there is an orm table explicitely created
		"""
		existing_object = self.session.query(db_class)
		for key, value in filters.items():
			existing_object = existing_object.filter(db_class.__dict__[key]==value)
		existing_object = existing_object.first()
		if existing_object is not None:
			return existing_object
		else:
			new_object = db_class()
			for key, value in filters.items():
				setattr(new_object, key, value) 
			return new_object

	def load_samples_from_text(self):
		"""
		Iterate through the dates and use load_content() to download everything to database in chunks
		"""

		db_session = self.session

		for filename in os.listdir('AllSamples'):
			if filename[:-4] not in self.EXCLUDED_DATASET_IDS:
				entry_details = demjson.decode_file(os.path.join('AllSamples', filename), encoding='utf-8')
				i = 1
				if 'sample' in entry_details['experiment']:
							if isinstance(entry_details['experiment']['sample'], list):
								for num, entry in enumerate(entry_details['experiment']['sample']):
									self.load_single_sample(entry_details['experiment']['sample'][num],filename[:-4])

							if isinstance(entry_details['experiment']['sample'], dict): #in case there is only one sample
								self.load_single_sample(entry_details['experiment']['sample'], filename[:-4])
				if i % 500 == 0:
					db_session.commit()
					print("COMMIT")
				i +=1
		db_session.commit()


def download_and_load_all_metadata():
	download_ax.download_all_metadata()
	src = array_express.ArrayExpress(download_backup=True)
	src.load_experiments_from_text()
	src.load_samples_from_text()
