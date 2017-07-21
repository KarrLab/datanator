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
import os
import sys
import download_ax
reload(sys)  
sys.setdefaultencoding('utf8')


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
	""" A local sqlite copy of the ECMDB database
	Attributes:
		DOWNLOAD_INDEX_URL (:obj:`str`): URL to download an index of array express experiments
		DOWNLOAD_COMPOUND_URL (:obj:`str`): URL pattern to download samples of a single experiment
	"""

	base_model = Base

	EXCLUSIONS = [u'E-MEXP-21', u'E-GEOD-714', u'E-MANP-2', u'E-NASC-6', u'E-NASC-10', u'E-GEOD-678', u'E-GEOD-652', u'E-GEOD-537', u'E-MAXD-1', u'E-GEOD-2109', u'E-MAXD-10', u'E-TABM-114', u'E-GEOD-7264', u'E-GEOD-7260', u'E-GEOD-8737', u'E-GEOD-8653', u'E-TABM-76', u'E-GEOD-8251', u'E-GEOD-7812', u'E-TABM-1', u'E-GEOD-10645', u'E-GEOD-9376', u'E-GEOD-8858', u'E-GEOD-7740', u'E-MTAB-27', u'E-GEOD-13317', u'E-GEOD-14559', u'E-MTAB-28', u'E-GEOD-16752', u'E-GEOD-16393', u'E-GEOD-15907', u'E-TABM-3', u'E-GEOD-11202', u'E-GEOD-15292', u'E-GEOD-26284', u'E-GEOD-21264', u'E-GEOD-25906', u'E-GEOD-25249', u'E-GEOD-25248', u'E-GEOD-21925', u'E-TABM-930', u'E-GEOD-18906', u'E-GEOD-25045', u'E-GEOD-25025', u'E-GEOD-22162', u'E-GEOD-24565', u'E-GEOD-20624', u'E-GEOD-23853', u'E-GEOD-23537', u'E-GEOD-19491', u'E-MTAB-62', u'E-GEOD-22778', u'E-TABM-887', u'E-TABM-475', u'E-MEXP-2637', u'E-GEOD-2415', u'E-GEOD-21478', u'E-GEOD-2295', u'E-GEOD-1981', u'E-GEOD-1714', u'E-GEOD-1366', u'E-GEOD-1361', u'E-GEOD-1356', u'E-GEOD-1342', u'E-GEOD-1290', u'E-GEOD-14774', u'E-GEOD-13870', u'E-GEOD-13056', u'E-GEOD-12821', u'E-GEOD-12727', u'E-GEOD-12717', u'E-GEOD-12354', u'E-GEOD-12322', u'E-GEOD-12263', u'E-GEOD-12249', u'E-GEOD-12238', u'E-GEOD-12157', u'E-GEOD-11500', u'E-GEOD-11444', u'E-GEOD-11255', u'E-GEOD-11248', u'E-GEOD-11247', u'E-GEOD-17762', u'E-GEOD-18700', u'E-MEXP-1995', u'E-MTAB-686', u'E-MTAB-622', u'E-MTAB-570', u'E-MTAB-510', u'E-MTAB-291', u'E-GEOD-27130', u'E-GEOD-26095', u'E-GEOD-34448', u'E-GEOD-26494', u'E-GEOD-33518', u'E-GEOD-32648', u'E-GEOD-31381', u'E-GEOD-32658', u'E-GEOD-31210', u'E-GEOD-33321', u'E-GEOD-25219', u'E-GEOD-33213', u'E-GEOD-25192', u'E-GEOD-32970', u'E-GEOD-32517', u'E-GEOD-32465', u'E-GEOD-27480', u'E-MTAB-365', u'E-GEOD-32019', u'E-GEOD-31946', u'E-GEOD-26331', u'E-GEOD-31908', u'E-GEOD-31548', u'E-GEOD-31863', u'E-GEOD-27317', u'E-GEOD-31477', u'E-GEOD-24837', u'E-GEOD-24693', u'E-GEOD-31038', u'E-GEOD-31039', u'E-GEOD-30725', u'E-GEOD-30567', u'E-GEOD-19301', u'E-GEOD-29619', u'E-GEOD-28919', u'E-GEOD-27491', u'E-GEOD-28791', u'E-GEOD-30263', u'E-GEOD-23143', u'E-GEOD-26338', u'E-GEOD-27923', u'E-GEOD-29174', u'E-GEOD-28884', u'E-GEOD-29692', u'E-GEOD-29611', u'E-TABM-1140', u'E-GEOD-26367', u'E-GEOD-26022', u'E-GEOD-25869', u'E-GEOD-28631', u'E-GEOD-24836', u'E-GEOD-25066', u'E-GEOD-25055', u'E-GEOD-28746', u'E-GEOD-24710', u'E-GEOD-26500', u'E-GEOD-20964', u'E-GEOD-20140', u'E-GEOD-39677', u'E-GEOD-34200', u'E-GEOD-33072', u'E-GEOD-37074', u'E-GEOD-41258', u'E-GEOD-40869', u'E-GEOD-35583', u'E-GEOD-35734', u'E-GEOD-36030', u'E-GEOD-35239', u'E-GEOD-31437', u'E-GEOD-48017', u'E-GEOD-47990', u'E-GEOD-47951', u'E-GEOD-47845', u'E-GEOD-47542', u'E-GEOD-34665', u'E-GEOD-53261', u'E-GEOD-47983', u'E-ERAD-186', u'E-GEOD-49905', u'E-GEOD-45463', u'E-GEOD-44874', u'E-GEOD-49530', u'E-GEOD-49527', u'E-GEOD-36369', u'E-GEOD-39672', u'E-GEOD-48417', u'E-GEOD-48415', u'E-GEOD-48405', u'E-GEOD-48377', u'E-GEOD-48310', u'E-GEOD-48279', u'E-GEOD-40967', u'E-GEOD-46712', u'E-GEOD-46517', u'E-GEOD-45892', u'E-GEOD-44777', u'E-GEOD-36889', u'E-GEOD-39518', u'E-GEOD-36796', u'E-GEOD-45468', u'E-GEOD-45159', u'E-GEOD-45149', u'E-GEOD-44944', u'E-GEOD-45480', u'E-GEOD-14217', u'E-ERAD-321', u'E-ERAD-319', u'E-GEOD-63341', u'E-GEOD-56047', u'E-GEOD-56045', u'E-GEOD-49417', u'E-GEOD-51341', u'E-GEOD-51338', u'E-GEOD-62992', u'E-GEOD-62564', u'E-GEOD-44722', u'E-ERAD-305', u'E-GEOD-54470', u'E-GEOD-60863', u'E-GEOD-60341', u'E-ERAD-287', u'E-GEOD-59923', u'E-GEOD-59913', u'E-GEOD-59905', u'E-GEOD-55347', u'E-GEOD-53348', u'E-GEOD-41119', u'E-GEOD-59150', u'E-GEOD-57611', u'E-GEOD-59097', u'E-GEOD-53643', u'E-GEOD-53080', u'E-GEOD-57530', u'E-GEOD-57822', u'E-GEOD-57815', u'E-GEOD-57542', u'E-MTAB-2067', u'E-GEOD-53165', u'E-GEOD-75685', u'E-GEOD-75268', u'E-GEOD-73103', u'E-GEOD-50410', u'E-MTAB-3732', u'E-ERAD-374', u'E-GEOD-39332', u'E-ERAD-412', u'E-GEOD-64844', u'E-GEOD-65858', u'E-GEOD-69597', u'E-GEOD-69498', u'E-GEOD-63120', u'E-GEOD-57739', u'E-GEOD-52903', u'E-GEOD-69180', u'E-MTAB-2919', u'E-GEOD-69004', u'E-GEOD-68984', u'E-GEOD-68972', u'E-MTAB-2617', u'E-GEOD-62625', u'E-GEOD-60836', u'E-MTAB-2325', u'E-GEOD-62372', u'E-GEOD-56749', u'E-GEOD-64763', u'E-GEOD-45218', u'E-GEOD-63429', u'E-GEOD-63246', u'E-GEOD-63042', u'E-GEOD-62734', u'E-GEOD-62292', u'E-GEOD-61635', u'E-GEOD-61628', u'E-GEOD-61626', u'E-GEOD-61582', u'E-MTAB-5214', u'E-GEOD-73518', u'E-GEOD-73515', u'E-GEOD-84422', u'E-MTAB-4888', u'E-GEOD-83951', u'E-GEOD-70936', u'E-GEOD-82549', u'E-GEOD-82545', u'E-GEOD-82543', u'E-GEOD-82539', u'E-GEOD-82537', u'E-GEOD-82534', u'E-GEOD-82532', u'E-GEOD-83160', u'E-GEOD-75330', u'E-GEOD-60690', u'E-GEOD-73290', u'E-GEOD-63467', u'E-MTAB-4032', u'E-GEOD-65391', u'E-GEOD-75220', u'E-GEOD-69979', u'E-GEOD-70185', u'E-GEOD-62044', u'E-GEOD-70774', u'E-GEOD-69872', u'E-GEOD-76809', u'E-MTAB-3947', u'E-GEOD-71585', u'E-MTAB-5522', u'E-MTAB-4547', u'E-GEOD-57362', u'E-MTAB-4388']

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
			self.load_experiments(json_object=json.loads(open("AllExperiments/{}".format(filename), 'r').read().encode('utf8')))
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
		for key, value in filters.iteritems():
			existing_object = existing_object.filter(db_class.__dict__[key]==value)
		existing_object = existing_object.first()
		if existing_object is not None:
			return existing_object
		else:
			new_object = db_class()
			for key, value in filters.iteritems():
				setattr(new_object, key, value) 
			return new_object

		


	def load_samples_from_text(self):
		"""
		Iterate through the dates and use load_content() to download everything to database in chunks
		"""

		db_session = self.session

		for filename in os.listdir('AllSamples'):
			if filename[:-4] not in self.EXCLUSIONS:
				entry_details =  json.loads(open("AllSamples/{}".format(filename), 'r').read())
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
