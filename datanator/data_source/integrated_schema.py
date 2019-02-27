from datanator.core import data_source
from datanator.util import molecule_util
from xml import etree
import Bio.Alphabet
import Bio.SeqUtils
import bs4
import csv
import datetime
import html
import openbabel	
import libsbml
import math
import os
import pint
import pubchempy
import re
import requests
import requests_cache
import six
import sqlalchemy
import sqlalchemy.ext.declarative
import sqlalchemy.orm
import sys
import time
import warnings
import wc_utils.util.list
import wc_utils.workbook.core
import wc_utils.workbook.io
import pronto
import enum
import sqlite3
from sqlite3 import Error

Base = sqlalchemy.ext.declarative.declarative_base()


entry_synonym = sqlalchemy.Table(
    'entry_synonym', Base.metadata,
    sqlalchemy.Column('entry__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('entry._id'), index=True),
    sqlalchemy.Column('synonym__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('synonym._id'), index=True),
)

entry_resource = sqlalchemy.Table(
    'entry_resource', Base.metadata,
    sqlalchemy.Column('entry__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('entry._id'), index=True),
    sqlalchemy.Column('resource__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('resource._id'), index=True),
)

compound_compound_structure = sqlalchemy.Table(
	'compound_compound_structure', Base.metadata,
	sqlalchemy.Column('compound__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('compound_id'), index=True),
	sqlalchemy.Column('compound_structure_id', sqlalchemy.Integer, sqlalchemy.ForeignKey('compound_structure_id'), index=True),
)

class Synonym(Base):

	_id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
	label = sqlalchemy.Column(sqlalchemy.String(), index=True)
	__tablename__ = 'synonym'

class Entry(Base):

	_id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
	name = sqlalchemy.Column(sqlalchemy.String())
	synonyms = sqlalchemy.orm.relationship('Synonym', secondary=entry_synonym, backref=sqlalchemy.orm.backref('entries'))
	comments = sqlalchemy.Column(sqlalchemy.String())
	db_refs = sqlalchemy.orm.relationship('Resource', secondary=entry_resouce, backref=sqlalchemy.orm.backref('entries'))
	_type = sqlalchemy.Column(sqlalchemy.String())

	sqlalchemy.schema.UniqueConstraint(_id, _type)

	__tablename__ = 'entry'
	__mapper_args__ = {'polymorphic_on': _type}

class Compartment(Entry):

	_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('entry._id'), primary_key=True)
	name = sqlalchemy.Column(sqlalchemy.String())

	__tablename__ = 'compartment'
	__mapper_args__ = {'polymorphic_identity': 'compartment'}

class Compound(Entry):

	_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('entry._id'), primary_key=True)
	structures = sqlalchemy.orm.relationship('BioPolymerForm', secondary=compound_compound_structure,
                                             backref=sqlalchemy.orm.backref('compounds'))
	_emprical_formula = sqlalchemy.Column(sqlalchemy.String())
	_mol_wt = sqlalchemy.Column(sqlalchemy.Float())
	_charge = sqlalchemy.Column(sqlalchemy.Integer())

	__tablename__ = 'compound'
	__mapper_args__ = {'polymorphic_identity': 'compound'}

	def get_inchi_structures(self):
		'''
		Get INCHI-formatted structures
		'''
		return [s.value for s in self.structures if s.format == 'inchi']

	def get_smiles_structures(self):
        """ 
        Get SMILES-formatted structures
        """
		return [s.value for s in self.structures if s.format == 'smiles']

class CompoundStructure(Base):
	_id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
	value = sqlalchemy.Column(sqlalchemy.Text())
	_format = sqlalchemy.Column(sqlalchemy.String())
	_value_inchi = sqlalchemy.Column(sqlalchemy.Text(), index=True)
	_value_inchi_formula_connectivity = sqlalchemy.Column(sqlalchemy.Text(), index=True)
	sqlalchemy.schema.UniqueConstraint(value, _format)

	__tablename__ = 'compound_structure'

	def calc_inchi_formula_connectivity(self):
		if self.format == 'inchi':
			self._value_inchi = self.value
		else:
			try:
				self._value_inchi = molecule_util.Molecule(structure=self.value).to_inchi() or None
			except ValueError:
				self._value_inchi = None
		if self._value_inchi:
			self._value_inchi_formula_connectivity = molecule_util.InchiMolecule(self._value_inchi).get_formula_and_connectivity()	


class BioPolymerForm:

	def __init__(self, type, structure):
		self.type = 'Dna'
		self.structure = 'InChi=1sACGU[]CGACC'
	
	def dna_structure(self, type, structure):
		if self.type.lower() == 'dna':
			# annotated sequence
			return structure
	
	def rna_structure(self, type, structure):
		if self.type.lower() == 'rna':
			# annotated sequence
			return structure
	
	def protein(self, type, structure):
		if self.type.lower() == 'protein':
			return structure

class Metabolite(Compound):
	structure = openbabel.OBMol()

class RNA(Compound):
	structure = BioPolymerForm('Rna')

class DNA(Compound):
	structure = BioPolymerForm('Dna')

class Protein(Compound):
	structure = BioPolymerForm('Protein')


class InteractionParticipant(Entry):
	_id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
	compound_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('compound._id'), index=True)
	compound = sqlalchemy.orm.relationship('Compound', backref=sqlalchemy.orm.backref('interaction_participants'), foreign_keys=[compound_id])
	compartment_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('compartment._id'), index=True)
	compartment = sqlalchemy.orm.relationship('Compartment', backref=sqlalchemy.orm.backref('interaction_participants'), foreign_keys=[compartment_id])
	coefficient = sqlalchemy.Column(sqlalchemy.Float())
	type = sqlalchemy.Column(sqlalchemy.String())

	__tablename__ = 'interation_participant'
	__mapper_args__ = {'polymorphic_identity': 'interaction_participant'}

class Complex(Compound):
	#
	structure = None

class Enzyme(Entry):
	"""
	represents enzymes 
	"""
	_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('entry._id'), primary_key=True)
	molecular_weight = sqlalchemy.Column(sqlalchemy.Float())
	__tablename__ = 'enzyme'
	__mapper_args__ = {'polymorphic_identity': 'enzyme'}

class Reaction(Entry):
	participants = sqlalchemy.orm.relationship('InteractionParticipant', backref=sqlalchemy.orm.backref('interaction_participants'),
                                            foreign_keys=[InteractionParticipant.compound_id],
                                            cascade='all, delete-orphan')
	#enzymes: things that affect the rate law but don't undergo permanent covalent transformations
	modifiers = sqlalchemy.orm.relationship('InteractionParticipant', backref=sqlalchemy.orm.backref('modifier_kinetic_law'),
                                            foreign_keys=[InteractionParticipant.compound_id],
                                            cascade='all, delete-orphan')
	__tablename__ = 'reaction'
	__mapper_args__ = {'polymorphic_identity': 'reaction'}

wcm_ontology = pronto.Ontology('./WCM.obo') #biontology API to switch to http link
'''
[Term]
id = quantitative_property
descrition = ''

[Term]
id = descriptive_property
descrition = ''
'''

class Observation(Base):
	participants = sqlalchemy.orm.relationship('InteractionParticipant', backref=sqlalchemy.orm.backref('interaction_participants'),
                                            foreign_keys=[InteractionParticipant.compound_id],
                                            cascade='all, delete-orphan')
	_property = Ontology.Term(wcm_notoloy)
	value = sqlalchemy.orm.relationship('BioPolymerForm', secondary=compound_compound_structure,
                                             backref=sqlalchemy.orm.backref('compounds'))
	# biological system observed 
	taxon = sqlalchemy.Column(sqlalchemy.Integer(), nullable = True)
	wildtype = sqlalchemy.Column(sqlalchemy.Boolean(), nullable = True)
	genetic_variant = sqlalchemy.Column(sqlalchemy.String(), nullable = True)
	tissue = sqlalchemy.Column(sqlalchemy.String())
	physiological_state = sqlalchemy.Column(sqlalchemy.String())

	# environmental conditions in which the measurement was made
	temp = sqlalchemy.Column(sqlalchemy.Float(), nullable = True)
	temp_units = sqlalchemy.Column(sqlalchemy.String(), nullable = True)
	ph = sqlalchemy.Column(sqlalchemy.Float(), nullable = True)	
	ph_units = sqlalchemy.Column(sqlalchemy.String(), nullable = True)
	growth_media = sqlalchemy.Column(sqlalchemy.String(), nullable = True)
	conditions = sqlalchemy.Column(sqlalchemy.String(), nullable = True)

	# how observed value was measured and reduced
	experiment_type = sqlalchemy.Column(sqlalchemy.String(), nullable = True) # Rna-seq ChIP-Seq
	experiment_design = sqlalchemy.Column(sqlalchemy.String(), nullable = True) # additional info
	experiment_method = sqlalchemy.Column(sqlalchemy.String(), nullable = True)
	analysis_method = sqlalchemy.Column(sqlalchemy.String(), nullable = True)

	# comments
	comments = sqlalchemy.Column(sqlalchemy.UnicodeText(), nullable = True)

	# provenance 
	source = sqlalchemy.Column(sqlalchemy.String()) # sabio, ecmdb, arrayexpress
	source_version = sqlalchemy.Column(sqlalchemy.String()) # source version of the database
	source_id = sqlalchemy.Column(sqlalchemy.String()) # identifier from the source
	source_url = sqlalchemy.Column(sqlalchemy.String())
	source_create = sqlalchemy.Column(sqlalchemy.DateTime,default=datetime.datetime.utcnow())
	source_modified = sqlalchemy.Column(sqlalchemy.DateTime,onupdate=datetime.datetime.utcnow())

	# time when data was retrieved from the source
	retrieval_method = sqlalchemy.Column(sqlalchemy.String()) # scripts used in data_source to retrieve the data
	retrieval_date = sqlalchemy.Column(sqlalchemy.DateTime,default=datetime.datetime.utcnow())

	pubs = ...

	__tablename__='observation'
	__mapper_args__ = {'polymorphic_on': type}


class QuantitativeObservation(Observation):
	value = sqlalchemy.Column(Float())
	error = sqlalchemy.Column(Float())
	units = sqlalchemy.Column(String())



class QualitativeObservation(Observation):
	value = sqlalchemy.Column(String())

class Resource(Base):
	namespace = sqlalchemy.Column(String()) # pubmed, pubchem, kegg.compound
	_id = sqlalchemy.Column(String(), primary_key=True) # unique identifier within namespace
	sqlalchemy.schema.UniqueConstraint(namespace, _id)
	__tablename__ = 'resource'

external_reference_publication = sqlalchemy.Table(
    'external_reference_publication', Base.metadata,
    sqlalchemy.Column('external_reference__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('external_references._id'), index=True),
    sqlalchemy.Column('publication__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('publications._id'), index=True),
)

class Publication(Base):
	type = sqlalchemy.Column(sqlalchemy.Enum(self, article, book, thesis, website)) # ???
	_id = sqlalchemy.Column(String(), primary_key=True)
	title = sqlalchemy.Column(sqlalchemy.String())
	author = sqlalchemy.Column(sqlalchemy.String())
	editor = sqlalchemy.Column(sqlalchemy.String())
	year = sqlalchemy.Column(sqlalchemy.Integer())
	publication = sqlalchemy.Column(sqlalchemy.String())
	publisher = sqlalchemy.Column(sqlalchemy.String())
	volume = sqlalchemy.Column(sqlalchemy.Integer())
	number = sqlalchemy.Column(sqlalchemy.Integer())
	issue = sqlalchemy.Column(sqlalchemy.Integer())
	edition = sqlalchemy.Column(sqlalchemy.Integer())
	chapter = sqlalchemy.Column(sqlalchemy.Integer())
	pages = sqlalchemy.Column(sqlalchemy.String())
	db_refs = sqlalchemy.orm.relationship("ExternalReference", secondary = external_references_publication, backref='publications')
	__tablename__ = 'publications'


class SabirrayEC(data_source.HttpDataSource):

	def __init__(self, name=None, cache_dirname=None, clear_content=False, load_content=False, max_entries=float('inf'),
		commit_intermediate_results=False, download_backups=True, verbose=False,
		clear_requests_cache=False, download_request_backup=False,
		webservice_batch_size=1, excel_batch_size=100,quilt_owner=None, quilt_package=None):

	    """
        Args:
            name (:obj:`str`, optional): name
            cache_dirname (:obj:`str`, optional): directory to store the local copy of the data source and the HTTP requests cache
            clear_content (:obj:`bool`, optional): if :obj:`True`, clear the content of the sqlite local copy of the data source
            load_content (:obj:`bool`, optional): if :obj:`True`, load the content of the local sqlite database from the external source
            max_entries (:obj:`float`, optional): maximum number of entries to save locally
            commit_intermediate_results (:obj:`bool`, optional): if :obj:`True`, commit the changes throughout the loading
                process. This is particularly helpful for restarting this method when webservices go offline.
            download_backups (:obj:`bool`, optional): if :obj:`True`, load the local copy of the data source from the Karr Lab server
            verbose (:obj:`bool`, optional): if :obj:`True`, print status information to the standard output
            clear_requests_cache (:obj:`bool`, optional): if :obj:`True`, clear the HTTP requests cache
            download_request_backup (:obj:`bool`, optional): if :obj:`True`, download the request backup
            webservice_batch_size (:obj:`int`, optional): default size of batches to download kinetic information from the SABIO webservice
            excel_batch_size (:obj:`int`, optional): default size of batches to download kinetic information from the SABIO
                Excel download service
            quilt_owner (:obj:`str`, optional): owner of Quilt package to save data
            quilt_package (:obj:`str`, optional): identifier of Quilt package to save data
        """

		self.webservice_batch_size = webservice_batch_size
		self.excel_batch_size = excel_batch_size

		super(SabioRk, self).__init__(name=name, cache_dirname=cache_dirname, clear_content=clear_content,
			load_content=load_content, max_entries=max_entries,
			commit_intermediate_results=commit_intermediate_results,
			download_backups=download_backups, verbose=verbose,
			clear_requests_cache=clear_requests_cache, download_request_backup=download_request_backup,
			quilt_owner=quilt_owner, quilt_package=quilt_package)

	#load three sqlite files (ArrayExpress, SabioRK, ECMDB)
# 	def load_content(self):

	def create_connection(db_file):

		""" create a database connection to the SQLite database
		    specified by the db_file
		:param db_file: database file
		:return: Connection object or None
		"""

		try:
			conn = sqlite3.connect(db_file)
			return conn
		except Error as e:
			print(e)
		
		return None

	def load_content(self):
		# directory where the three databases are
		array_express = './cache/ArrayExpress.sqlite'
		ecmdb = './cache/Ecmdb.sqlite'
		sabio_rk = './cache/SabioRk.sqlite'

		conn_array = create_connection(array_express)
		conn_ecmdb = create_connection(ecmdb)
		conn_sabio = create_connection(sabio_rk)



def main():


if __name__ = '__main__':
	main()