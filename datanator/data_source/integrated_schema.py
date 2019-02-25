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

class BioPolymerForm(Base):
	_id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
	value = sqlalchemy.Column(sqlalchemy.Text())
	_format = sqlalchemy.Column(sqlalchemy.String())
	_value_inchi = sqlalchemy.Column(sqlalchemy.Text(), index=True)
	_value_inchi_formula_connectivity = sqlalchemy.Column(sqlalchemy.Text(), index=True)
	sqlalchemy.schema.UniqueConstraint(value, format)

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

class Entry(Base):
	_id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
	name = sqlalchemy.Column(sqlalchemy.String())
	synonyms = sqlalchemy.orm.relationship('Synonym', secondary=entry_synonym, backref=sqlalchemy.orm.backref('entries'))
	comments = sqlalchemy.Column(sqlalchemy.String())
	db_refs = sqlalchemy.orm.relationship('Resource', secondary=external_reference, backref=sqlalchemy.orm.backref('external_references'))

	__tablename__ = 'entry'
	__mapper_args__ = {'polymorphic_on': _type}

class Compartment(Entry):

	_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('entry._id'), primary_key=True)
	name = sqlalchemy.Column(sqlalchemy.String())

	__tablename__ = 'compartment'
	__mapper_args__ = {'polymorphic_identity': 'compartment'}


class Compound(Entry):
	#?????
	structures = sqlalchemy.orm.relationship('BioPolymerForm', secondary=compound_compound_structure,
                                             backref=sqlalchemy.orm.backref('compounds'))
	_emprical_formula = sqlalchemy.Column(sqlalchemy.String())
	_mol_wt = sqlalchemy.Column(sqlalchemy.Float())
	_charge = sqlalchemy.Column(sqlalchemy.Integer())
	__tablename__ = 'compound'
	__mapper_args__ = {'polymorphic_identity': 'compound'}

	def get_inchi_structures(self):
		return [s.value for s in self.structures if s.format == 'inchi']	


class Metabolite(Compound):
	structure = openbabel.OBMol()

class RNA(Compound):
	structure = None

class DNA(Compound):
	structure = None

class Protein(Compound):
	structure = None

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
	#20192020
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

class ExternalReference(Base):
	namespace = sqlalchemy.Column(String()) # pubmed, pubchem, kegg.compound
	_id = sqlalchemy.Column(String(), primary_key=True) # unique identifier within namespace
	sqlalchemy.schema.UniqueConstraint(namespace, _id)
	__tablename__ = 'external_references'

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


'''
	TODO : BioPolymerForm()
'''