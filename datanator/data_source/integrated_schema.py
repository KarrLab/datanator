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

class Compartment(Entry):

	_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('entry._id'), primary_key=True)
	name = sqlalchemy.Column(sqlalchemy.String())

	__tablename__ = 'compartment'
	# __mapper_args__ = {'polymorphic_identity': 'compartment'}


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
	compound = sqlalchemy.orm.relationship('Compound', backref=sqlalchemy.orm.backref('reaction_participants'), foreign_keys=[compound_id])
	compartment_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('compartment._id'), index=True)
	compartment = sqlalchemy.orm.relationship('Compartment', backref=sqlalchemy.orm.backref('reaction_participants'), foreign_keys=[compartment_id])
	coefficient = sqlalchemy.Column(sqlalchemy.Float())
	_type = sqlalchemy.Column(sqlalchemy.String())

class Complex(Compound):
	#20192020
	structure = None