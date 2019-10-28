# -*- coding: utf-8 -*-

"""
:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-04
:Copyright: 2017, Karr Lab
:License: MIT
"""

from datanator.core import data_source
from datanator.util import molecule_util
from xml import etree
import Bio.Alphabet
import Bio.SeqUtils
import bs4
import csv
import datetime
import html
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
# :obj:`Base`: base model for local sqlite database

entry_synonym = sqlalchemy.Table(
    'entry_synonym', Base.metadata,
    sqlalchemy.Column('entry__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('entry._id'), index=True),
    sqlalchemy.Column('synonym__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('synonym._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Entry:Synonym many-to-many association table

entry_resource = sqlalchemy.Table(
    'entry_resource', Base.metadata,
    sqlalchemy.Column('entry__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('entry._id'), index=True),
    sqlalchemy.Column('resource__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('resource._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Entry:Resource many-to-many association table

compound_compound_structure = sqlalchemy.Table(
    'compound_compound_structure', Base.metadata,
    sqlalchemy.Column('compound__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('compound._id'), index=True),
    sqlalchemy.Column('compound_structure__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('compound_structure._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Compound:CompoundStructure many-to-many association table

kinetic_law_resource = sqlalchemy.Table(
    'kinetic_law_resource', Base.metadata,
    sqlalchemy.Column('kinetic_law__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('kinetic_law._id'), index=True),
    sqlalchemy.Column('resource__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('resource._id'), index=True),
)
# :obj:`sqlalchemy.Table`: KineticLaw:Resource many-to-many association table


class Synonym(Base):
    """ Represents a synonym to a SABIO-RK entry

    Attributes:
        name (:obj:`str`): name of the synonym
        entries (:obj:`list` of :obj:`Entry`): list of entries with the synonym
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    name = sqlalchemy.Column(sqlalchemy.String(), index=True)

    __tablename__ = 'synonym'


class Resource(Base):
    """ Represents an external resource

    Attributes:
        namespace (:obj:`str`): external namespace
        id (:obj:`str`): external identifier
        entries (:obj:`list` of :obj:`Entry`): entries
        kinetic_laws (:obj:`list` of :obj:`KineticLaw`): kinetic laws
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    namespace = sqlalchemy.Column(sqlalchemy.String())
    id = sqlalchemy.Column(sqlalchemy.String())

    sqlalchemy.schema.UniqueConstraint(namespace, id)

    __tablename__ = 'resource'


class Entry(Base):
    """ Represents a compartment in the SABIO-RK database

    Attributes:
        id (:obj:`int`): external identifier
        name (:obj:`str`): name
        synonyms (:obj:`list` of :obj:`Synonym`): list of synonyms
        cross_references (:obj:`list` of :obj:`Resource`): list of cross references
        created (:obj:`datetime.datetime`): date that the sqlite object was created
        updated (:obj:`datetime.datetime`): date that the sqlite object was last updated
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    _type = sqlalchemy.Column(sqlalchemy.String())

    id = sqlalchemy.Column(sqlalchemy.Integer(), index=True)
    name = sqlalchemy.Column(sqlalchemy.String(), index=True)
    synonyms = sqlalchemy.orm.relationship('Synonym', secondary=entry_synonym, backref=sqlalchemy.orm.backref('entries'))
    cross_references = sqlalchemy.orm.relationship('Resource', secondary=entry_resource, backref=sqlalchemy.orm.backref('entries'))
    created = sqlalchemy.Column(sqlalchemy.DateTime, default=datetime.datetime.utcnow())
    modified = sqlalchemy.Column(sqlalchemy.DateTime, onupdate=datetime.datetime.utcnow())

    sqlalchemy.schema.UniqueConstraint(id, _type)

    __tablename__ = 'entry'
    __mapper_args__ = {'polymorphic_on': _type}


class ReactionParticipant(Base):
    """ Represents a participant in a SABIO-RK reaction

    Attributes:
        compound (:obj:`Compound`): compound
        compartment (:obj:`Compartment`): compartment
        coefficient (:obj:`float`): coefficient
        type (:obj:`str`): type
        reactant_kinetic_law (:obj:`KineticLaw`): kinetic law in which the participant appears as a reactant
        product_kinetic_law (:obj:`KineticLaw`): kinetic law in which the participant appears as a product
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    compound_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('compound._id'), index=True)
    compound = sqlalchemy.orm.relationship('Compound', backref=sqlalchemy.orm.backref('reaction_participants'), foreign_keys=[compound_id])
    compartment_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('compartment._id'), index=True)
    compartment = sqlalchemy.orm.relationship('Compartment', backref=sqlalchemy.orm.backref(
        'reaction_participants'), foreign_keys=[compartment_id])
    coefficient = sqlalchemy.Column(sqlalchemy.Float())
    type = sqlalchemy.Column(sqlalchemy.String())
    reactant_kinetic_law_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('kinetic_law._id'), index=True)
    product_kinetic_law_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('kinetic_law._id'), index=True)
    modifier_kinetic_law_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('kinetic_law._id'), index=True)

    __tablename__ = 'reaction_participant'


class Parameter(Entry):
    """ Represents a parameter in the SABIO-RK database

    Attributes:
        kinetic_law (:obj:`KineticLaw`): kinetic law
        type (:obj:`int`): SBO term
        compound (:obj:`Compound`): compound
        enzyme (:obj:`Enzyme`): enzyme
        compartment (:obj:`Compartment`): compartment
        value (:obj:`float`): normalized value
        error (:obj:`float`): normalized error
        units (:obj:`str`): normalized units
        observed_name (:obj:`str`): name
        observed_type (:obj:`int`): SBO term
        observed_value (:obj:`float`): observed value
        observed_error (:obj:`float`): observed error
        observed_units (:obj:`str`): observed units

        TYPES (:obj:`dict` of :obj:`int`: :obj:`str`): dictionary of SBO terms and their canonical string symbols
        UNITS (:obj:`dict` of :obj:`int`: :obj:`str`): dictionary of SBO terms and their canonical units
    """

    TYPES = {
        # 9: 'k',
        # 16: 'k',
        # 17: 'k',
        25: 'k_cat',
        27: 'k_m',
        # 156: 'k_rev',
        186: 'v_max',
        # 190: n,
        # 191: n,
        261: 'k_i',
        # 281: 'k_eq',
        # 282: 'k_d',
        # 338: 'k_d',
        # 283: 'k_a',
        # 337: 'k_a',
        # 320: 'k_cat_p',
        # 321: 'k_cat_s',
        # 322: 'k_m_s',
        # 323: 'k_m_p',
        # 349: 'k_inact',
        # 363: 'k_x',
    }

    _id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('entry._id'), primary_key=True)

    kinetic_law_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('kinetic_law._id'), index=True)
    type = sqlalchemy.Column(sqlalchemy.Integer(), index=True)
    compound_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('compound._id'), index=True)
    compound = sqlalchemy.orm.relationship(
        'Compound', uselist=False, 
        backref=sqlalchemy.orm.backref('parameters'), 
        foreign_keys=[compound_id])
    enzyme_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('enzyme._id'), index=True)
    enzyme = sqlalchemy.orm.relationship(
        'Enzyme', uselist=False,
        backref=sqlalchemy.orm.backref('parameters'),
        foreign_keys=[enzyme_id])
    compartment_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('compartment._id'), index=True)
    compartment = sqlalchemy.orm.relationship(
        'Compartment', uselist=False,
        backref=sqlalchemy.orm.backref('parameters'),
        foreign_keys=[compartment_id])
    value = sqlalchemy.Column(sqlalchemy.Float())
    error = sqlalchemy.Column(sqlalchemy.Float())
    units = sqlalchemy.Column(sqlalchemy.String(), index=True)

    observed_name = sqlalchemy.Column(sqlalchemy.String())
    observed_type = sqlalchemy.Column(sqlalchemy.Integer())
    observed_value = sqlalchemy.Column(sqlalchemy.Float())
    observed_error = sqlalchemy.Column(sqlalchemy.Float())
    observed_units = sqlalchemy.Column(sqlalchemy.String())

    __tablename__ = 'parameter'
    __mapper_args__ = {'polymorphic_identity': 'parameter'}


class KineticLaw(Entry):
    """ Represents a kinetic law in the SABIO-RK database

    Attributes:
        reactants (:obj:`list` of :obj:`ReactionParticipant`): list of reactants
        products (:obj:`list` of :obj:`ReactionParticipant`): list of products
        enzyme (:obj:`Enzyme`): enzyme
        enzyme_compartment (:obj:`Compartment`): compartment
        enzyme_type (:obj:`str`): type of the enzyme (e.g. Modifier-Catalyst)
        tissue (:obj:`str`): tissue
        mechanism (:obj:`str`): mechanism of enzymatic catalysis (e.g. Michaelis-Menten)
        equation (:obj:`str`): equation
        parameters (:obj:`list` of :obj:`Parameter`): list of parameters
        modifiers (:obj:`list` of :obj:`ReactionParticipant`): list of modifiers
        taxon (:obj:`str`): taxon
        taxon_wildtype (:obj:`bool`): if :obj:`True`, the taxon represent the wild type
        taxon_variant (:obj:`str`): variant of the taxon
        temperature (:obj:`float`): temperature in C
        ph (:obj:`float`): pH
        media (:obj:`str`): media
        references (:obj:`list` of :obj:`Resource`): list of PubMed references
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('entry._id'), primary_key=True)
    reactants = sqlalchemy.orm.relationship('ReactionParticipant', backref=sqlalchemy.orm.backref('reactant_kinetic_law'),
                                            foreign_keys=[ReactionParticipant.reactant_kinetic_law_id],
                                            cascade='all, delete-orphan')
    products = sqlalchemy.orm.relationship('ReactionParticipant', backref=sqlalchemy.orm.backref('product_kinetic_law'),
                                           foreign_keys=[ReactionParticipant.product_kinetic_law_id],
                                           cascade='all, delete-orphan')
    enzyme_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('enzyme._id'), index=True)
    enzyme = sqlalchemy.orm.relationship('Enzyme', uselist=False, backref=sqlalchemy.orm.backref('kinetic_laws'), foreign_keys=[enzyme_id])
    enzyme_compartment_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('compartment._id'), index=True)
    enzyme_compartment = sqlalchemy.orm.relationship(
        'Compartment', uselist=False, backref=sqlalchemy.orm.backref('kinetic_laws'), foreign_keys=[enzyme_compartment_id])
    enzyme_type = sqlalchemy.Column(sqlalchemy.String())
    tissue = sqlalchemy.Column(sqlalchemy.String())
    mechanism = sqlalchemy.Column(sqlalchemy.String())
    equation = sqlalchemy.Column(sqlalchemy.Text())
    parameters = sqlalchemy.orm.relationship('Parameter', backref=sqlalchemy.orm.backref('kinetic_law'),
                                             foreign_keys=[Parameter.kinetic_law_id], cascade='all, delete-orphan')
    modifiers = sqlalchemy.orm.relationship('ReactionParticipant', backref=sqlalchemy.orm.backref('modifier_kinetic_law'),
                                            foreign_keys=[ReactionParticipant.modifier_kinetic_law_id],
                                            cascade='all, delete-orphan')
    taxon = sqlalchemy.Column(sqlalchemy.Integer())
    taxon_wildtype = sqlalchemy.Column(sqlalchemy.Boolean())
    taxon_variant = sqlalchemy.Column(sqlalchemy.UnicodeText())
    temperature = sqlalchemy.Column(sqlalchemy.Float())
    ph = sqlalchemy.Column(sqlalchemy.Float())
    media = sqlalchemy.Column(sqlalchemy.UnicodeText())
    references = sqlalchemy.orm.relationship('Resource', secondary=kinetic_law_resource, backref=sqlalchemy.orm.backref('kinetic_laws'))

    __tablename__ = 'kinetic_law'
    __mapper_args__ = {'polymorphic_identity': 'kinetic_law'}


class CompoundStructure(Base):
    """ Represents the structure of a compound and its format

    Attributes:
        compounds (:obj:`list` of :obj:`Compound`): list of compounds
        value (:obj:`str`): the structure in InChI, SMILES, etc. format
        format (:obj:`str`): format (InChI, SMILES, etc.) of the structure
        _value_inchi (:obj:`str`): structure in InChI format
        _value_inchi_formula_connectivity (:obj:`str`): empiral formula (without hydrogen) and connectivity InChI layers; used to
            quickly search for compound structures
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    value = sqlalchemy.Column(sqlalchemy.Text())
    format = sqlalchemy.Column(sqlalchemy.String())
    _value_inchi = sqlalchemy.Column(sqlalchemy.Text(), index=True)
    _value_inchi_formula_connectivity = sqlalchemy.Column(sqlalchemy.Text(), index=True)

    sqlalchemy.schema.UniqueConstraint(value, format)

    __tablename__ = 'compound_structure'

    def calc_inchi_formula_connectivity(self):
        """ Calculate a searchable structures

        * InChI format
        * Core InChI format

            * Formula layer (without hydrogen)
            * Connectivity layer
        """

        # if necessary, convert structure to InChI
        if self.format == 'inchi':
            self._value_inchi = self.value
        else:
            try:
                self._value_inchi = molecule_util.Molecule(structure=self.value).to_inchi() or None
            except ValueError:
                self._value_inchi = None

        # calculate formula (without hydrogen) and connectivity
        if self._value_inchi:
            self._value_inchi_formula_connectivity = molecule_util.InchiMolecule(self._value_inchi) \
                .get_formula_and_connectivity()


class Compound(Entry):
    """ Represents a compound in the SABIO-RK database

    Attributes:
        _is_name_ambiguous (:obj:`bool`): if :obj:`True`, the currently stored compound name should not be trusted
            because multiple names for the same compound have been discovered. The consensus name must be obtained
            using :obj:`download_compounds`
        structures (:obj:`list` of :obj:`CompoundStructure`): structures
        reaction_participants (:obj:`list` of :obj:`ReactionParticipant`): list of reaction participants
        parameters (:obj:`list` of :obj:`Parameter`): list of parameters
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('entry._id'), primary_key=True)

    _is_name_ambiguous = sqlalchemy.Column(sqlalchemy.Boolean(), default=False)
    structures = sqlalchemy.orm.relationship('CompoundStructure', secondary=compound_compound_structure,
                                             backref=sqlalchemy.orm.backref('compounds'))

    __tablename__ = 'compound'
    __mapper_args__ = {'polymorphic_identity': 'compound'}

    def get_inchi_structures(self):
        """ Get InChI-formatted structures

        Returns:
            :obj:`list` of :obj:`str`: list of structures in InChI format
        """
        return [s.value for s in self.structures if s.format == 'inchi']

    def get_smiles_structures(self):
        """ Get SMILES-formatted structures

        Returns:
            :obj:`list` of :obj:`str`: list of structures in SMILES format
        """
        return [s.value for s in self.structures if s.format == 'smiles']


class Enzyme(Entry):
    """ Represents an enzyme in the SABIO-RK database

    Attributes:
        subunits (:obj:`list` of :obj:`EnzymeSubunit`): list of subunits
        kinetic_laws (:obj:`list` of :obj:`KineticLaw`): list of kinetic laws
        molecular_weight (:obj:`float`): molecular weight in Daltons
        parameters (:obj:`list` of :obj:`Parameter`): list of parameters
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('entry._id'), primary_key=True)
    molecular_weight = sqlalchemy.Column(sqlalchemy.Float())

    __tablename__ = 'enzyme'
    __mapper_args__ = {'polymorphic_identity': 'enzyme'}


class EnzymeSubunit(Entry):
    """ Represents an enzyme in the SABIO-RK database

    Attributes:
        enzyme (:obj:`Enzyme`): enzyme
        coefficient (:obj:`int`): stoichiometry of the subunit in the enzyme
        sequence (:obj:`str`): amino acid sequence
        molecular_weight (:obj:`float`): molecular weight in Daltons
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('entry._id'), primary_key=True)
    enzyme_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('enzyme._id'), index=True)
    enzyme = sqlalchemy.orm.relationship('Enzyme', uselist=False, backref=sqlalchemy.orm.backref(
        'subunits', cascade='all, delete-orphan'), foreign_keys=[enzyme_id])
    coefficient = sqlalchemy.Column(sqlalchemy.Integer())
    sequence = sqlalchemy.Column(sqlalchemy.Text())
    molecular_weight = sqlalchemy.Column(sqlalchemy.Float())

    __tablename__ = 'enzyme_subunit'
    __mapper_args__ = {'polymorphic_identity': 'enzyme_subunit'}


class Compartment(Entry):
    """ Represents a compartment in the SABIO-RK database

    Attributes:
        kinetic_laws (:obj:`list` of :obj:`KineticLaw`): list of kinetic laws
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('entry._id'), primary_key=True)

    __tablename__ = 'compartment'
    __mapper_args__ = {'polymorphic_identity': 'compartment'}


class SabioRk(data_source.HttpDataSource):
    """ A local sqlite copy of the SABIO-RK database

    Attributes:
        webservice_batch_size (:obj:`int`): default size of batches to download kinetic information from the SABIO webservice. Note:
            this should be set to one because SABIO exports units incorrectly when multiple kinetic laws are requested
        excel_batch_size (:obj:`int`): default size of batches to download kinetic information from the SABIO
            Excel download service

        ENDPOINT_KINETIC_LAWS_SEARCH (:obj:`str`): URL to obtain a list of the ids of all of the kinetic laws in SABIO-Rk
        ENDPOINT_WEBSERVICE (:obj:`str`): URL for the SABIO-RK webservice
        ENDPOINT_EXCEL_EXPORT (:obj:`str`): URL to download kinetic data as a table in TSV format
        ENDPOINT_COMPOUNDS_PAGE (:obj:`str`): URL to download information about a SABIO-RK compound
        SKIP_KINETIC_LAW_IDS (:obj:`tuple` of :obj:`int`): IDs of kinetic laws that should be skipped (because they cannot contained
            errors and can't be downloaded from SABIO)
        PUBCHEM_MAX_TRIES (:obj:`int`): maximum number of times to time querying PubChem before failing
        PUBCHEM_TRY_DELAY (:obj:`float`): delay in seconds between PubChem queries (to delay overloading the server)
    """

    base_model = Base
    ENDPOINT_DOMAINS = {
        'sabio_rk': 'http://sabiork.h-its.org',
        'uniprot': 'http://www.uniprot.org',
    }
    ENDPOINT_KINETIC_LAWS_SEARCH = ENDPOINT_DOMAINS['sabio_rk'] + '/sabioRestWebServices/searchKineticLaws/entryIDs'
    ENDPOINT_WEBSERVICE = ENDPOINT_DOMAINS['sabio_rk'] + '/sabioRestWebServices/kineticLaws'
    ENDPOINT_EXCEL_EXPORT = ENDPOINT_DOMAINS['sabio_rk'] + '/entry/exportToExcelCustomizable'
    ENDPOINT_COMPOUNDS_PAGE = ENDPOINT_DOMAINS['sabio_rk'] + '/compdetails.jsp'
    ENDPOINT_KINETIC_LAWS_PAGE = ENDPOINT_DOMAINS['sabio_rk'] + '/kindatadirectiframe.jsp'
    SKIP_KINETIC_LAW_IDS = (51286,)
    PUBCHEM_MAX_TRIES = 10
    PUBCHEM_TRY_DELAY = 0.25

    def __init__(self, name=None, cache_dirname=None, clear_content=False, load_content=False, max_entries=float('inf'),
                 commit_intermediate_results=False, download_backups=True, verbose=False,
                 clear_requests_cache=False, download_request_backup=False,
                 webservice_batch_size=1, excel_batch_size=100,quilt_owner=None, quilt_package=None): #
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
                                      quilt_owner=quilt_owner, quilt_package=quilt_package) #

    def load_content(self):
        """ Download the content of SABIO-RK and store it to a local sqlite database. """

        ##################################
        ##################################
        # determine ids of kinetic laws
        if self.verbose:
            print('Downloading the IDs of the kinetic laws ...')

        ids = self.load_kinetic_law_ids()

        if self.verbose:
            print('  Downloaded {} IDs'.format(len(ids)))

        ##################################
        ##################################
        # remove bad IDs
        ids = list(filter(lambda id: id not in self.SKIP_KINETIC_LAW_IDS, ids))

        # sort ids
        ids.sort()

        # load only `max_entries` IDs
        if len(ids) > self.max_entries:
            ids = ids[0:self.max_entries]

        ##################################
        ##################################
        # download kinetic laws
        new_ids = list(set(ids).difference(set(l.id for l in self.session.query(KineticLaw).all())))
        new_ids.sort()

        if self.verbose:
            print('Downloading {} kinetic laws ...'.format(len(new_ids)))

        self.load_kinetic_laws(new_ids)

        if self.verbose:
            print('  done')

        ##################################
        ##################################
        # download compounds
        compounds = self.session.query(Compound).filter(~Compound.structures.any()).order_by(Compound.id).all()

        if self.verbose:
            print('Downloading {} compounds ...'.format(len(compounds)))

        self.load_compounds(compounds)

        if self.verbose:
            print('  done')

        ##################################
        ##################################
        # fill in missing information from Excel export
        loaded_new_ids = list(set(new_ids).intersection(set(l.id for l in self.session.query(KineticLaw).all())))
        loaded_new_ids.sort()

        if self.verbose:
            print('Updating {} kinetic laws ...'.format(len(loaded_new_ids)))

        self.load_missing_kinetic_law_information_from_tsv(loaded_new_ids)

        if self.verbose:
            print('  done')

        ##################################
        ##################################
        # infer structures for compounds with no provided structure
        compounds = self.session.query(Compound).filter(~Compound.structures.any()).all()

        if self.verbose:
            print('Inferring structures for {} compounds ...'.format(len(compounds)))

        self.infer_compound_structures_from_names(compounds)

        if self.verbose:
            print('  done')

        ##################################
        ##################################
        # normalize compound structures to facilitate seaching. retain only
        # - InChI formula layer (without hydrogen)
        # - InChI connectivity layer
        compound_structures = self.session.query(CompoundStructure).filter_by(_value_inchi_formula_connectivity=None).all()

        if self.verbose:
            print('Calculating searchable structures for {} structures ...'.format(len(compound_structures)))

        for i_compound_structure, compound_structure in enumerate(compound_structures):
            if self.verbose and (i_compound_structure % 100 == 0):
                print('  Calculating searchable structure for compound {} of {}'.format(
                    i_compound_structure + 1, len(compound_structures)))
            compound_structure.calc_inchi_formula_connectivity()

        if self.verbose:
            print('  done')

        ##################################
        ##################################
        # fill in missing information from HTML pages
        if self.verbose:
            print('Updating {} kinetic laws ...'.format(len(loaded_new_ids)))

        self.load_missing_enzyme_information_from_html(loaded_new_ids)

        if self.verbose:
            print('  done')

        ##################################
        ##################################
        # calculate enzyme molecular weights
        enzymes = self.session \
            .query(Enzyme) \
            .filter_by(molecular_weight=None) \
            .all()

        if self.verbose:
            print('Calculating {} enzyme molecular weights ...'.format(len(enzymes)))

        self.calc_enzyme_molecular_weights(enzymes)

        if self.verbose:
            print('  done')

        ##################################
        ##################################
        if self.verbose:
            print('Normalizing {} parameter values ...'.format(len(loaded_new_ids)))

        self.normalize_kinetic_laws(loaded_new_ids)

        if self.verbose:
            print('  done')

        ##################################
        ##################################
        # commit changes
        self.session.commit()

        # calculate statistics
        self.export_stats(self.calc_stats())

    def load_kinetic_law_ids(self):
        """ Download the IDs of all of the kinetic laws stored in SABIO-RK

        Returns:
            :obj:`list` of :obj:`int`: list of kinetic law IDs

        Raises:
            :obj:`Error`: if an HTTP request fails or the expected number of kinetic laws is not returned
        """
        # create session
        session = self.requests_session
        response = session.get(self.ENDPOINT_KINETIC_LAWS_SEARCH, params={
            'q': 'DateSubmitted:01/01/2000',
        })
        response.raise_for_status()

        # get IDs of kinetic laws
        root = etree.ElementTree.fromstring(response.text)
        ids = [int(float(node.text)) for node in root.findall('SabioEntryID')]

        # sort ids
        ids.sort()

        # return IDs
        return ids

    def load_kinetic_laws(self, ids):
        """ Download kinetic laws from SABIO-RK

        Args:
            ids (:obj:`list` of :obj:`int`): list of IDs of kinetic laws to download

        Raises:
            :obj:`Error`: if an HTTP request fails
        """
        # todo: scrape strain, recombinant, experiment type information from web pages
        session = self.requests_session

        batch_size = self.webservice_batch_size

        for i_batch in range(int(math.ceil(float(len(ids)) / batch_size))):
            if self.verbose and (i_batch % max(1, 100. / batch_size) == 0):
                print('  Downloading kinetic laws {}-{} of {} in SBML format'.format(
                    i_batch * batch_size + 1,
                    min(len(ids), i_batch * batch_size + max(100, batch_size)),
                    len(ids)))

            batch_ids = ids[i_batch * batch_size:min((i_batch + 1) * batch_size, len(ids))]
            response = session.get(self.ENDPOINT_WEBSERVICE, params={
                'kinlawids': ','.join(str(id) for id in batch_ids),
            })

            response.raise_for_status()
            if not response.text:
                cache = session.cache
                key = cache.create_key(response.request)
                cache.delete(key)
                raise Exception('Unable to download kinetic laws with ids {}'.format(', '.join([str(id) for id in batch_ids])))

            self.create_kinetic_laws_from_sbml(batch_ids, response.content if six.PY2 else response.text)

            if self.commit_intermediate_results:
                self.session.commit()

        # print warning with list of unidentified ids
        loaded_ids = [l.id for l in self.session.query(KineticLaw).order_by(KineticLaw.id)]

        not_loaded_ids = list(set(ids).difference(loaded_ids))
        if not_loaded_ids:
            not_loaded_ids.sort()
            warnings.warn('Several kinetic laws were not found:\n  {}'.format(
                '\n  '.join([str(id) for id in not_loaded_ids])), data_source.DataSourceWarning)

    def create_kinetic_laws_from_sbml(self, ids, sbml):
        """ Add kinetic laws defined in an SBML file to the local sqlite database

        Args:
            ids (:obj:`list` of :obj:`int`): list kinetic law IDs
            sbml (:obj:`str`): SBML representation of one or more kinetic laws

        Returns:
            :obj:`tuple`:

                * :obj:`list` of :obj:`KineticLaw`: list of kinetic laws
                * :obj:`list` of :obj:`Compound` or :obj:`Enzyme`: list of species (compounds or enzymes)
                * :obj:`list` of :obj:`Compartment`: list of compartments
        """
        reader = libsbml.SBMLReader()
        doc = reader.readSBMLFromString(sbml)
        model = doc.getModel()

        functions = {}
        functions_sbml = model.getListOfFunctionDefinitions()
        for i_function in range(functions_sbml.size()):
            function_sbml = functions_sbml.get(i_function)
            math_sbml = function_sbml.getMath()
            if math_sbml.isLambda() and math_sbml.getNumChildren():
                eq = libsbml.formulaToL3String(math_sbml.getChild(math_sbml.getNumChildren() - 1))
            else:
                eq = None
            if eq in ('', 'NaN'):
                eq = None
            functions[function_sbml.getId()] = eq

        units = {}
        units_sbml = model.getListOfUnitDefinitions()
        for i_unit in range(units_sbml.size()):
            unit_sbml = units_sbml.get(i_unit)
            units[unit_sbml.getId()] = unit_sbml.getName()

        # compartments
        compartments_sbml = model.getListOfCompartments()
        compartments = []
        for i_compartment in range(compartments_sbml.size()):
            compartment_sbml = compartments_sbml.get(i_compartment)
            compartments.append(self.create_compartment_from_sbml(compartment_sbml))

        # species
        specie_properties = {}
        species_sbml = model.getListOfSpecies()
        species = []
        for i_specie in range(species_sbml.size()):
            specie_sbml = species_sbml.get(i_specie)
            specie, properties = self.create_specie_from_sbml(specie_sbml)
            species.append(specie)
            specie_properties[specie_sbml.getId()] = properties

        # kinetic laws
        reactions_sbml = model.getListOfReactions()
        if reactions_sbml.size() != len(ids):
            raise ValueError('{} reactions {} is different from the expected {}'.format(reaction_sbml.size(), len(ids)))
        kinetic_laws = []
        for i_reaction, id in enumerate(ids):
            reaction_sbml = reactions_sbml.get(i_reaction)
            kinetic_law = self.create_kinetic_law_from_sbml(
                id, reaction_sbml, specie_properties, functions, units)
            kinetic_laws.append(kinetic_law)

        return (kinetic_laws, species, compartments)

    def create_compartment_from_sbml(self, sbml):
        """ Add a compartment to the local sqlite database

        Args:
            sbml (:obj:`libsbml.Compartment`): SBML-representation of a compartment

        Returns:
            :obj:`Compartment`: compartment
        """
        name = sbml.getName()
        if name == 'Cell':
            return None

        query = self.session.query(Compartment).filter_by(name=name)
        if query.count():
            compartment = query.first()
        else:
            compartment = Compartment(name=name)
            self.session.add(compartment)
            compartment.modified = datetime.datetime.utcnow()

        return compartment

    def create_specie_from_sbml(self, sbml):
        """ Add a species to the local sqlite database

        Args:
            sbml (:obj:`libsbml.Species`): SBML-representation of a compound or enzyme

        Returns:
            :obj:`tuple`:

                * :obj:`Compound`: or :obj:`Enzyme`: compound or enzyme
                * :obj:`dict`: additional properties of the compound/enzyme

                    * `is_wildtype` (:obj:`bool`): indicates if the enzyme is wildtype or mutant
                    * `variant` (:obj:`str`): description of the variant of the eznyme
                    * `modifier_type` (:obj:`str`): type of the enzyme (e.g. Modifier-Catalyst)

        Raises:
            :obj:`ValueError`: if a species is of an unsupported type (i.e. not a compound or enzyme)
        """
        # id, name
        type_id_compartment = sbml.getId().split('_')
        type = type_id_compartment[0]
        id = int(float(type_id_compartment[1]))

        # modifier type
        modifier_type = ''
        if sbml.isSetAnnotation() \
                and sbml.getAnnotation().hasChild('sabiork') \
                and sbml.getAnnotation().getChild('sabiork').hasChild('modifierType'):
            modifier_type = sbml \
                .getAnnotation() \
                .getChild('sabiork') \
                .getChild('modifierType') \
                .getChild(0) \
                .getCharacters()

        # create object or return existing object
        if type == 'SPC':
            name = sbml.getName()
            properties = {'modifier_type': modifier_type}
            query = self.session.query(Compound).filter_by(id=id)
            if query.count():
                specie = query.first()
            else:
                specie = Compound()
                self.session.add(specie)
        elif type == 'ENZ':
            name, is_wildtype, variant = self.parse_enzyme_name(sbml.getName())
            if six.PY2:
                variant = unicode(variant.decode('utf-8'))
            properties = {'is_wildtype': is_wildtype, 'variant': variant, 'modifier_type': modifier_type}

            query = self.session.query(Enzyme).filter_by(id=id)
            if query.count():
                specie = query.first()
            else:
                specie = Enzyme()
                self.session.add(specie)
        else:
            raise ValueError('Unsupported species type: {}'.format(type))

        # set properties
        specie.id = id
        if specie.name is None:
            specie.name = name
        elif specie.name != name:
            specie._is_name_ambiguous = True

        # cross references
        cross_references = self.create_cross_references_from_sbml(sbml)
        if type == 'SPC':
            specie.cross_references = cross_references
        elif type == 'ENZ':
            specie.subunits = []
            specie.cross_references = []
            for cross_reference in cross_references:
                if cross_reference.namespace == 'uniprot':
                    specie.subunits.append(EnzymeSubunit(cross_references=[cross_reference]))
                else:
                    specie.cross_references.append(cross_reference)

        # updated
        specie.modified = datetime.datetime.utcnow()

        return (specie, properties)

    def parse_enzyme_name(self, sbml):
        """ Parse the name of an enzyme in SBML for the enzyme name, wild type status, and variant
        description that it contains.

        Args:
            sbml (:obj:`str`): enzyme name in SBML

        Returns:
            :obj:`tuple`:

                * :obj:`str`: name
                * :obj:`bool`: if :obj:`True`, the enzyme is wild type
                * :obj:`str`: variant

        Raises:
            :obj:`ValueError`: if the enzyme name is formatted in an unsupport format
        """
        match = re.match(r'^(.*?)\(Enzyme\) (wildtype|mutant),?(.*?)$', sbml, re.IGNORECASE)
        if match:
            name = match.group(1)
            is_wildtype = match.group(2).lower() == 'wildtype'
            variant = match.group(3).strip()
            return (name, is_wildtype, variant)

        match = re.match(r'^Enzyme (wildtype|mutant),?( (.*?))*$', sbml, re.IGNORECASE)
        if match:
            if match.group(3):
                name = match.group(3).strip()
            else:
                name = None
            is_wildtype = match.group(1).lower() == 'wildtype'
            variant = None
            return (name, is_wildtype, variant)

        match = re.match(r'^Enzyme - *$', sbml, re.IGNORECASE)
        if match:
            name = None
            is_wildtype = True
            variant = None
            return (name, is_wildtype, variant)

        raise ValueError('Cannot parse enzyme name: {}'.format(sbml))

    def create_kinetic_law_from_sbml(self, id, sbml, specie_properties, functions, units):
        """ Add a kinetic law to the local sqlite database

        Args:
            id (:obj:`int`): identifier
            sbml (:obj:`libsbml.KineticLaw`): SBML-representation of a reaction
            specie_properties (:obj:`dict`): additional properties of the compounds/enzymes

                * `is_wildtype` (:obj:`bool`): indicates if the enzyme is wildtype or mutant
                * `variant` (:obj:`str`): description of the variant of the eznyme
                * `modifier_type` (:obj:`str`): type of the enzyme (e.g. Modifier-Catalyst)

            functions (:obj:`dict` of :obj:`str`: :obj:`str`): dictionary of rate law equations (keys = IDs in SBML, values = equations)
            units (:obj:`dict` of :obj:`str`: :obj:`str`): dictionary of units (keys = IDs in SBML, values = names)

        Returns:
            :obj:`KineticLaw`: kinetic law

        Raises:
            :obj:`ValueError`: if the temperature is expressed in an unsupported unit
        """
        law = sbml.getKineticLaw()
        x_refs = self.create_cross_references_from_sbml(law)
        reaction_x_refs = self.create_cross_references_from_sbml(sbml)

        # stop if kinetic law entry is empty
        if not law.getMetaId():
            return None

        # ID
        annotated_id = next((int(float(x_ref.id)) for x_ref in x_refs if x_ref.namespace == 'sabiork.kineticrecord'), None)
        if annotated_id is not None and annotated_id != id:
            raise ValueError('Annotated ID {} is different from expected ID {}'.format(annotated_id, id))
        query = self.session.query(KineticLaw).filter_by(id=id)
        if query.count():
            kinetic_law = query.first()
        else:
            kinetic_law = KineticLaw(id=id)
            self.session.add(kinetic_law)

        """ participants """
        kinetic_law.reactants[:] = []
        reactants = sbml.getListOfReactants()
        for i_part in range(reactants.size()):
            part_sbml = reactants.get(i_part)
            compound, compartment = self.get_specie_reference_from_sbml(part_sbml.getSpecies())
            part = ReactionParticipant(
                compound=compound,
                compartment=compartment,
                coefficient=part_sbml.getStoichiometry())
            self.session.add(part)
            kinetic_law.reactants.append(part)

        kinetic_law.products[:] = []
        products = sbml.getListOfProducts()
        for i_part in range(products.size()):
            part_sbml = products.get(i_part)
            compound, compartment = self.get_specie_reference_from_sbml(part_sbml.getSpecies())
            part = ReactionParticipant(
                compound=compound,
                compartment=compartment,
                coefficient=part_sbml.getStoichiometry())
            self.session.add(part)
            kinetic_law.products.append(part)

        """ cross references """
        # Note: these are stored KineticLaws rather than under Reactions because this seems to how SABIO-RK stores this information.
        # For example, kinetic laws 16016 and 28003 are associated with reaction 9930, but they have different EC numbers 1.1.1.52 and
        # 1.1.1.50, respectively.
        kinetic_law.cross_references = list(filter(lambda x_ref: x_ref.namespace not in ['taxonomy'], reaction_x_refs))

        # rate_law
        kinetic_law.equation = functions[law.getMetaId()[5:]]

        # parameters
        kinetic_law.parameters = []
        params = law.getListOfLocalParameters()
        for i_param in range(params.size()):
            param = params.get(i_param)

            match = re.match(r'^(.*?)_((SPC|ENZ)_([0-9]+)_(.*?))$', param.getId(), re.IGNORECASE)
            if match:
                observed_name = match.group(1)
                species, compartment = self.get_specie_reference_from_sbml(match.group(2))
                if isinstance(species, Compound):
                    compound = species
                    enzyme = None
                else:
                    compound = None
                    enzyme = species
            else:
                observed_name = param.getId()
                compound = None
                enzyme = None
                compartment = None
            observed_name = observed_name.replace('div', '/')

            observed_type = param.getSBOTerm()

            observed_units_id = param.getUnits()
            if observed_units_id:
                if observed_units_id in units:
                    observed_units = units[observed_units_id]
                else:
                    observed_units = observed_units_id
            else:
                observed_units = None

            observed_value = param.getValue()

            parameter = Parameter(
                compound=compound,
                enzyme=enzyme,
                compartment=compartment,
                observed_name=observed_name,
                observed_type=observed_type,
                observed_value=observed_value,
                observed_units=observed_units,
                modified=datetime.datetime.utcnow(),
            )
            self.session.add(parameter)
            kinetic_law.parameters.append(parameter)

        # modifiers to kinetic law
        kinetic_law.modifiers[:] = []
        modifiers = sbml.getListOfModifiers()
        for i_modifier in range(modifiers.size()):
            modifier = modifiers.get(i_modifier)
            modifier_id = modifier.getSpecies()
            specie, compartment = self.get_specie_reference_from_sbml(modifier_id)
            type = specie_properties[modifier.getSpecies()]['modifier_type']
            if modifier.getSpecies()[0:3] == 'SPC':
                part = ReactionParticipant(
                    compound=specie,
                    compartment=compartment,
                    type=type,
                )
                self.session.add(part)
                kinetic_law.modifiers.append(part)
            elif modifier_id[0:3] == 'ENZ':
                kinetic_law.enzyme, kinetic_law.enzyme_compartment = self.get_specie_reference_from_sbml(modifier_id)
                kinetic_law.enzyme_type = specie_properties[modifier.getSpecies()]['modifier_type']
                kinetic_law.taxon_wildtype = specie_properties[modifier_id]['is_wildtype']
                kinetic_law.taxon_variant = specie_properties[modifier_id]['variant']

        # taxon
        kinetic_law.taxon = next((int(float(x_ref.id)) for x_ref in reaction_x_refs if x_ref.namespace == 'taxonomy'), None)

        """ conditions """
        conditions = law \
            .getAnnotation() \
            .getChild('sabiork') \
            .getChild('experimentalConditions')

        # temperature
        if conditions.hasChild('temperature'):
            temperature = conditions \
                .getChild('temperature') \
                .getChild('startValueTemperature') \
                .getChild(0) \
                .getCharacters()
            temperature = float(temperature)
            temperature_units = conditions \
                .getChild('temperature') \
                .getChild('temperatureUnit') \
                .getChild(0) \
                .getCharacters()
            if temperature_units not in ['°C', '��C']:
                raise ValueError('Unsupported temperature units: {}'.format(temperature_units))
            kinetic_law.temperature = temperature

        # pH
        if conditions.hasChild('pH'):
            ph = conditions \
                .getChild('pH') \
                .getChild('startValuepH') \
                .getChild(0) \
                .getCharacters()
            kinetic_law.ph = float(ph)

        # media
        if conditions.hasChild('buffer'):
            media = conditions \
                .getChild('buffer') \
                .getChild(0) \
                .getCharacters() \
                .strip()
            if six.PY2:
                media = unicode(media.decode('utf-8'))
            kinetic_law.media = media

        """ references """
        kinetic_law.references = list(filter(lambda x_ref: x_ref.namespace != 'sabiork.kineticrecord', x_refs))

        """ updated """
        kinetic_law.modified = datetime.datetime.utcnow()

        return kinetic_law

    def normalize_kinetic_laws(self, ids):
        """ Normalize parameter values

        Args:
            ids (:obj:`list` of :obj:`int`): list of IDs of kinetic laws to download
        """
        for i_law, law in enumerate(self.session.query(KineticLaw).filter(KineticLaw.id.in_(ids)).all()):
            if self.verbose and (i_law % 100 == 0):
                print('  Normalizing kinetic law {} of {}'.format(i_law + 1, len(ids)))

            if law.enzyme:
                enzyme_molecular_weight = law.enzyme.molecular_weight
            else:
                enzyme_molecular_weight = None

            for p in law.parameters:
                p.name, p.type, p.value, p.error, p.units = self.normalize_parameter_value(
                    p.observed_name, p.observed_type, p.observed_value, p.observed_error, p.observed_units,
                    enzyme_molecular_weight)

        if self.commit_intermediate_results:
            self.session.commit()

    def normalize_parameter_value(self, name, type, value, error, units, enzyme_molecular_weight):
        """
        Args:
            name (:obj:`str`): parameter name
            type (:obj:`int`) parameter type (SBO term id)
            value (:obj:`float`): observed value
            error (:obj:`float`): observed error
            units (:obj:`str`): observed units
            enzyme_molecular_weight (:obj:`float`): enzyme molecular weight

        Returns:
            :obj:`tuple` of :obj:`str`, :obj:`int`, :obj:`float`, :obj:`float`, :obj:`str`: normalized name and
                its type (SBO term), value, error, and units

        Raises:
            :obj:`ValueError`: if :obj:`units` is not a supported unit of :obj:`type`
        """

        if type not in Parameter.TYPES:
            return (None, None, None, None, None)

        if value is None:
            return (None, None, None, None, None)

        type_name = Parameter.TYPES[type]

        if type_name == 'k_cat':
            if units in ['s^(-1)', 'mol*s^(-1)*mol^(-1)']:
                return ('k_cat', 25, value, error, 's^(-1)')

            if units in ['katal', 'katal_base']:
                # cannot be converted without knowing the enzyme amount in the measurement
                return (None, None, None, None, None)

            if units in ['M^(-1)*s^(-1)']:
                # off by mol^(2) * liter^(-1)
                return (None, None, None, None, None)

            if units in ['s^(-1)*g^(-1)', 'mol*s^(-1)*g^(-1)', 'M']:
                return (None, None, None, None, None)

            if units is None:
                return (None, None, None, None, None)

        elif type_name == 'k_m':
            if units in ['M']:
                return ('k_m', 27, value, error, 'M')

            if units in ['mol']:
                # off by liter^(-1)
                return (None, None, None, None, None)

            if units in ['mg/ml', 'M^2', 'mol/mol', 'katal*g^(-1)', 's^(-1)', 'mol*s^(-1)*g^(-1)',
                'l*g^(-1)*s^(-1)', 'M^(-1)*s^(-1)', 'M^(-1)']:
                return (None, None, None, None, None)

            if units is None:
                return (None, None, None, None, None)

        elif type_name == 'v_max':
            if units in ['mol*s^(-1)*g^(-1)', 'katal*g^(-1)']:
                if enzyme_molecular_weight:
                    f = enzyme_molecular_weight
                    return ('k_cat', 25, value * float(f) if value else None, error * float(f) if error else None, 's^(-1)')
                else:
                    return ('v_max', 186, value, error, 'mol*s^(-1)*g^(-1)')

            if units in ['katal*mol^(-1)', 'mol*s^(-1)*mol^(-1)', 'Hz', 'M*s^(-1)*M^(-1)', 's^(-1)', 'g/(s*g)']:
                return ('k_cat', 25, value, error, 's^(-1)')

            if units in ['mol/s', 'katal', 'katal_base']:
                # cannot be converted without knowing the enzyme amount in the measurement
                return (None, None, None, None, None)

            if units in ['M*s^(-1)*g^(-1)']:
                # has incorrect dimensionality from v_max by factor of liter^(-1)
                return (None, None, None, None, None)

            if units in ['M*s^(-1)', 'katal*l^(-1)']:
                # has incorrect dimensionality from k_cat by factor of liter^(-1)
                return (None, None, None, None, None)

            if units in ['mol*s^(-1)*m^(-1)', 'M', 'g', 'g/(l*s)', 'M^2', 'katal*s^(-1)', 'mol*g^(-1)', 'mol/(sec*m^2)', 'mg/ml',
                         'l*g^(-1)*s^(-1)', 'mol/(s*M)', 'katal*M^(-1)*g^(-1)', 'M^(-1)*s^(-1)', 'mol*s^',
                         's^(-1)*g^(-1)']:
                return (None, None, None, None, None)

            if units is None:
                return (None, None, None, None, None)

        elif type_name == 'k_i':
            if units in ['M']:
                return ('k_i', 261, value, error, 'M')

            if units in ['mol/mol', 'M^2', 'g', 'M^(-1)*s^(-1)', 'mol*s^(-1)*g^(-1)']:
                return (None, None, None, None, None)

            if units is None:
                return (None, None, None, None, None)

        raise ValueError('Unsupported units "{}" for parameter type {}'.format(units, type_name))

    def create_cross_references_from_sbml(self, sbml):
        """ Add cross references to the local sqlite database for an SBML object

        Args:
            sbml (:obj:`libsbml.SBase`): object in an SBML documentation

        Returns:
            :obj:`list` of :obj:`Resource`: list of resources
        """
        if not sbml.isSetAnnotation():
            return []

        xml = sbml.getAnnotation().getChild('RDF').getChild('Description')

        attr = libsbml.XMLTriple('resource', 'http://www.w3.org/1999/02/22-rdf-syntax-ns#', 'rdf')

        x_refs = []
        for i_child in range(xml.getNumChildren()):
            bag = xml.getChild(i_child).getChild('Bag')
            for i_li in range(bag.getNumChildren()):
                li = bag.getChild(i_li)

                val = li.getAttrValue(attr)
                if val.startswith('http://identifiers.org/ec-code/Complex '):
                    _, _, ids = val.partition(' ')
                    resources = [('ec-code', id) for id in ids.split('/')]
                else:
                    parsed_url = val.split('/')
                    if len(parsed_url) == 5:
                        resources = [(parsed_url[3], parsed_url[4])]
                    else:
                        resources = [(parsed_url[2], parsed_url[3])]

                for namespace, id in resources:
                    query = self.session.query(Resource).filter_by(namespace=namespace, id=id)
                    if query.count():
                        resource = query.first()
                    else:
                        resource = Resource(namespace=namespace, id=id)
                        self.session.add(resource)

                if resource not in x_refs:
                    x_refs.append(resource)

        return x_refs

    def get_specie_reference_from_sbml(self, specie_id):
        """ Get the compound/enzyme associated with an SBML species by its ID

        Args:
            specie_id (:obj:`str`): ID of an SBML species

        Returns:
            :obj:`tuple`:

                * :obj:`Compound` or :obj:`Enzyme`: compound or enzyme
                * :obj:`Compartment`: compartment

        Raises:
            :obj:`ValueError`: if the species is not a compound or enzyme, no species
                with `id` = `specie_id` exists, or no compartment with `name` = `compartment_name`
                exists
        """
        tmp = specie_id.split('_')
        type = tmp[0]
        specie_id = int(float(tmp[1]))
        compartment_name = '_'.join(tmp[2:])

        if type == 'SPC':
            q = self.session.query(Compound).filter_by(id=specie_id)
        elif type == 'ENZ':
            q = self.session.query(Enzyme).filter_by(id=specie_id)
        else:
            raise ValueError('Unsupported species type {}'.format(type))

        if q.count() != 1:
            raise ValueError('Could not find species with id {}'.format(specie_id))
        specie = q.first()

        if compartment_name != 'Cell':
            q = self.session.query(Compartment).filter_by(name=compartment_name)
            if q.count() != 1:
                raise ValueError('Could not find compartment with name "{}"'.format(compartment_name))
            compartment = q.first()
        else:
            compartment = None

        return (specie, compartment)

    def load_missing_kinetic_law_information_from_tsv(self, ids):
        """ Update the properties of kinetic laws in the local sqlite database based on content downloaded
        from SABIO in TSV format.

        Args:
            ids (:obj:`list` of :obj:`int`): list of IDs of kinetic laws to download
        """
        session = self.requests_session

        batch_size = self.excel_batch_size

        for i_batch in range(int(math.ceil(float(len(ids)) / batch_size))):
            if self.verbose:
                print('  Downloading kinetic laws {}-{} of {} in Excel format'.format(
                    i_batch * batch_size + 1,
                    min(len(ids), (i_batch + 1) * batch_size),
                    len(ids)))

            batch_ids = ids[i_batch * batch_size:min((i_batch + 1) * batch_size, len(ids))]
            response = session.get(self.ENDPOINT_EXCEL_EXPORT, params={
                'entryIDs[]': batch_ids,
                'fields[]': [
                    'EntryID',
                    'KineticMechanismType',
                    'Tissue',
                    'Parameter',
                ],
                'preview': False,
                'format': 'tsv',
                'distinctRows': 'false',
            })
            response.raise_for_status()
            if not response.text:
                cache = session.cache
                key = cache.create_key(response.request)
                cache.delete(key)
                raise Exception('Unable to download kinetic laws with ids {}'.format(', '.join([str(id) for id in batch_ids])))

            self.load_missing_kinetic_law_information_from_tsv_helper(response.text)

            if self.commit_intermediate_results:
                self.session.commit()

    def load_missing_kinetic_law_information_from_tsv_helper(self, tsv):
        """ Update the properties of kinetic laws in the local sqlite database based on content downloaded
        from SABIO in TSV format.

        Note: this method is necessary because neither of SABIO's SBML and Excel export methods provide
        all of the SABIO's content.

        Args:
            tsv (:obj:`str`): TSV-formatted table

        Raises:
            :obj:`ValueError`: if a kinetic law or compartment is not contained in the local sqlite database
        """
        # group properties
        tsv = tsv.split('\n')
        law_properties = {}
        for row in csv.DictReader(tsv, delimiter='\t'):
            entry_id = int(float(row['EntryID']))
            if entry_id not in law_properties:
                law_properties[entry_id] = {
                    'KineticMechanismType': row['KineticMechanismType'],
                    'Tissue': row['Tissue'],
                    'Parameters': [],
                }

            parameter = {}

            # type
            parameter['type'] = row['parameter.type']

            if row['parameter.type'] == 'kcat':
                parameter['type_code'] = 25
            elif row['parameter.type'] == 'Vmax':
                parameter['type_code'] = 186
            elif row['parameter.type'] == 'Km':
                parameter['type_code'] = 27
            elif row['parameter.type'] == 'Ki':
                parameter['type_code'] = 261
            else:
                parameter['type_code'] = None

            # associatated species
            if row['parameter.associatedSpecies'] in ['', '-']:
                parameter['associatedSpecies'] = None
            else:
                parameter['associatedSpecies'] = row['parameter.associatedSpecies']

            # start value
            if row['parameter.startValue'] in ['', '-']:
                parameter['startValue'] = None
            else:
                parameter['startValue'] = float(row['parameter.startValue'])

            # end value
            if row['parameter.endValue'] in ['', '-']:
                parameter['endValue'] = None
            else:
                parameter['endValue'] = float(row['parameter.endValue'])

            # error
            if row['parameter.standardDeviation'] in ['', '-']:
                parameter['standardDeviation'] = None
            else:
                parameter['standardDeviation'] = float(row['parameter.standardDeviation'])

            # units
            if row['parameter.unit'] in ['', '-']:
                parameter['unit'] = None
            else:
                parameter['unit'] = row['parameter.unit']

            law_properties[entry_id]['Parameters'].append(parameter)

        # update properties
        law_q = self.session.query(KineticLaw)
        for id, properties in law_properties.items():
            # get kinetic law
            q = law_q.filter_by(id=id)
            if not q.count():
                raise ValueError('No Kinetic Law with id {}'.format(id))
            law = q.first()

            # mechanism
            if properties['KineticMechanismType'] == 'unknown':
                law.mechanism = None
            else:
                law.mechanism = properties['KineticMechanismType']

            # tissue
            if properties['Tissue'] in ['', '-']:
                law.tissue = None
            else:
                law.tissue = properties['Tissue']

            # parameter
            for param_properties in properties['Parameters']:
                param = self.get_parameter_by_properties(law, param_properties)
                if param is None:
                    if param_properties['type'] != 'concentration':
                        warnings.warn('Unable to find parameter `{}:{}` for law {}'.format(
                            param_properties['type'], param_properties['associatedSpecies'], law.id))
                    continue

                param.observed_value = param_properties['startValue']
                param.observed_error = param_properties['standardDeviation']
                param.observed_units = param_properties['unit']

            # updated
            law.modified = datetime.datetime.utcnow()

    def get_parameter_by_properties(self, kinetic_law, parameter_properties):
        """ Get the parameter of :obj:`kinetic_law` whose attribute values are equal to that of :obj:`parameter_properties`

        Args:
            kinetic_law (:obj:`KineticLaw`): kinetic law to find parameter of
            parameter_properties (:obj:`dict`): properties of parameter to find

        Returns:
            :obj:`Parameter`: parameter with attribute values equal to values of :obj:`parameter_properties`
        """
        if parameter_properties['type'] == 'concentration':
            return None

        # match observed name and compound
        def func(parameter):
            return parameter.observed_type == parameter_properties['type_code'] and \
                ((parameter.compound is None and parameter_properties['associatedSpecies'] is None) or
                 (parameter.compound is not None and parameter.compound.name == parameter_properties['associatedSpecies']))
        parameters = list(filter(func, kinetic_law.parameters))
        if len(parameters) == 1:
            return parameters[0]

        # match observed name
        def func(parameter):
            return parameter.observed_type == parameter_properties['type_code']
        parameters = list(filter(func, kinetic_law.parameters))
        if len(parameters) == 1:
            return parameters[0]

        # match compound
        def func(parameter):
            return (parameter.compound is None and parameter_properties['associatedSpecies'] is None) or \
                (parameter.compound is not None and parameter.compound.name == parameter_properties['associatedSpecies'])
        parameters = list(filter(func, kinetic_law.parameters))
        if len(parameters) == 1:
            return parameters[0]

        # match value
        def func(parameter):
            return parameter.observed_value == parameter_properties['startValue']
        parameters = list(filter(func, kinetic_law.parameters))
        if len(parameters) == 1:
            return parameters[0]

        # return :obj:`None` if there is no matching parameter
        return None

    def load_compounds(self, compounds=None):
        """ Download information from SABIO-RK about all of the compounds stored in the local sqlite copy of SABIO-RK

        Args:
            compounds (:obj:`list` of :obj:`Compound`): list of compounds to download

        Raises:
            :obj:`Error`: if an HTTP request fails
        """
        session = self.requests_session

        if compounds is None:
            compounds = self.session.query(Compound).order_by(Compound.id).all()
        n_compounds = len(compounds)

        for i_compound, c in enumerate(compounds):
            # print status
            if self.verbose and (i_compound % 100 == 0):
                print('  Downloading compound {} of {}'.format(i_compound + 1, n_compounds))

            # download info
            response = session.get(self.ENDPOINT_COMPOUNDS_PAGE, params={'cid': c.id})
            response.raise_for_status()

            # parse info
            doc = bs4.BeautifulSoup(response.text, 'html.parser')
            table = doc.find('table')

            # name
            node = table.find('span', id='commonName')
            if node:
                c.name = node.get_text()

            # synonyms
            c.synonyms = []
            synonym_label_node = table.find('b', text='Synonyms')
            if synonym_label_node:
                for node in list(synonym_label_node.parents)[1].find_all('span'):
                    name = node.get_text()

                    q = self.session.query(Synonym).filter_by(name=name)
                    if q.count():
                        synonym = q[0]
                    else:
                        synonym = Synonym(name=name)
                        self.session.add(synonym)

                    c.synonyms.append(synonym)

            # structure
            c.structures = []

            inchi_label_node = table.find('b', text='InChI')
            if inchi_label_node:
                for node in list(inchi_label_node.parents)[1].find_all('span'):
                    value = node.get_text()
                    q = self.session.query(CompoundStructure).filter_by(format='inchi', value=value)
                    if q.count():
                        compound_structure = q.first()
                    else:
                        compound_structure = CompoundStructure(format='inchi', value=value)
                    c.structures.append(compound_structure)

            smiles_label_node = table.find('b', text='SMILES')
            if smiles_label_node:
                for node in list(smiles_label_node.parents)[1].find_all('span'):
                    value = node.get_text()
                    q = self.session.query(CompoundStructure).filter_by(format='smiles', value=value)
                    if q.count():
                        compound_structure = q.first()
                    else:
                        compound_structure = CompoundStructure(format='smiles', value=value)
                    c.structures.append(compound_structure)

            # cross references
            c.cross_references = []
            for node in table.find_all('a'):

                url = node.get('href')

                id = node.get_text()

                if url.startswith('http://www.genome.jp/dbget-bin/www_bget?cpd:'):
                    namespace = 'kegg.compound'
                elif url.startswith('http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?sid='):
                    namespace = 'pubchem.substance'
                elif url.startswith('http://www.ebi.ac.uk/chebi/searchId.do?chebiId='):
                    namespace = 'chebi'
                    id = 'CHEBI:' + id
                elif url.startswith('https://reactome.org/content/detail/'):
                    namespace = 'reactome'
                elif url.startswith('https://biocyc.org/compound?orgid=META&id='):
                    namespace = 'biocyc'
                elif url.startswith('https://www.metanetx.org/chem_info/'):
                    namespace = 'metanetx.chemical'
                elif url.startswith('http://www.chemspider.com/Chemical-Structure.'):
                    namespace = 'chemspider'
                elif url.startswith('http://sabiork.h-its.org/newSearch?q=sabiocompoundid:'):
                    continue
                else:
                    namespace = html.unescape(node.parent.parent.parent.find_all('td')[0].get_text()).strip()
                    warnings.warn('Compound {} has unkonwn cross reference type to namespace {}'.format(c.id, namespace), 
                        data_source.DataSourceWarning)

                q = self.session.query(Resource).filter_by(namespace=namespace, id=id)
                if q.count():
                    resource = q[0]
                else:
                    resource = Resource(namespace=namespace, id=id)
                    self.session.add(resource)

                c.cross_references.append(resource)

            # udated
            c.modified = datetime.datetime.utcnow()

            if self.commit_intermediate_results and (i_compound % 100 == 99):
                self.session.commit()

    def infer_compound_structures_from_names(self, compounds):
        """ Try to use PubChem to infer the structure of compounds from their names

        Notes: we don't try look up structures from their cross references because SABIO has already gathered
        all structures from their cross references to ChEBI, KEGG, and PubChem

        Args:
            compounds (:obj:`list` of :obj:`Compound`): list of compounds
        """
        resource_query = self.session.query(Resource)
        structure_query = self.session.query(CompoundStructure)

        for i_compound, compound in enumerate(compounds):
            if self.verbose and (i_compound % 100 == 0):
                print('  Trying to infer the structure of compound {} of {}'.format(i_compound + 1, len(compounds)))

            if compound.name == 'Unknown':
                continue

            for i_try in range(self.PUBCHEM_MAX_TRIES):
                try:
                    p_compounds = pubchempy.get_compounds(compound.name, 'name')
                    break
                except pubchempy.PubChemHTTPError:
                    if i_try < self.PUBCHEM_MAX_TRIES - 1:
                        # sleep to avoid getting overloading PubChem server and then try again
                        time.sleep(self.PUBCHEM_TRY_DELAY)
                    else:
                        raise

            for p_compound in p_compounds:
                q = resource_query.filter_by(namespace='pubchem.compound', id=str(p_compound.cid))
                if q.count():
                    resource = q.first()
                else:
                    resource = Resource(namespace='pubchem.compound', id=str(p_compound.cid))
                    self.session.add(resource)

                compound.cross_references.append(resource)

                q = structure_query.filter_by(value=p_compound.inchi, format='inchi')
                if q.count():
                    structure = q.first()
                else:
                    structure = CompoundStructure(value=p_compound.inchi, format='inchi')
                compound.structures.append(structure)

            if self.commit_intermediate_results and (i_compound % 100 == 99):
                self.session.commit()

    def load_missing_enzyme_information_from_html(self, ids):
        """ Loading enzyme subunit information from html

        Args:
            ids (:obj:`list` of :obj:`int`): list of IDs of kinetic laws to download
        """
        db_session = self.session
        req_session = self.requests_session

        kinetic_laws = db_session \
            .query(KineticLaw) \
            .filter(KineticLaw.enzyme_id.isnot(None), KineticLaw.id.in_(ids)) \
            .all()

        for i_kinetic_law, kinetic_law in enumerate(kinetic_laws):
            if self.verbose and (i_kinetic_law % 100 == 0):
                print('  Loading enzyme information for {} of {} kinetic laws'.format(i_kinetic_law + 1, len(kinetic_laws)))

            response = req_session.get(self.ENDPOINT_KINETIC_LAWS_PAGE, params={'kinlawid': kinetic_law.id, 'newinterface': 'true'})
            response.raise_for_status()

            enzyme = kinetic_law.enzyme
            subunits = enzyme.subunits
            for subunit in subunits:
                subunit.coefficient = None

            doc = bs4.BeautifulSoup(response.text, 'html.parser')
            td = doc.find('td', text='Modifier-Catalyst')
            tr = td.parent
            td = tr.find_all('td')[-1]
            inner_html = td.decode_contents(formatter='html').strip() + ' '
            if inner_html == '- ':
                continue

            try:
                subunit_coefficients = self.parse_complex_subunit_structure(inner_html)
            except Exception as error:
                six.reraise(
                    ValueError,
                    ValueError('Subunit structure for kinetic law {} could not be parsed: {}\n\t{}'.format(
                        kinetic_law.id, inner_html, str(error).replace('\n', '\n\t'))),
                    sys.exc_info()[2])

            enzyme.subunits = []
            for subunit_id, coefficient in subunit_coefficients.items():
                q = db_session.query(Resource).filter_by(namespace='uniprot', id=subunit_id)
                if q.count():
                    xref = q.first()
                else:
                    xref = Resource(namespace='uniprot', id=subunit_id)
                subunit = EnzymeSubunit(coefficient=coefficient, cross_references=[xref])
                enzyme.subunits.append(subunit)

    def parse_complex_subunit_structure(self, text):
        """ Parse the subunit structure of complex into a dictionary of subunit coefficients

        Args:
            text (:obj:`str`): subunit structure described with nested parentheses

        Returns:
            :obj:`dict` of :obj:`str`, :obj:`int`: dictionary of subunit coefficients
        """
        # try adding missing parentheses
        n_open = text.count('(')
        n_close = text.count(')')
        if n_open > n_close:
            text += ')' * (n_open - n_close)
        elif n_open < n_close:
            text = '(' * (n_close - n_open) + text

        # for convenenice, add parenthesis at the beginning and end of the string and before and after each subunit
        if text[0] != '(':
            text = '(' + text + ')'
        text = text.replace('<a ', '(<a ').replace('</a>', '</a>)')

        # parse the nested subunit structure
        i = 0
        stack = [{'subunits': {}}]
        while i < len(text):
            if text[i] == '(':
                stack.append({'start': i, 'subunits': {}})
                i += 1
            elif text[i] == ')':
                tmp = stack.pop()
                start = tmp['start']
                end = i

                subunits = tmp['subunits']
                if not subunits:
                    matches = re.findall(r'<a href="http://www\.uniprot\.org/uniprot/(.*?)" target="?_blank"?>.*?</a>', text[start+1:end])
                    for match in matches:
                        subunits[match] = 1

                match = re.match(r'^\)\*(\d+)', text[i:])
                if match:
                    i += len(match.group(0))
                    coefficient = int(float(match.group(1)))
                else:
                    i += 1
                    coefficient = 1

                for id in subunits.keys():
                    if id not in stack[-1]['subunits']:
                        stack[-1]['subunits'][id] = 0
                    stack[-1]['subunits'][id] += subunits[id] * coefficient

            else:
                i += 1

        # check that all subunits were extracted
        matches = re.findall(r'<a href="http://www\.uniprot\.org/uniprot/(.*?)" target="?_blank"?>.*?</a>', text)
        if len(set(matches)) != len(stack[0]['subunits'].keys()):
            raise ValueError('Subunit structure could not be parsed: {}'.format(text))

        return stack[0]['subunits']

    def calc_enzyme_molecular_weights(self, enzymes):
        """ Calculate the molecular weight of each enzyme

        Args:
            enzymes (:obj:`list` of :obj:`Enzyme`): list of enzymes
        """
        letters = Bio.Alphabet.IUPAC.IUPACProtein.letters
        mean_aa_mw = Bio.SeqUtils.molecular_weight(letters, seq_type='protein') / len(letters)

        for i_enzyme, enzyme in enumerate(enzymes):
            if self.verbose and i_enzyme % 100 == 0:
                print('  Calculating molecular weight of enzyme {} of {}'.format(i_enzyme + 1, len(enzymes)))

            enzyme_molecular_weight = 0
            for subunit in enzyme.subunits:
                for xref in subunit.cross_references:
                    if xref.namespace == 'uniprot':
                        response = self.requests_session.get(
                            self.ENDPOINT_DOMAINS['uniprot'] + '/uniprot/?query={}&columns=id,sequence&format=tab'.format(xref.id))
                        response.raise_for_status()
                        seqs = list(csv.DictReader(response.text.split('\n'), delimiter='\t'))
                        if seqs:
                            subunit.sequence = next((seq['Sequence'] for seq in seqs if seq['Entry'] == xref.id), seqs[0]['Sequence'])
                            iupac_seq = re.sub(r'[^' + Bio.Alphabet.IUPAC.IUPACProtein.letters + r']', '', subunit.sequence)
                            subunit.molecular_weight = \
                                + Bio.SeqUtils.molecular_weight(iupac_seq, seq_type='protein') \
                                + (len(subunit.sequence) - len(iupac_seq)) * mean_aa_mw
                            enzyme_molecular_weight += (subunit.coefficient or float('nan')) * subunit.molecular_weight
                        else:
                            subunit.sequence = None
                            subunit.molecular_weight = None
                            enzyme_molecular_weight = float('nan')

            if not enzyme_molecular_weight or math.isnan(enzyme_molecular_weight):
                enzyme.molecular_weight = None
            else:
                enzyme.molecular_weight = enzyme_molecular_weight

    def calc_stats(self):
        """ Calculate statistics about SABIO-RK

        Returns:
            :obj:`list` of :obj:`list` of :obj:`obj`: list of list of statistics
        """
        session = self.session

        units = session \
            .query(Parameter.units, Parameter.observed_units, sqlalchemy.func.count(Parameter._id)) \
            .filter(Parameter.observed_type.in_(Parameter.TYPES.keys())) \
            .group_by(Parameter.units, Parameter.observed_units) \
            .order_by(Parameter.units, Parameter.observed_units) \
            .all()

        row_labels = [
            [None, None, None, None, None]
            + [u[0] for u in units]
            + [None, None],
            ['SBO name', 'SBO ID', 'Entries', 'Entries with values', 'Entries without values']
            + [u[1] for u in units]
            + ['Total usable', 'Percent usable'],
        ]

        columns = []
        for sbo_id, sbo_name in Parameter.TYPES.items():
            column = [sbo_name, sbo_id]
            columns.append(column)

            n_entries = session \
                .query(Parameter) \
                .filter(Parameter.observed_type == sbo_id) \
                .count()
            n_value_entries = session \
                .query(Parameter) \
                .filter(Parameter.observed_type == sbo_id, Parameter.observed_value.isnot(None)) \
                .count()
            n_none_entries = session \
                .query(Parameter) \
                .filter(Parameter.observed_type == sbo_id, Parameter.observed_value == None) \
                .count()
            column.extend([n_entries, n_value_entries, n_none_entries])

            for unit in units:
                column.append(session
                              .query(Parameter)
                              .filter(Parameter.observed_type == sbo_id, Parameter.units == unit[0], Parameter.observed_units == unit[1])
                              .count())

            usable = session \
                .query(Parameter) \
                .filter(Parameter.observed_type == sbo_id, Parameter.value.isnot(None)) \
                .count()
            column.append(usable)
            column.append(usable / n_value_entries * 100 if n_value_entries else float('nan'))

            for i, val in enumerate(column):
                if val == 0:
                    column[i] = None

        return wc_utils.util.list.transpose(row_labels + columns)

    def export_stats(self, stats, filename=None):
        """ Export statistics to an Excel workbook

        Args:
            stats (:obj:`list` of :obj:`list` of :obj:`obj`): list of list of statistics
            filename (:obj:`str`, optional): path to export statistics
        """
        wb = wc_utils.workbook.core.Workbook()
        ws = wb['Stats'] = wc_utils.workbook.core.Worksheet()

        style = wc_utils.workbook.io.WorkbookStyle()
        style['Stats'] = wc_utils.workbook.io.WorksheetStyle(
            head_rows=2, head_columns=2,
            head_row_font_bold=True, head_row_fill_fgcolor='CCCCCC', row_height=15)

        for row in stats:
            ws.append(wc_utils.workbook.core.Row(row))

        if not filename:
            filename = os.path.join(os.path.dirname(self.filename), os.path.splitext(self.filename)[0] + '.summary.xlsx')
        wc_utils.workbook.io.write(filename, wb, style=style)
