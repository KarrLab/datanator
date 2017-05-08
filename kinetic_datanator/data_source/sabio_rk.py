# -*- coding: utf-8 -*-

"""
:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-04
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.util import molecule_util
from requests.packages.urllib3.exceptions import InsecureRequestWarning
from wc_utils.backup import BackupManager
import bs4
import csv
import datetime
import libsbml
import math
import os
import pubchempy
import re
import requests
import requests_cache
import six
import sqlalchemy
import sqlalchemy.ext.declarative
import sqlalchemy.orm
import sys

_, _, NAME = __name__.rpartition('.')

DEFAULT_CACHE_DIRNAME = os.path.join(os.path.dirname(__file__), 'cache')
# :obj:`str`: default path for the sqlite database

SqlalchemyBase = sqlalchemy.ext.declarative.declarative_base()
# :obj:`Base`: SQL Alchemy base class for class definitions


def get_engine(filename=None):
    """ Get an engine for the sqlite database

    Args:
        filename (:obj:`str`, optional): path to sqlite file

    Returns:
        :obj:`sqlalchemy.engine.Engine`: database engine
    """
    if not filename:
        filename = os.path.join(DEFAULT_CACHE_DIRNAME, NAME + '.sqlite')
    if not os.path.isdir(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename))
    return sqlalchemy.create_engine('sqlite:///' + filename)


def get_session(engine=None, auto_download=True, auto_update=False, force_download=False, force_update=False,
                requests_cache_filename=None, max_entries=float('inf'), arcname=None):
    """ Get a session for the sqlite database

    Args:
        engine (:obj:`sqlalchemy.engine.Engine`, optional): database engine
        auto_download (:obj:`bool`, optional): if :obj:`True` and there is no local sqlite database, download an archive of
            the database from the Karr Lab server and then update the local database from SABIO
        auto_update (:obj:`bool`, optional): if :obj:`True` and there is no local sqlite database, update the local database
            from SABIO
        force_download (:obj:`bool`, optional): if :obj:`True`, download an archive of the database from the Karr Lab server
            Note: this setting will be overriden to :obj:`True` if there is no local sqlite data and :obj:`auto_download` is
            :obj:`True`
        force_update (:obj:`bool`, optional): if :obj:`True`, update the local database from SABIO
            Note: this setting will be overriden to :obj:`True` if there is no local sqlite database and :obj:`auto_download`
            or :obj:`auto_update` is :obj:`True`
        requests_cache_filename (:obj:`str`, optional): filename of the sqlite database to cache SABIO-RK HTTP requests
        max_entries (:obj:`int`, optional): maximum number of laws to download
        arcname (:obj:`str`, optional): filename to download sqlite database from remote server

    Returns:
        :obj:`sqlalchemy.orm.session.Session`: database session
    """
    if not engine:
        engine = get_engine()
    filename = str(engine.url).replace('sqlite:///', '')

    if not os.path.isfile(filename):
        if auto_download:
            force_download = True
        elif auto_update:
            force_update = True

    if force_download and os.getenv('CODE_SERVER_TOKEN'):
        download(filename=filename, arcname=arcname)

    session = sqlalchemy.orm.sessionmaker(bind=engine)()

    if force_update:
        Downloader(session=session, requests_cache_filename=requests_cache_filename, max_entries=max_entries).download()

    return session


def backup(filename=None, arcname=None):
    """ Backup the local sqlite database to the Karr Lab server

    Args:
        filename (:obj:`str`, optional): path to sqlite file
        arcname (:obj:`str`, optional): filename to save sqlite database on remote server
    """
    if not filename:
        filename = os.path.join(DEFAULT_CACHE_DIRNAME, NAME + '.sqlite')
    if not arcname:
        arcname = NAME + '.sqlite'
    BackupManager(filename, arcname=arcname) \
        .create() \
        .upload() \
        .cleanup()


def download(filename=None, arcname=None):
    """ Download the local sqlite database from the Karr Lab server

    Args:
        filename (:obj:`str`, optional): path to sqlite file
        arcname (:obj:`str`, optional): filename to download sqlite database from remote server
    """
    if not filename:
        filename = os.path.join(DEFAULT_CACHE_DIRNAME, NAME + '.sqlite')
    if not arcname:
        arcname = NAME + '.sqlite'
    BackupManager(filename, arcname=arcname) \
        .download() \
        .extract() \
        .cleanup()

entry_synonym = sqlalchemy.Table(
    'entry_synonym', SqlalchemyBase.metadata,
    sqlalchemy.Column('entry__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('entry._id')),
    sqlalchemy.Column('synonym__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('synonym._id')),
)
# :obj:`sqlalchemy.Table`: Entry:Synonym many-to-many association table

entry_resource = sqlalchemy.Table(
    'entry_resource', SqlalchemyBase.metadata,
    sqlalchemy.Column('entry__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('entry._id')),
    sqlalchemy.Column('resource__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('resource._id')),
)
# :obj:`sqlalchemy.Table`: Entry:Resource many-to-many association table

compound_compound_structure = sqlalchemy.Table(
    'compound_compound_structure', SqlalchemyBase.metadata,
    sqlalchemy.Column('compound__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('compound._id')),
    sqlalchemy.Column('compound_structure__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('compound_structure._id')),
)
# :obj:`sqlalchemy.Table`: Compound:CompoundStructure many-to-many association table


kinetic_law_resource = sqlalchemy.Table(
    'kinetic_law_resource', SqlalchemyBase.metadata,
    sqlalchemy.Column('kinetic_law__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('kinetic_law._id')),
    sqlalchemy.Column('resource__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('resource._id')),
)
# :obj:`sqlalchemy.Table`: KineticLaw:Resource many-to-many association table


class Synonym(SqlalchemyBase):
    """ Represents a synonym to a SABIO-RK entry

    Attributes:
        name (:obj:`str`): name of the synonym
        entries (:obj:`list` of :obj:`Entry`): list of entries with the synonym
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    name = sqlalchemy.Column(sqlalchemy.String(), index=True)
    entries = sqlalchemy.orm.relationship('Entry', secondary=entry_synonym, back_populates='synonyms')

    __tablename__ = 'synonym'


class Resource(SqlalchemyBase):
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
    entries = sqlalchemy.orm.relationship('Entry', secondary=entry_resource, back_populates='cross_references')
    kinetic_laws = sqlalchemy.orm.relationship('KineticLaw', secondary=kinetic_law_resource, back_populates='references')

    sqlalchemy.schema.UniqueConstraint(namespace, id)

    __tablename__ = 'resource'


class Entry(SqlalchemyBase):
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
    synonyms = sqlalchemy.orm.relationship('Synonym', secondary=entry_synonym, back_populates='entries')
    cross_references = sqlalchemy.orm.relationship('Resource', secondary=entry_resource, back_populates='entries')
    created = sqlalchemy.Column(sqlalchemy.DateTime, default=datetime.datetime.utcnow())
    modified = sqlalchemy.Column(sqlalchemy.DateTime, onupdate=datetime.datetime.utcnow())

    sqlalchemy.schema.UniqueConstraint(id, _type)

    __tablename__ = 'entry'
    __mapper_args__ = {'polymorphic_on': _type}


class ReactionParticipant(SqlalchemyBase):
    """ Represents a participant in a SABIO-RK reaction

    Attributes:
        compound (:obj:`Compound`): compound
        compartment (:obj:`Compartment`): compartment
        coefficient (:obj:`float`): coefficient
        type (:obj:`str`): type
        reactant_reaction (:obj:`Reaction`): reaction in which the participant appears as a reactant
        product_reaction (:obj:`Reaction`): reaction in which the participant appears as a product
        modifier_reaction (:obj:`Reaction`): reaction in which the participant appears as a modifier
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    compound_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('compound._id'))
    compound = sqlalchemy.orm.relationship('Compound', back_populates='reaction_participants', foreign_keys=[compound_id])
    compartment_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('compartment._id'))
    compartment = sqlalchemy.orm.relationship('Compartment', back_populates='reaction_participants', foreign_keys=[compartment_id])
    coefficient = sqlalchemy.Column(sqlalchemy.Float())
    type = sqlalchemy.Column(sqlalchemy.String())
    reactant_reaction_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('reaction._id'))
    reactant_reaction = sqlalchemy.orm.relationship(
        'Reaction', uselist=False, back_populates='reactants', foreign_keys=[reactant_reaction_id])
    product_reaction_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('reaction._id'))
    product_reaction = sqlalchemy.orm.relationship(
        'Reaction', uselist=False, back_populates='products', foreign_keys=[product_reaction_id])
    kinetic_law_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('kinetic_law._id'))
    kinetic_law = sqlalchemy.orm.relationship(
        'KineticLaw', uselist=False, back_populates='modifiers', foreign_keys=[kinetic_law_id])

    __tablename__ = 'reaction_participant'


class Parameter(Entry):
    """ Represents a parameter in the SABIO-RK database

    Attributes:
        kinetic_law (:obj:`KineticLaw`): kinetic law
        type (:obj:`str`): type (kcat, Km, etc.)
        compound (:obj:`Compound`): compound
        value (:obj:`float`): value
        units (:obj:`str`): units
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('entry._id'), primary_key=True)

    kinetic_law_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('kinetic_law._id'))
    kinetic_law = sqlalchemy.orm.relationship('KineticLaw', uselist=False, back_populates='parameters', foreign_keys=[kinetic_law_id])
    type = sqlalchemy.Column(sqlalchemy.String())
    compound_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('compound._id'))
    compound = sqlalchemy.orm.relationship('Compound', uselist=False, back_populates='parameters', foreign_keys=[compound_id])
    compartment_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('compartment._id'))
    compartment = sqlalchemy.orm.relationship('Compartment', uselist=False, back_populates='parameters', foreign_keys=[compartment_id])
    value = sqlalchemy.Column(sqlalchemy.Float())
    units = sqlalchemy.Column(sqlalchemy.String())

    __tablename__ = 'parameter'
    __mapper_args__ = {'polymorphic_identity': 'parameter'}


class KineticLaw(Entry):
    """ Represents a kinetic law in the SABIO-RK database

    Attributes:
        reaction (:obj:`Reaction`): reaction
        enzyme (:obj:`Enzyme`): enzyme
        enzyme_compartment (:obj:`Compartment): compartment
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

    reaction_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('reaction._id'))
    reaction = sqlalchemy.orm.relationship('Reaction', uselist=False, back_populates='kinetic_law', foreign_keys=[reaction_id])
    enzyme_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('enzyme._id'))
    enzyme = sqlalchemy.orm.relationship('Enzyme', uselist=False, back_populates='kinetic_laws', foreign_keys=[enzyme_id])
    enzyme_compartment_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('compartment._id'))
    enzyme_compartment = sqlalchemy.orm.relationship(
        'Compartment', uselist=False, back_populates='kinetic_laws', foreign_keys=[enzyme_compartment_id])
    enzyme_type = sqlalchemy.Column(sqlalchemy.String())
    tissue = sqlalchemy.Column(sqlalchemy.String())
    mechanism = sqlalchemy.Column(sqlalchemy.String())
    equation = sqlalchemy.Column(sqlalchemy.Text())
    parameters = sqlalchemy.orm.relationship('Parameter', back_populates='kinetic_law', foreign_keys=[Parameter.kinetic_law_id])
    modifiers = sqlalchemy.orm.relationship('ReactionParticipant', back_populates='kinetic_law',
                                            foreign_keys=[ReactionParticipant.kinetic_law_id])
    taxon = sqlalchemy.Column(sqlalchemy.Integer())
    taxon_wildtype = sqlalchemy.Column(sqlalchemy.Boolean())
    taxon_variant = sqlalchemy.Column(sqlalchemy.UnicodeText())
    temperature = sqlalchemy.Column(sqlalchemy.Float())
    ph = sqlalchemy.Column(sqlalchemy.Float())
    media = sqlalchemy.Column(sqlalchemy.UnicodeText())
    references = sqlalchemy.orm.relationship('Resource', secondary=kinetic_law_resource, back_populates='kinetic_laws')

    __tablename__ = 'kinetic_law'
    __mapper_args__ = {'polymorphic_identity': 'kinetic_law'}


class Reaction(Entry):
    """ Represents a reaction in the SABIO-RK database

    Attributes:
        reactants (:obj:`list` of :obj:`ReactionParticipant`): list of reactants
        products (:obj:`list` of :obj:`ReactionParticipant`): list of products
        kinetic_law (:obj:`KineticLaw`): kinetic law
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('entry._id'), primary_key=True)

    reactants = sqlalchemy.orm.relationship('ReactionParticipant', back_populates='reactant_reaction',
                                            foreign_keys=[ReactionParticipant.reactant_reaction_id])
    products = sqlalchemy.orm.relationship('ReactionParticipant', back_populates='product_reaction',
                                           foreign_keys=[ReactionParticipant.product_reaction_id])
    kinetic_law = sqlalchemy.orm.relationship('KineticLaw', uselist=False, back_populates='reaction',
                                              foreign_keys=[KineticLaw.reaction_id])

    __tablename__ = 'reaction'
    __mapper_args__ = {'polymorphic_identity': 'reaction'}


class CompoundStructure(SqlalchemyBase):
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

    compounds = sqlalchemy.orm.relationship('Compound', secondary=compound_compound_structure, back_populates='structures')
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
    structures = sqlalchemy.orm.relationship('CompoundStructure', secondary=compound_compound_structure, back_populates='compounds')
    reaction_participants = sqlalchemy.orm.relationship('ReactionParticipant', back_populates='compound')
    parameters = sqlalchemy.orm.relationship('Parameter', back_populates='compound', foreign_keys=[Parameter.compound_id])

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
        kinetic_laws (:obj:`list` of :obj:`KineticLaw`): list of kinetic laws
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('entry._id'), primary_key=True)

    kinetic_laws = sqlalchemy.orm.relationship('KineticLaw', back_populates='enzyme', foreign_keys=[KineticLaw.enzyme_id])

    __tablename__ = 'enzyme'
    __mapper_args__ = {'polymorphic_identity': 'enzyme'}


class Compartment(Entry):
    """ Represents a compartment in the SABIO-RK database

    Attributes:
        kinetic_laws (:obj:`list` of :obj:`KineticLaw`): list of kinetic laws
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('entry._id'), primary_key=True)

    reaction_participants = sqlalchemy.orm.relationship('ReactionParticipant', back_populates='compartment')
    kinetic_laws = sqlalchemy.orm.relationship('KineticLaw', back_populates='enzyme_compartment',
                                               foreign_keys=[KineticLaw.enzyme_compartment_id])
    parameters = sqlalchemy.orm.relationship('Parameter', back_populates='compartment',
                                             foreign_keys=[Parameter.compartment_id])

    __tablename__ = 'compartment'
    __mapper_args__ = {'polymorphic_identity': 'compartment'}


class Downloader(object):
    """ Downloads the content of the SABIO-RK database into a local sqlite database

    Attributes:
        session (:obj:`sqlalchemy.orm.session.Session`): database session
        index_batch_size (:obj:`int`): size of batches to download the IDs of the kinetic laws
        webservice_batch_size (:obj:`int`): size of batches to download kinetic information from the SABIO webservice
        excel_batch_size (:obj:`int`): size of batches to download kinetic information from the SABIO Excel download service
        compound_batch_size (:obj:`int`): size of batches to download compound information
        requests_cache_filename (:obj:`str`): filename of the sqlite database to cache SABIO-RK HTTP requests
        max_entries (:obj:`int`): maximum number of laws to download
        verbose (:obj:`bool`): if :obj:`True`, print status information to the standard output

        ENDPOINT_KINETIC_LAWS_SEARCH (:obj:`str`): URL to obtain a list of the ids of all of the kinetic laws in SABIO-Rk
        ENDPOINT_WEBSERVICE (:obj:`str`): URL for the SABIO-RK webservice
        ENDPOINT_EXCEL_EXPORT (:obj:`str`): URL to download kinetic data as a table in TSV format
        ENDPOINT_COMPOUNDS_PAGE (:obj:`str`): URL to download information about a SABIO-RK compound
        DEFAULT_INDEX_BATCH_SIZE (:obj:`int`): size of batches to download the IDs of the kinetic laws
        DEFAULT_WEBSERVICE_BATCH_SIZE (:obj:`int`): default size of batches to download kinetic information from the SABIO webservice
        DEFAULT_EXCEL_BATCH_SIZE (:obj:`int`): default size of batches to download kinetic information from the SABIO
            Excel download service
        DEFAULT_COMPOUND_BATCH_SIZE (:obj:`int`): size of batches to download compound information
        SKIP_KINETIC_LAW_IDS (:obj:`tuple` of :obj:`int`): IDs of kinetic laws that should be skipped (because they cannot contained
            errors and can't be downloaded from SABIO)
    """

    ENDPOINT_KINETIC_LAWS_SEARCH = 'http://sabio.villa-bosch.de/newSearch/search'
    ENDPOINT_WEBSERVICE = 'http://sabiork.h-its.org/sabioRestWebServices/kineticLaws'
    ENDPOINT_EXCEL_EXPORT = 'http://sabio.villa-bosch.de/entry/exportToExcelCustomizable'
    ENDPOINT_COMPOUNDS_PAGE = 'http://sabio.villa-bosch.de/compdetails.jsp'
    DEFAULT_INDEX_BATCH_SIZE = 1000
    DEFAULT_WEBSERVICE_BATCH_SIZE = 250
    DEFAULT_EXCEL_BATCH_SIZE = 100
    DEFAULT_COMPOUND_BATCH_SIZE = 100
    SKIP_KINETIC_LAW_IDS = (51286,)

    def __init__(self, session,
                 index_batch_size=DEFAULT_INDEX_BATCH_SIZE,
                 webservice_batch_size=DEFAULT_WEBSERVICE_BATCH_SIZE,
                 excel_batch_size=DEFAULT_EXCEL_BATCH_SIZE,
                 compound_batch_size=DEFAULT_COMPOUND_BATCH_SIZE,
                 requests_cache_filename=None, max_entries=float('inf'), verbose=False):
        """
        Args:
            session (:obj:`sqlalchemy.orm.session.Session`): database session
            index_batch_size (:obj:`int`, optional): size of batches to download the IDs of the kinetic laws
            webservice_batch_size (:obj:`int`, optional): size of batches to download kinetic information from the SABIO webservice
            compound_batch_size (:obj:`int`, optional): size of batches to download compound information
            excel_batch_size (:obj:`int`, optional): size of batches to download kinetic information from the SABIO Excel download service
            requests_cache_filename (:obj:`str`, optional): filename of the sqlite database to cache SABIO-RK HTTP requests
            max_entries (:obj:`int`, optional): maximum number of laws to download
            verbose (:obj:`bool`, optional): if :obj:`True`, print status information to the standard output
        """
        self.session = session

        self.index_batch_size = index_batch_size
        self.webservice_batch_size = webservice_batch_size
        self.excel_batch_size = excel_batch_size
        self.compound_batch_size = compound_batch_size

        if requests_cache_filename is None:
            requests_cache_filename = os.path.join(DEFAULT_CACHE_DIRNAME, NAME + '.requests.py{}.sqlite'.format(sys.version_info[0]))
        elif not requests_cache_filename.endswith('.sqlite'):
            raise ValueError('Request cache filename must have the extension .sqlite')
        if not os.path.isdir(os.path.dirname(requests_cache_filename)):
            os.makedirs(os.path.dirname(requests_cache_filename))
        self.requests_cache_filename = requests_cache_filename

        self.max_entries = max_entries
        self.verbose = verbose

        # initialize cache, database
        self.init_requests_cache()
        self.init_database()

    def init_requests_cache(self):
        """ Setup the cache for SABIO-RK HTTP requests """
        pass

    def clear_requests_cache(self):
        """ Clear the cahce for SABIO-RK HTTP requests """
        name, _, _ = self.requests_cache_filename.partition('.')
        session = requests_cache.core.CachedSession(name, backend='sqlite', expire_after=None)
        session.cache.clear()

    def init_database(self):
        """ Initialize the local sqlite database to store the content of SABIO-RK """
        engine = self.session.get_bind()
        if not os.path.isfile(str(engine.url).replace('sqlite:///', '')):
            SqlalchemyBase.metadata.create_all(engine)

    def clear_database(self):
        """ Clear the content of the local sqlite database to store the content of SABIO-RK """
        engine = self.session.get_bind()
        SqlalchemyBase.metadata.drop_all(engine)
        SqlalchemyBase.metadata.create_all(engine)

    def download(self, update_database=False, update_requests=False):
        """ Download the content of SABIO-RK and store it to a local sqlite database.
        Optionally, clear any existing database and/or HTTP request cache.

        Args:
            update_database (:obj:`bool`): if :obj:`True`, clear the existing database
            update_requests (:obj:`bool`): if :obj:`True`, clear HTTP request cache
        """
        if update_database:
            self.clear_database()

        if update_requests:
            self.clear_requests_cache()

        # disable requests warning
        requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

        # determine ids of kinetic laws
        if self.verbose:
            print('Downloading the IDs of the kinetic laws ...')

        ids = self.download_kinetic_law_ids()

        if self.verbose:
            print('  Downloaded {} IDs'.format(len(ids)))

        # remove bad IDs
        ids = filter(lambda id: id not in self.SKIP_KINETIC_LAW_IDS, ids)

        # download kinetic laws
        new_ids = list(set(ids).difference(set(l.id for l in self.session.query(KineticLaw).all())))
        new_ids.sort()

        if self.verbose:
            print('Downloading {} kinetic laws ...'.format(len(new_ids)))

        self.download_kinetic_laws(new_ids)

        if self.verbose:
            print('  done')

        # download compounds
        compounds = self.session.query(Compound).filter(~Compound.structures.any()).order_by(Compound.id).all()

        if self.verbose:
            print('Downloading {} compounds ...'.format(len(compounds)))

        self.download_compounds(compounds)

        if self.verbose:
            print('  done')

        # infer structures for compounds with no provided structure
        compounds = self.session.query(Compound).filter(~Compound.structures.any()).all()

        if self.verbose:
            print('Inferring structures for {} compounds ...'.format(len(compounds)))

        self.infer_compound_structures_from_names(compounds)

        if self.verbose:
            print('  done')

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

        # commit changes
        self.session.commit()

    def download_kinetic_law_ids(self):
        """ Download the IDs of all of the kinetic laws stored in SABIO-RK

        Raises:
            :obj:`Error`: if an HTTP request fails or the expected number of kinetic laws is not returned
        """
        batch_size = self.index_batch_size

        # create session
        session = requests.Session()
        response = session.post(self.ENDPOINT_KINETIC_LAWS_SEARCH, params={
            'q': '',
            'ontologySearch': 'false',
            'wildtype': 'true',
            'mutant': 'true',
            'recombinant': 'false',
            'kineticData': 'false',
            'transportReaction': 'false',
            'phValues': '0 - 14',
            'temperatures': '-10 C° - 115 C°',
            'directSubmission': 'true',
            'journal': 'true',
            'biomodel': 'true',
            'date': 'false',
            'entryDate': '14/10/2008',
            'ipAddress': '127.0.1.1',
            'ipAddress2': '127.0.0.1',
            'remoteHost': '127.0.0.1',
            'view': 'entry',
        })
        response.raise_for_status()

        # get IDs of kinetic laws
        ids = []
        i_batch = 0
        while True:
            # get page
            response = session.post(self.ENDPOINT_KINETIC_LAWS_SEARCH, params={
                'view': 'entry',
                'offset': i_batch * batch_size,
                'max': batch_size,
            })
            response.raise_for_status()
            if response.text.find('Sorry, we found no results for your query...') != -1:
                raise ValueError('Error creating session with SABIO-RK')

            # parse page
            doc = bs4.BeautifulSoup(response.text, 'html.parser')

            # get number of laws
            n_laws = int(float(doc.find('div', id='numberofKinLaw').find('span').get_text()))

            # print status
            if self.verbose:
                print('  Downloaded ids {}-{} of {}'.format(
                    i_batch * batch_size + 1,
                    min(n_laws, (i_batch + 1) * batch_size),
                    min(self.max_entries, n_laws),
                ))

            # get ids of kinetic laws
            more_ids = [int(float(node.get('kinlawentryid'))) for node in doc.find_all('img', attrs={'name': 'kinlawEntryId'})]
            if len(more_ids) != min(batch_size, min(n_laws, (i_batch + 1) * batch_size) - i_batch * batch_size):
                raise ValueError('Batch {} failed to download the expected number of kinetic law ids'.format(i_batch))
            ids.extend(more_ids)

            # increment batch number
            i_batch += 1

            # stop if all IDs have been collected
            if len(ids) >= min(self.max_entries, n_laws):
                if len(ids) > self.max_entries:
                    ids = ids[0:self.max_entries]
                break

        # return IDs
        ids.sort()
        return ids

    def download_kinetic_laws(self, ids):
        """ Download kinetic laws from SABIO-RK

        Args:
            ids (:obj:`list` of :obj:`int`): list of IDs of kinetic laws to download

        Raises:
            :obj:`Error`: if an HTTP request fails
        """
        # todo: scrape strain, recombinant, experiment type information from web pages
        name, _, _ = self.requests_cache_filename.partition('.')
        session = requests_cache.core.CachedSession(name, backend='sqlite', expire_after=None)

        batch_size = self.webservice_batch_size
        for i_batch in range(int(math.ceil(float(len(ids)) / batch_size))):
            if self.verbose:
                print('  Downloading kinetic laws {}-{} of {} in SBML format'.format(
                    i_batch * batch_size + 1,
                    min(len(ids), (i_batch + 1) * batch_size),
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

            self.create_kinetic_laws_from_sbml(response.content if six.PY2 else response.text)

        # download information from custom Excel export
        batch_size = self.excel_batch_size
        for i_batch in range(int(math.ceil(float(len(ids)) / batch_size))):
            if self.verbose:
                print('  Downloading kinetic laws {}-{} of {} in Excel format'.format(
                    i_batch * batch_size + 1,
                    min(len(ids), (i_batch + 1) * batch_size),
                    len(ids)))

            response = session.get(self.ENDPOINT_EXCEL_EXPORT, params={
                'entryIDs[]': batch_ids,
                'fields[]': [
                    'EntryID',
                    'KineticMechanismType',
                    'Tissue',
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

            self.update_kinetic_laws_from_tsv(response.text)

    def create_kinetic_laws_from_sbml(self, sbml):
        """ Add kinetic laws defined in an SBML file to the local sqlite database

        Args:
            sbml (:obj:`str`): SBML representation of one or more kinetic laws
        """
        reader = libsbml.SBMLReader()
        doc = reader.readSBMLFromString(sbml)
        model = doc.getModel()

        functions = {}
        functions_sbml = model.getListOfFunctionDefinitions()
        for i_function in range(functions_sbml.size()):
            function_sbml = functions_sbml.get(i_function)
            eq = libsbml.formulaToString(function_sbml.getMath())
            _, _, eq = eq.rpartition(', ')
            eq = eq[0:-1]
            if eq == 'NaN':
                eq = ''
            functions[function_sbml.getId()] = eq or None

        units = {}
        units_sbml = model.getListOfUnitDefinitions()
        for i_unit in range(units_sbml.size()):
            unit_sbml = units_sbml.get(i_unit)
            units[unit_sbml.getId()] = unit_sbml.getName()

        # compartments
        compartments_sbml = model.getListOfCompartments()
        for i_compartment in range(compartments_sbml.size()):
            compartment_sbml = compartments_sbml.get(i_compartment)
            self.create_compartment_from_sbml(compartment_sbml)

        # species
        specie_properties = {}
        species_sbml = model.getListOfSpecies()
        for i_specie in range(species_sbml.size()):
            specie_sbml = species_sbml.get(i_specie)
            _, properties = self.create_specie_from_sbml(specie_sbml)
            specie_properties[specie_sbml.getId()] = properties

        # reactions
        reactions_sbml = model.getListOfReactions()
        for i_reaction in range(reactions_sbml.size()):
            reaction_sbml = reactions_sbml.get(i_reaction)
            self.create_reaction_from_sbml(reaction_sbml, specie_properties)

        # kinetic laws
        reactions_sbml = model.getListOfReactions()
        for i_reaction in range(reactions_sbml.size()):
            reaction_sbml = reactions_sbml.get(i_reaction)
            self.create_kinetic_law_from_sbml(reaction_sbml, specie_properties, functions, units)

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
        if specie.name is not None and specie.name != name:
            specie._is_name_ambiguous = True
        specie.name = name
        specie.cross_references = self.create_cross_references_from_sbml(sbml)

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
        match = re.match('^(.*?)\(Enzyme\) (wildtype|mutant),?(.*?)$', sbml, re.IGNORECASE)
        if match:
            name = match.group(1)
            is_wildtype = match.group(2).lower() == 'wildtype'
            variant = match.group(3).strip()
            return (name, is_wildtype, variant)

        match = re.match('^Enzyme (wildtype|mutant),?( (.*?))*$', sbml, re.IGNORECASE)
        if match:
            if match.group(3):
                name = match.group(3).strip()
            else:
                name = None
            is_wildtype = match.group(1).lower() == 'wildtype'
            variant = None
            return (name, is_wildtype, variant)

        match = re.match('^Enzyme - *$', sbml, re.IGNORECASE)
        if match:
            name = None
            is_wildtype = True
            variant = None
            return (name, is_wildtype, variant)

        raise ValueError('Cannot parse enzyme name: {}'.format(sbml))

    def create_reaction_from_sbml(self, sbml, specie_properties):
        """ Add a reaction to the local sqlite database

        Args:
            sbml (:obj:`libsbml.Reaction`): SBML-representation of a reaction
            specie_properties (:obj:`dict`): additional properties of the compounds/enzymes

                * `is_wildtype` (:obj:`bool`): indicates if the enzyme is wildtype or mutant
                * `variant` (:obj:`str`): description of the variant of the eznyme
                * `modifier_type` (:obj:`str`): type of the enzyme (e.g. Modifier-Catalyst)

        Returns:
            :obj:`Reaction`: reaction
        """
        """ annotation """
        x_refs = self.create_cross_references_from_sbml(sbml)

        # id
        id = next(int(float(x_ref.id)) for x_ref in x_refs if x_ref.namespace == 'sabiork.reaction')
        query = self.session.query(Reaction).filter_by(id=id)
        if query.count():
            reaction = query.first()
        else:
            reaction = Reaction(id=id)
            self.session.add(reaction)

        # cross references
        x_refs = list(filter(lambda x_ref: x_ref.namespace not in ('sabiork.reaction', 'taxonomy'), x_refs))
        reaction.cross_references = x_refs

        """ participants """
        reaction.reactants[:] = []
        reactants = sbml.getListOfReactants()
        for i_part in range(reactants.size()):
            part_sbml = reactants.get(i_part)
            compound, compartment = self.get_specie_reference_from_sbml(part_sbml.getSpecies())
            part = ReactionParticipant(
                compound=compound,
                compartment=compartment,
                coefficient=part_sbml.getStoichiometry())
            self.session.add(part)
            reaction.reactants.append(part)

        reaction.products[:] = []
        products = sbml.getListOfProducts()
        for i_part in range(products.size()):
            part_sbml = products.get(i_part)
            compound, compartment = self.get_specie_reference_from_sbml(part_sbml.getSpecies())
            part = ReactionParticipant(
                compound=compound,
                compartment=compartment,
                coefficient=part_sbml.getStoichiometry())
            self.session.add(part)
            reaction.products.append(part)

        """ updated """
        reaction.modified = datetime.datetime.utcnow()

        """ return reaction """
        return reaction

    def create_kinetic_law_from_sbml(self, sbml, specie_properties, functions, units):
        """ Add a kinetic law to the local sqlite database

        Args:
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
        id = next(int(float(x_ref.id)) for x_ref in x_refs if x_ref.namespace == 'sabiork.kineticrecord')
        query = self.session.query(KineticLaw).filter_by(id=id)
        if query.count():
            kinetic_law = query.first()
        else:
            kinetic_law = KineticLaw(id=id)
            self.session.add(kinetic_law)

        # reaction id
        reaction_id = next(int(float(x_ref.id)) for x_ref in reaction_x_refs if x_ref.namespace == 'sabiork.reaction')
        kinetic_law.reaction = self.session.query(Reaction).filter_by(id=reaction_id).first()

        # rate_law
        kinetic_law.equation = functions[law.getMetaId()[5:]]

        # parameters
        kinetic_law.parameters = []
        params = law.getListOfLocalParameters()
        for i_param in range(params.size()):
            param = params.get(i_param)

            match = re.match('^(.*?)_((SPC|ENZ)_([0-9]+)_(.*?))$', param.getId(), re.IGNORECASE)
            if match:
                type = match.group(1)
                compound, compartment = self.get_specie_reference_from_sbml(match.group(2))
            else:
                type = param.getId()
                compound = None
                compartment = None

            if param.getUnits():
                if param.getUnits() in units:
                    param_units = units[param.getUnits()]
                else:
                    param_units = param.getUnits()
            else:
                param_units = None
            parameter = Parameter(
                type=type.replace('div', '/'),
                compound=compound,
                compartment=compartment,
                value=param.getValue(),
                units=param_units,
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
            if temperature_units != '°C':
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
        x_refs = list(filter(lambda x_ref: x_ref.namespace != 'sabiork.kineticrecord', x_refs))
        kinetic_law.references = x_refs

        """ updated """
        kinetic_law.modified = datetime.datetime.utcnow()

        return kinetic_law

    def create_cross_references_from_sbml(self, sbml):
        """ Add cross references to the local sqlite database for an SBML object

        Args:
            sbml (:obj:`libsbml.SBase`): object in an SBML documentation

        Returns:
            :obj:`list` of `Resource`: list of resources
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
                    _, _, _, namespace, id = val.split('/')
                    resources = [(namespace, id)]

                for namespace, id in resources:
                    query = self.session.query(Resource).filter_by(namespace=namespace, id=id)
                    if query.count():
                        resource = query.first()
                    else:
                        resource = Resource(namespace=namespace, id=id)
                        self.session.add(resource)

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
            :obj:`ValueError`: if the species is not a compound or enzyme
        """
        tmp = specie_id.split('_')
        type = tmp[0]
        specie_id = int(float(tmp[1]))
        compartment_name = '_'.join(tmp[2:])

        if type == 'SPC':
            specie = self.session.query(Compound).filter_by(id=specie_id).first()
        elif type == 'ENZ':
            specie = self.session.query(Enzyme).filter_by(id=specie_id).first()
        else:
            raise ValueError('Unsupported species type {}'.format(type))

        if compartment_name != 'Cell':
            compartment = self.session.query(Compartment).filter_by(name=compartment_name).first()
        else:
            compartment = None

        return (specie, compartment)

    def update_kinetic_laws_from_tsv(self, tsv):
        """ Update the properties of kinetic laws in the local sqlite database based on content downloaded
        from SABIO in TSV format.

        Note: this method is necessary because neither of SABIO's SBML and Excel export methods provide
        all of the SABIO's content.

        Args:
            tsv (:obj:`str`): TSV-formatted table

        Raises:
            :obj:`ValueError`: if a kinetic law or compartment is not contained in the local sqlite database
        """
        law_q = self.session.query(KineticLaw)

        tsv = tsv.split('\n')
        for row in csv.DictReader(tsv, delimiter='\t'):
            # get kinetic law
            l = law_q.filter_by(id=int(float(row['EntryID'])))
            if l.count() == 0:
                raise ValueError('No Kinetic Law with id {}'.format(row['EntryID']))
            l = l.first()

            # mechanism
            if row['KineticMechanismType'] == 'unknown':
                l.mechanism = None
            else:
                l.mechanism = row['KineticMechanismType']

            # tissue
            if row['Tissue'] == '-':
                l.tissue = None
            else:
                l.tissue = row['Tissue']

            # updated
            l.modified = datetime.datetime.utcnow()

    def download_compounds(self, compounds=None):
        """ Download information from SABIO-RK about all of the compounds stored in the local sqlite copy of SABIO-RK

        Args:
            compounds (:obj:`list` of :obj:`Compound`): list of compounds to download

        Raises:
            :obj:`Error`: if an HTTP request fails
        """
        name, _, _ = self.requests_cache_filename.partition('.')
        session = requests_cache.core.CachedSession(name, backend='sqlite', expire_after=None)

        batch_size = self.compound_batch_size

        if compounds is None:
            compounds = self.session.query(Compound).order_by(Compound.id).all()
        n_compounds = len(compounds)

        for i_compound, c in enumerate(compounds):
            # print status
            if self.verbose and (i_compound % batch_size == 0):
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
                else:
                    ValueError('Compound {} has unkonwn cross reference type to namespace {}'.format(c.id, url))

                q = self.session.query(Resource).filter_by(namespace=namespace, id=id)
                if q.count():
                    resource = q[0]
                else:
                    resource = Resource(namespace=namespace, id=id)
                    self.session.add(resource)

                c.cross_references.append(resource)

            # udated
            c.modified = datetime.datetime.utcnow()

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

            p_compounds = pubchempy.get_compounds(compound.name, 'name')
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
