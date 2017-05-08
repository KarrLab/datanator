# -*- coding: utf-8 -*-

"""
:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-04
:Copyright: 2017, Karr Lab
:License: MIT
"""

from io import BytesIO
from kinetic_datanator.util import molecule_util
from wc_utils.backup import BackupManager
import datetime
import dateutil.parser
import json
import jxmlease
import os
import requests
import requests_cache
import sqlalchemy
import sqlalchemy.ext.declarative
import sqlalchemy.orm
import sys
import urllib
import zipfile


DEFAULT_DATABASE_FILENAME = os.path.join(os.path.dirname(__file__), 'cache', 'ecmdb.sqlite')
# :obj:`str`: default path for the sqlite database

DEFAULT_DATABASE_ARCNAME = 'ecmdb.sqlite'
# :obj:`str`: default name for the sqlite database within a gzip backup

SqlalchemyBase = sqlalchemy.ext.declarative.declarative_base()
# :obj:`Base`: SQL Alchemy base class for class definitions


def get_engine(filename=DEFAULT_DATABASE_FILENAME):
    """ Get an engine for the sqlite database

    Args:
        filename (:obj:`str`): path to sqlite database

    Returns:
        :obj:`sqlalchemy.engine.Engine`: database engine
    """
    if not os.path.isdir(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename))
    return sqlalchemy.create_engine('sqlite:///' + filename)


def get_session(engine=None, auto_download=True, auto_update=False, force_download=False, force_update=False,
                arcname=DEFAULT_DATABASE_ARCNAME, max_laws=float('inf')):
    """ Get a session for the sqlite database

    Args:
        engine (:obj:`sqlalchemy.engine.Engine`, optional): database engine
        arcname (:obj:`str`, optional): name for the sqlite database within a gzip backup
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
        arcname (:obj:`str`, optional): name for the sqlite database within a gzip backup
        max_laws (:obj:`int`, optional): maximum number of laws to download

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
        Downloader(session=session, max_laws=max_laws).download()

    return session

def backup(filename=DEFAULT_DATABASE_FILENAME, arcname=DEFAULT_DATABASE_ARCNAME):
    """ Backup the local sqlite database to the Karr Lab server

    Args:
        filename (:obj:`str`, optional): path to sqlite database
        arcname (:obj:`str`, optional): name for the sqlite database within a gzip backup
    """
    BackupManager(filename, arcname=arcname) \
        .create() \
        .upload() \
        .cleanup()


def download(filename=DEFAULT_DATABASE_FILENAME, arcname=DEFAULT_DATABASE_ARCNAME):
    """ Download the local sqlite database from the Karr Lab server

    Args:
        filename (:obj:`str`, optional): path to sqlite database
        arcname (:obj:`str`, optional): name for the sqlite database within a gzip backup
    """
    BackupManager(filename, arcname=arcname) \
        .download() \
        .extract() \
        .cleanup()

compound_compartment = sqlalchemy.Table(
    'compound_compartment', SqlalchemyBase.metadata,
    sqlalchemy.Column('compound__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('compound._id')),
    sqlalchemy.Column('compartment__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('compartment._id')),
)
# :obj:`sqlalchemy.Table`: Compound:Compartment many-to-many association table

compound_synonym = sqlalchemy.Table(
    'compound_synonym', SqlalchemyBase.metadata,
    sqlalchemy.Column('compound__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('compound._id')),
    sqlalchemy.Column('synonym__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('synonym._id')),
)
# :obj:`sqlalchemy.Table`: Compound:Synonym many-to-many association table

compound_resource = sqlalchemy.Table(
    'compound_resource', SqlalchemyBase.metadata,
    sqlalchemy.Column('compound__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('compound._id')),
    sqlalchemy.Column('resource__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('resource._id')),
)
# :obj:`sqlalchemy.Table`: Compound:Resource many-to-many association table

concentration_resource = sqlalchemy.Table(
    'concentration_resource', SqlalchemyBase.metadata,
    sqlalchemy.Column('concentration__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('concentration._id')),
    sqlalchemy.Column('resource__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('resource._id')),
)
# :obj:`sqlalchemy.Table`: Concentration:Resource many-to-many association table


class Synonym(SqlalchemyBase):
    """
    Args:
        name
        compouds
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    name = sqlalchemy.Column(sqlalchemy.String(), unique=True, index=True)
    compounds = sqlalchemy.orm.relationship('Compound', secondary=compound_synonym, back_populates='synonyms')

    __tablename__ = 'synonym'


class Compartment(SqlalchemyBase):
    """
    Attributes:
        name
        compounds
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    name = sqlalchemy.Column(sqlalchemy.String(), unique=True, index=True)
    compounds = sqlalchemy.orm.relationship('Compound', secondary=compound_compartment, back_populates='compartments')

    __tablename__ = 'compartment'


class Concentration(SqlalchemyBase):
    """
    Attributes:
        value
        error
        strain
        state
        media
        temperature
        system
        compound
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    value = sqlalchemy.Column(sqlalchemy.Float())
    error = sqlalchemy.Column(sqlalchemy.Float())
    strain = sqlalchemy.Column(sqlalchemy.String())
    growth_status = sqlalchemy.Column(sqlalchemy.String())
    media = sqlalchemy.Column(sqlalchemy.String())
    temperature = sqlalchemy.Column(sqlalchemy.Float())
    growth_system = sqlalchemy.Column(sqlalchemy.String())
    compound_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('compound._id'))
    compound = sqlalchemy.orm.relationship('Compound', back_populates='concentrations', foreign_keys=[compound_id])
    references = sqlalchemy.orm.relationship('Resource', secondary='concentration_resource', back_populates='concentrations')

    __tablename__ = 'concentration'


class Resource(SqlalchemyBase):
    """ Represents an external resource

    Attributes:
        namespace (:obj:`str`): external namespace
        id (:obj:`str`): external identifier
        compounds (:obj:`list` of :obj:`Compound`): compounds
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    namespace = sqlalchemy.Column(sqlalchemy.String())
    id = sqlalchemy.Column(sqlalchemy.String())
    compounds = sqlalchemy.orm.relationship('Compound', secondary=compound_resource, back_populates='cross_references')
    concentrations = sqlalchemy.orm.relationship('Concentration', secondary='concentration_resource', back_populates='references')

    sqlalchemy.schema.UniqueConstraint(namespace, id)

    __tablename__ = 'resource'


class Compound(SqlalchemyBase):
    """ Represents an ECMDB entry

    Attributes:
        id
        name
        synonyms
        description
        structure
        _structure_formula_connectivity
        compartments
        concentrations
        cross_references
        comment
        created
        updated
        downloaded
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    id = sqlalchemy.Column(sqlalchemy.String())
    name = sqlalchemy.Column(sqlalchemy.String(), index=True)
    synonyms = sqlalchemy.orm.relationship('Synonym', secondary=compound_synonym, back_populates='compounds')
    description = sqlalchemy.Column(sqlalchemy.Text())
    structure = sqlalchemy.Column(sqlalchemy.Text())
    _structure_formula_connectivity = sqlalchemy.Column(sqlalchemy.Text(), index=True)
    compartments = sqlalchemy.orm.relationship('Compartment', secondary=compound_compartment, back_populates='compounds')
    concentrations = sqlalchemy.orm.relationship('Concentration', back_populates='compound', foreign_keys=[Concentration.compound_id])
    cross_references = sqlalchemy.orm.relationship('Resource', secondary=compound_resource, back_populates='compounds')
    comment = sqlalchemy.Column(sqlalchemy.Text())
    created = sqlalchemy.Column(sqlalchemy.DateTime)
    updated = sqlalchemy.Column(sqlalchemy.DateTime)
    downloaded = sqlalchemy.Column(sqlalchemy.DateTime, default=datetime.datetime.utcnow())

    __tablename__ = 'compound'


class Downloader(object):
    DOWNLOAD_INDEX_URL = 'http://ecmdb.ca/download/ecmdb.json.zip'
    # :obj:`str`: URL to download an index of ECMDB

    DOWNLOAD_COMPOUND_URL = 'http://ecmdb.ca/compounds/{}.xml'
    # :obj:`str`: URL pattern to download an ECMDB compound entry

    def __init__(self, session, requests_cache_name=None, max_entries=float('inf'), verbose=False):
        self.session = session
        self.max_entries = max_entries

        if requests_cache_name is None:
            requests_cache_dir = os.path.join(os.path.dirname(__file__), 'cache')
            if not os.path.isdir(requests_cache_dir):
                os.makedirs(requests_cache_dir)
            requests_cache_name = os.path.join(requests_cache_dir, 'ecmdb.requests.py{}'.format(sys.version_info[0]))
        self.requests_cache_name = requests_cache_name
        self.verbose = verbose

        # initialize cache, database
        self.init_requests_cache()
        self.init_database()

    def init_requests_cache(self):
        """ Setup the cache for SABIO-RK HTTP requests """
        pass

    def clear_requests_cache(self):
        """ Clear the cahce for SABIO-RK HTTP requests """
        session = requests_cache.core.CachedSession(self.requests_cache_name, backend='sqlite', expire_after=None)
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

        db_session = self.session
        req_session = requests_cache.core.CachedSession(self.requests_cache_name, backend='sqlite', expire_after=None)

        # download content from server
        if self.verbose:
            print('Downloading compound IDs ...')

        response = req_session.get(self.DOWNLOAD_INDEX_URL)
        response.raise_for_status()

        if self.verbose:
            print('  done')

        # unzip and parse content
        if self.verbose:
            print('Parsing compound IDs ...')

        with zipfile.ZipFile(BytesIO(response.content), 'r') as zip_file:
            with zip_file.open('ecmdb.json', 'r') as json_file:
                entries = json.load(json_file)

        if self.verbose:
            print('  found {} compounds'.format(len(entries)))

        # sort entires
        entries.sort(key=lambda e: e['m2m_id'])

        # limit number of processed entries
        if len(entries) > self.max_entries:
            entries = entries[0:self.max_entries]

        # load content into sqlite database
        if self.verbose:
            print('Downloading {} compounds ...'.format(len(entries)))

        xml_parser = jxmlease.Parser()
        for i_entry, entry in enumerate(entries):
            if self.verbose and (i_entry % 100 == 0):
                print('  Downloading compound {} of {}'.format(i_entry + 1, len(entries)))
            print(entry['m2m_id'])  # todo: remove

            # get details
            response = req_session.get(self.DOWNLOAD_COMPOUND_URL.format(entry['m2m_id']))
            response.raise_for_status()
            entry_details = xml_parser(response.text)['compound']

            compound = self.get_or_create_object(Compound, id=self.get_node_text(entry_details['m2m_id']))

            if 'name' in entry_details:
                compound.name = self.get_node_text(entry_details['name'])

            if 'description' in entry_details:
                compound.description = self.get_node_text(entry_details['description'])

            compound.structure = self.get_node_text(entry_details['inchi'])

            compound.comment = entry['comment']

            compound.created = dateutil.parser.parse(self.get_node_text(entry_details['creation_date'])).replace(tzinfo=None)
            compound.updated = dateutil.parser.parse(self.get_node_text(entry_details['update_date'])).replace(tzinfo=None)

            # calculate core InChI layers to facilitate searching
            compound._structure_formula_connectivity = molecule_util.InchiMolecule(compound.structure) \
                .get_formula_and_connectivity(hydrogen=False)

            # synonyms
            compound.synonyms = []

            if 'iupac_name' in entry_details:
                node = entry_details['iupac_name']
                name = self.get_node_text(node)
                compound.synonyms.append(self.get_or_create_object(Synonym, name=name))

            if 'traditional_iupac' in entry_details:
                node = entry_details['traditional_iupac']
                name = self.get_node_text(node)
                compound.synonyms.append(self.get_or_create_object(Synonym, name=name))

            parent_node = entry_details['synonyms']
            if 'synonym' in parent_node:
                nodes = self.get_node_children(parent_node, 'synonym')
                for node in nodes:
                    name = self.get_node_text(node)
                    compound.synonyms.append(self.get_or_create_object(Synonym, name=name))

            # locations
            compound.compartments = []
            parent_node = entry_details['cellular_locations']
            if 'cellular_location' in parent_node:
                nodes = self.get_node_children(parent_node, 'cellular_location')
                for node in nodes:
                    name = self.get_node_text(node)
                    compound.compartments.append(self.get_or_create_object(Compartment, name=name))

            # todo: experimental properties
            # * state
            # * melting_point
            # * water_solubility
            # * logp_hydrophobicity

            # concentrations
            compound.concentrations = []
            parent_node = entry_details['concentrations']
            if 'concentration' in parent_node:
                values = self.get_node_children(parent_node, 'concentration')
                errors = self.get_node_children(parent_node, 'error')
                units = self.get_node_children(parent_node, 'concentration_units')
                strains = self.get_node_children(parent_node, 'strain')
                statuses = self.get_node_children(parent_node, 'growth_status')
                medias = self.get_node_children(parent_node, 'growth_media')
                temperatures = self.get_node_children(parent_node, 'temperature')
                systems = self.get_node_children(parent_node, 'growth_system')
                references = self.get_node_children(parent_node, 'reference')

                for i in range(len(values)):
                    value = float(self.get_node_text(values[i]))
                    error = float(self.get_node_text(errors[i]) or 'nan')
                    unit = self.get_node_text(units[i])
                    if unit == 'uM':
                        pass
                    else:
                        raise ValueError('Unsupport units: {}'.format(unit))

                    if temperatures[i]:
                        temperature, unit = self.get_node_text(temperatures[i]).split(' ')
                        temperature = float(temperature)
                        if unit != 'oC':
                            raise ValueError('Unsupport units: {}'.format(unit))
                    else:
                        temperature = None

                    concentration = Concentration(
                        value=value,
                        error=error,
                        strain=self.get_node_text(strains[i]) or None,
                        growth_status=self.get_node_text(statuses[i]) or None,
                        media=self.get_node_text(medias[i]) or None,
                        temperature=temperature,
                        growth_system=self.get_node_text(systems[i]) or None,
                    )
                    db_session.add(concentration)

                    if 'pubmed_id' in references[i]:
                        pmid_nodes = self.get_node_children(references[i], 'pubmed_id')
                        for node in pmid_nodes:
                            id = self.get_node_text(node)
                            concentration.references.append(self.get_or_create_object(Resource, namespace='pubmed', id=id))

                    compound.concentrations.append(concentration)

            # cross references
            compound.cross_references = []

            id = self.get_node_text(entry_details['biocyc_id'])
            if id:
                compound.cross_references.append(self.get_or_create_object(Resource, namespace='biocyc', id=id))

            id = self.get_node_text(entry_details['cas_registry_number'])
            if id:
                compound.cross_references.append(self.get_or_create_object(Resource, namespace='cas', id=id))

            id = self.get_node_text(entry_details['chebi_id'])
            if id:
                compound.cross_references.append(self.get_or_create_object(Resource, namespace='chebi', id='CHEBI:' + id))

            id = self.get_node_text(entry_details['chemspider_id'])
            if id:
                compound.cross_references.append(self.get_or_create_object(Resource, namespace='chemspider', id=id))

            id = self.get_node_text(entry_details['foodb_id'])
            if id:
                compound.cross_references.append(self.get_or_create_object(Resource, namespace='foodb.compound', id=id))

            id = self.get_node_text(entry_details['het_id'])
            if id:
                compound.cross_references.append(self.get_or_create_object(Resource, namespace='ligandexpo', id=id))

            id = self.get_node_text(entry_details['hmdb_id'])
            if id:
                compound.cross_references.append(self.get_or_create_object(Resource, namespace='hmdb', id=id))

            id = self.get_node_text(entry_details['kegg_id'])
            if id:
                compound.cross_references.append(self.get_or_create_object(Resource, namespace='kegg.compound', id=id))

            id = self.get_node_text(entry_details['msds_url'])
            if id:
                compound.cross_references.append(self.get_or_create_object(Resource, namespace='msds.url', id=id))

            id = self.get_node_text(entry_details['pubchem_compound_id'])
            if id:
                compound.cross_references.append(self.get_or_create_object(Resource, namespace='pubchem.compound', id=id))

            id = self.get_node_text(entry_details['wikipidia'])
            if id:
                compound.cross_references.append(self.get_or_create_object(Resource, namespace='wikipedia.en', id=id))

            # add to session
            db_session.add(compound)

        if self.verbose:
            print('  done')

        # commit changes to database
        if self.verbose:
            print('Saving database ...')

        db_session.commit()

        if self.verbose:
            print('  done')

    def get_or_create_object(self, cls, **kwargs):
        q = self.session.query(cls).filter_by(**kwargs)
        if q.count():
            return q.first()
        else:
            obj = cls(**kwargs)
            self.session.add(obj)
            return obj

    def get_node_children(self, node, children_name):
        nodes = node[children_name]
        if isinstance(nodes, jxmlease.listnode.XMLListNode):
            return nodes
        return [nodes]

    def get_node_text(self, node):
        if isinstance(node, jxmlease.cdatanode.XMLCDATANode):
            return node.get_cdata()
        return node
