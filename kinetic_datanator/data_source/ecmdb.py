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


DEFAULT_CACHE_DIRNAME = os.path.join(os.path.dirname(__file__), 'cache')
# :obj:`str`: default path for the cached content

DEFAULT_DATABASE_FILENAME = os.path.join(DEFAULT_CACHE_DIRNAME, 'ecmdb.sqlite')
# :obj:`str`: default path for the sqlite database

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


def get_session(engine=None):
    """ Get a session for the sqlite database

    Args:
        engine (:obj:`sqlalchemy.engine.Engine`, optional): database engine

    Returns:
        :obj:`sqlalchemy.orm.session.Session`: database session
    """
    if not engine:
        engine = get_engine()
    filename = str(engine.url).replace('sqlite:///', '')

    if not os.path.isfile(filename):
        downloader = Downloader(session)
        downloader.download()

    return session


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


class Synonym(SqlalchemyBase):
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    name = sqlalchemy.Column(sqlalchemy.String(), unique=True, index=True)
    compounds = sqlalchemy.orm.relationship('Compound', secondary=compound_synonym, back_populates='synonyms')

    __tablename__ = 'synonym'


class Compartment(SqlalchemyBase):
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    name = sqlalchemy.Column(sqlalchemy.String(), unique=True, index=True)
    compounds = sqlalchemy.orm.relationship('Compound', secondary=compound_compartment, back_populates='compartments')

    __tablename__ = 'compartment'


class Concentration(SqlalchemyBase):
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    compound_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('compound._id'))
    compound = sqlalchemy.orm.relationship('Compound', back_populates='concentrations', foreign_keys=[compound_id])

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

    sqlalchemy.schema.UniqueConstraint(namespace, id)

    __tablename__ = 'resource'


class Compound(SqlalchemyBase):
    """ Represents an ECMDB entry

    Attributes:
        id
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
    # :obj:`str`: URL to download ECMDB content

    DOWNLOAD_COMPOUND_URL = 'http://ecmdb.ca/compounds/{}.xml'

    def __init__(self, session, cache_dirname=DEFAULT_CACHE_DIRNAME, max_compounds=float('inf')):
        self.session = session
        self.cache_dirname = cache_dirname
        self.max_compounds = max_compounds

    def download(self):
        db_session = self.session
        req_session = requests_cache.core.CachedSession(os.path.join(
            self.cache_dirname, 'ecmdb.requests.py{}'.format(sys.version_info[0])),
            backend='sqlite', expire_after=None)

        # download content from server
        response = req_session.get(self.DOWNLOAD_URL)
        response.raise_for_status()

        # unzip and parse content
        with zipfile.ZipFile(BytesIO(response.content), 'r') as zip_file:
            with zip_file.open('ecmdb.json', 'r') as json_file:
                entries = json.load(json_file)

        # initialize sqlite database
        engine = db_session.get_bind()
        SqlalchemyBase.metadata.drop_all(engine)
        SqlalchemyBase.metadata.create_all(engine)

        # load content into sqlite database
        xml_parser = jxmlease.Parser()
        for entry in entries[0:max(self.max_compounds, len(entries))]:
            # create compound
            compound = Compound(
                reference_synthesis=entry['synthesis_reference'],
            )

            # get details
            response = req_session.get(self.DOWNLOAD_COMPOUND_URL.format(entry['m2m_id']))
            response.raise_for_status()
            data = xml_parser(response.text)['compound']

            compound = Compound(
                id=data['m2m_id'].get_cdata(),
                name=data['name'].get_cdata(),
                description=data['description'].get_cdata(),
                structure=data['inchi'].get_cdata(),
                comment=entry['comment'],
                created=dateutil.parser.parse(data['creation_date']),
                updated=dateutil.parser.parse(data['update_date']),
            )

            # calculate core InChI layers to facilitate searching
            compound._structure_formula_connectivity = molecule_util.InchiMolecule(compound.structure) \
                .get_formula_and_connectivity(hydrogen=False)

            # synonyms
            compound.synonyms = []
            for node in data['synonyms']['synonym']:
                compound.synonyms.append(self.get_or_create_object(Synonym, name=node.get_cdata()))

            # locations
            compound.compartments = []
            for node in data['cellular_locations']['cellular_location']:
                compound.compartments.append(self.get_or_create_object(Compartment, name=node.get_cdata()))

            # todo: experimental properties
            # * state
            # * melting_point
            # * water_solubility
            # * logp_hydrophobicity

            # concentrations
            compound.concentrations = []
            tmp = data['concentrations']
            values = tmp['concentration']
            errors = tmp['error']
            units = tmp['concentration_units']
            strains = tmp['strain']
            states = tmp['growth_state']
            medias = tmp['growth_media']
            temperatures = tmp['temperature']
            systems = tmp['growth_system']
            references = tmp['reference']
            for i in range(len(growth_medias)):
                value = float(values[i].get_cdata())
                error = float(errors[i].get_cdata() or 'nan')
                unit = units[i].get_cdata()
                if unit == 'uM':
                    pass
                else:
                    raise ValueError('Unsupport units: {}'.format(unit))

                if temperatures[i].get_cdata():
                    temperature, unit = split(temperatures[i].get_cdata())
                    temperature = float(temperature)
                    if unit != 'oC':
                        raise ValueError('Unsupport units: {}'.format(unit))
                else:
                    temperature = None

                concentration = Concentration(
                    value=value,
                    error=error,
                    strain=strains[i].get_cdata() or None,
                    state=states[i].get_cdata() or None,
                    media=medias[i].get_cdata() or None,
                    temperature=temperature,
                    system=systems[i].get_cdata() or None)

                for node in references[i]['pubmed_id']:
                    concentration.references.append(self.get_or_create_object(Resource, namespace='pubmed', id=node.get_cdata()))

                compound.concentrations.append(concentration)
                db_session.add(concentration)

            # cross references
            compound.cross_references = []
            if data['biocyc_id']:
                compound.cross_references.append(self.get_or_create_object(Resource, namespace='biocyc', id=data['biocyc_id'].get_cdata()))
            if data['cas_registry_number']:
                compound.cross_references.append(self.get_or_create_object(
                    Resource, namespace='cas', id=data['cas_registry_number'].get_cdata()))
            if data['chebi_id']:
                compound.cross_references.append(self.get_or_create_object(
                    Resource, namespace='chebi', id='CHEBI:' + data['chebi_id'].get_cdata()))
            if data['chemspider_id']:
                compound.cross_references.append(self.get_or_create_object(
                    Resource, namespace='chemspider', id=data['chemspider_id'].get_cdata()))
            if data['foodb_id']:
                compound.cross_references.append(self.get_or_create_object(
                    Resource, namespace='foodb.compound', id=data['foodb_id'].get_cdata()))
            if data['het_id']:
                compound.cross_references.append(self.get_or_create_object(
                    Resource, namespace='ligandexpo', id=data['het_id'].get_cdata()))
            if data['hmdb_id']:
                compound.cross_references.append(self.get_or_create_object(Resource, namespace='hmdb', id=data['hmdb_id'].get_cdata()))
            if data['kegg_id']:
                compound.cross_references.append(self.get_or_create_object(
                    Resource, namespace='kegg.compound', id=data['kegg_id'].get_cdata()))
            if data['msds_url']:
                compound.cross_references.append(self.get_or_create_object(
                    Resource, namespace='msds.url', id=data['msds_url'].get_cdata()))
            if data['pubchem_compound_id']:
                compound.cross_references.append(self.get_or_create_object(
                    Resource, namespace='pubchem.compound', id=data['pubchem_compound_id'].get_cdata()))
            if data['wikipidia']:
                compound.cross_references.append(self.get_or_create_object(
                    Resource, namespace='wikipedia.en', id=data['wikipidia'].get_cdata()))

            # add to session
            db_session.add(compound)

        # commit changes to database
        db_session.commit()

    def get_or_create_object(self, cls, **kwargs):
        q = self.session.query(Resource).filter_by(**kwargs)
        if q.count():
            return q.first()
        else:
            obj = cls(**kwargs)
            self.session.add(obj)
