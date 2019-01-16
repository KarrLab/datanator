# -*- coding: utf-8 -*-

"""
:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Zhouyang Lian <zhouyang.lian@familian.life>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2019-01-08
:Copyright: 2017, Karr Lab
:License: MIT
"""

from datanator.core import data_source
from datanator.util import molecule_util
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
import sdf
import numpy as np

Base = sqlalchemy.ext.declarative.declarative_base()
# :obj:`Base`: base model for local sqlite database

compound_compartment = sqlalchemy.Table(
    'compound_compartment', Base.metadata,
    sqlalchemy.Column('compound__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('compound._id'), index=True),
    sqlalchemy.Column('compartment__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('compartment._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Compound:Compartment many-to-many association table

compound_synonym = sqlalchemy.Table(
    'compound_synonym', Base.metadata,
    sqlalchemy.Column('compound__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('compound._id'), index=True),
    sqlalchemy.Column('synonym__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('synonym._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Compound:Synonym many-to-many association table

compound_resource = sqlalchemy.Table(
    'compound_resource', Base.metadata,
    sqlalchemy.Column('compound__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('compound._id'), index=True),
    sqlalchemy.Column('resource__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('resource._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Compound:Resource many-to-many association table

concentration_resource = sqlalchemy.Table(
    'concentration_resource', Base.metadata,
    sqlalchemy.Column('concentration__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('concentration._id'), index=True),
    sqlalchemy.Column('resource__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('resource._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Concentration:Resource many-to-many association table


class Synonym(Base):
    """ Represents a synonym

    Args:
        name (:obj:`str`): name
        compounds (:obj:`list` of :obj:`Compound`): list of compounds
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    name = sqlalchemy.Column(sqlalchemy.String(), unique=True, index=True)

    __tablename__ = 'synonym'


class Compartment(Base):
    """ Represents a compartment

    Attributes:
        name (:obj:`str`): name
        compounds (:obj:`list` of :obj:`Compound`): list of compounds
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    name = sqlalchemy.Column(sqlalchemy.String(), unique=True, index=True)

    __tablename__ = 'compartment'


class Concentration(Base):
    """ Represents an observed concentration

    Attributes:
        compound (:obj:`Compound`): compound
        value (:obj:`float`): value in uM
        error (:obj:`float`): error in uM
        strain (:obj:`str`): observed strain
        growth_status (:obj:`str`): observed growth status (e.g. exponential phase, log phase, etc.)
        media (:obj:`str`): observed media
        temperaturer (:obj:`float`): temperature in C
        growth_system (:obj:`str`): observed growth system (e.g. chemostat, 384 well plate, etc.)
        references (:obj:`list` of :obj:`Resource`): list of references
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    compound_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('compound._id'), index=True)
    value = sqlalchemy.Column(sqlalchemy.Float())
    error = sqlalchemy.Column(sqlalchemy.Float())
    strain = sqlalchemy.Column(sqlalchemy.String())
    growth_status = sqlalchemy.Column(sqlalchemy.String())
    media = sqlalchemy.Column(sqlalchemy.String())
    temperature = sqlalchemy.Column(sqlalchemy.Float())
    growth_system = sqlalchemy.Column(sqlalchemy.String())
    references = sqlalchemy.orm.relationship('Resource', secondary='concentration_resource',
                                             backref=sqlalchemy.orm.backref('concentrations'))

    __tablename__ = 'concentration'


class Resource(Base):
    """ Represents an external resource

    Attributes:
        namespace (:obj:`str`): external namespace
        id (:obj:`str`): external identifier
        compounds (:obj:`list` of :obj:`Compound`): compounds
        concentrations (:obj:`list` of :obj:`Concentration`): concentrations
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    namespace = sqlalchemy.Column(sqlalchemy.String())
    id = sqlalchemy.Column(sqlalchemy.String())

    sqlalchemy.schema.UniqueConstraint(namespace, id)

    __tablename__ = 'resource'


class Compound(Base):
    """ Represents an YMDB entry
	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5210545/
	
    Attributes:
        id (:obj:`str`): YMDB identifier
        name (:obj:`str`): name
        synonyms (:obj:`list` of :obj:`Synonym`): synonyms
        description (:obj:`str`): description
        structure (:obj:`str`): structure in InChI format
        _structure_formula_connectivity (:obj:`str`): empiral formula and connectivity InChI layers; used to
            quickly search for compound structures
        compartments (:obj:`list` of :obj:`Compartment`): compartments
        concentrations (:obj:`list` of :obj:`Concentration`): concentrations
        cross_references (:obj:`list` of :obj:`Resources`): cross references
        comment (:obj:`str`): internal YMDB comments about the entry
        created (:obj:`datetime.datetime`): time that the entry was created in YMDB
        updated (:obj:`datetime.datetime`): time that the entry was last updated in YMDB
        downloaded (:obj:`datetime.datetime`): time that the entry was downloaded from YMDB
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    id = sqlalchemy.Column(sqlalchemy.String())
    name = sqlalchemy.Column(sqlalchemy.String(), index=True)
    synonyms = sqlalchemy.orm.relationship('Synonym', secondary=compound_synonym, backref=sqlalchemy.orm.backref('compounds'))
    description = sqlalchemy.Column(sqlalchemy.Text())
    structure = sqlalchemy.Column(sqlalchemy.Text())
    _structure_formula_connectivity = sqlalchemy.Column(sqlalchemy.Text(), index=True)
    compartments = sqlalchemy.orm.relationship('Compartment', secondary=compound_compartment, backref=sqlalchemy.orm.backref('compounds'))
    concentrations = sqlalchemy.orm.relationship('Concentration',
                                                 foreign_keys=[Concentration.compound_id],
                                                 backref=sqlalchemy.orm.backref('compound'),
                                                 cascade="all, delete-orphan")
    cross_references = sqlalchemy.orm.relationship('Resource', secondary=compound_resource, backref=sqlalchemy.orm.backref('compounds'))
    comment = sqlalchemy.Column(sqlalchemy.Text())
    created = sqlalchemy.Column(sqlalchemy.DateTime)
    updated = sqlalchemy.Column(sqlalchemy.DateTime)
    downloaded = sqlalchemy.Column(sqlalchemy.DateTime, default=datetime.datetime.utcnow())

    __tablename__ = 'compound'


class Ymdb(data_source.HttpDataSource):
    """ A local sqlite copy of the YMDB database

    Attributes:
        DOWNLOAD_FULL_DB_URL (:obj:`str`): URL to download the full database of YMDB
        DOWNLOAD_STRUCTURES_URL (:obj:`str`): URL pattern to download all metabolite structures of YMDB
        DOWNLOAD_PROTEIN_SEQUENCE_URL (:obj:`str`): URL pattern to download all protein sequence of YMDB
        DOWNLOAD_GENE_SEQUENCE_URL (:obj:`str`): URL pattern to download all genetic sequence of YMDB
    """

    base_model = Base
    ENDPOINT_DOMAINS = {
        # 'ecmdb': 'http://ecmdb.ca',
        'ymdb': 'http://ymdb.ca',
    }
    DOWNLOAD_FULL_DB_URL = ENDPOINT_DOMAINS['ymdb'] + '/system/downloads/current/ymdb.json.zip'
    DOWNLOAD_STRUCTURES_URL = ENDPOINT_DOMAINS['ymdb'] + '/system/downloads/current/ymdb.sdf.zip'
    DOWNLOAD_COMPOUND_URL = ENDPOINT_DOMAINS['ymdb'] + '/compounds/{}.xml'
    DOWNLOAD_COMPOUND_STRUCTURE_URL = ENDPOINT_DOMAINS['ymdb'] + '/structures/compounds/{}.inchi'
    DOWNLOAD_PROTEIN_SEQUENCE_URL = ENDPOINT_DOMAINS['ymdb'] + '/system/downloads/current/protein_sequences.fasta.zip'
    DOWNLOAD_GENE_SEQUENCE_URL = ENDPOINT_DOMAINS['ymdb'] + '/system/downloads/current/gene_sequences.fasta.zip'

    def load_content(self):
        """ Download the content of YMDB and store it to a local sqlite database. """
        db_session = self.session
        req_session = self.requests_session

        # download content from server
        if self.verbose:
            print('Downloading compound IDs ...')

        response = req_session.get(self.DOWNLOAD_FULL_DB_URL)
        response.raise_for_status()

        if self.verbose:
            print('  done')

        # unzip and parse content
        if self.verbose:
            print('Parsing compound IDs ...')

        with zipfile.ZipFile(io.BytesIO(response.content), 'r') as zip_file:
            with zip_file.open('ymdb.json', 'r') as json_file:
                entries = json.load(json_file)

        if self.verbose:
            print('  found {} compounds'.format(len(entries)))

        # sort entires
        entries.sort(key=lambda e: e['ymdb_id'])

        # limit number of processed entries
        if len(entries) > self.max_entries:
            entries = entries[0:self.max_entries]

        # load content into sqlite database
        if self.verbose:
            print('Downloading {} compounds ...'.format(len(entries)))

        xml_parser = jxmlease.Parser()
        for i_entry, entry in enumerate(entries):
            if self.verbose and (i_entry % 10 == 0):
                print('  Downloading compound {} of {}'.format(i_entry + 1, len(entries)))

            # get details
            response = req_session.get(self.DOWNLOAD_COMPOUND_URL.format(entry['ymdb_id']))
            try:
                response.raise_for_status()
            except requests.exceptions.HTTPError:
                warnings.warn('Unable to download data for compound {}'.format(entry['ymdb_id']), data_source.DataSourceWarning)
                continue

            entry_details = xml_parser(response.text)['compound']

            compound = self.get_or_create_object(Compound, id=self.get_node_text(entry_details['ymdb_id']))

            if 'name' in entry_details:
                compound.name = self.get_node_text(entry_details['name'])

            if 'description' in entry_details:
                compound.description = self.get_node_text(entry_details['description'])

            compound.structure = self.get_node_text(entry_details['inchi'])
            if not compound.structure:
                response2 = req_session.get(self.DOWNLOAD_COMPOUND_STRUCTURE_URL.format(entry['ymdb_id']))
                response2.raise_for_status()
                compound.structure = response2.text

            compound.comment = entry.get('comment','default')

            compound.created = dateutil.parser.parse(self.get_node_text(entry_details['creation_date'])).replace(tzinfo=None)
            compound.updated = dateutil.parser.parse(self.get_node_text(entry_details['update_date'])).replace(tzinfo=None)

            # calculate core InChI layers to facilitate searching
            try:
                compound._structure_formula_connectivity = molecule_util.InchiMolecule(compound.structure) \
                    .get_formula_and_connectivity()
            except ValueError:
                warnings.warn('Unable to encode structure for {} in InChI'.format(entry['ymdb_id']), data_source.DataSourceWarning)
                compound._structure_formula_connectivity = None

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

            # todo (enhancement): parse experimental properties
            # * state
            # * melting_point
            # * water_solubility
            # * logp_hydrophobicity

            # concentrations
            compound.concentrations = []
            parent_node = entry_details['concentrations']
            if 'concentration' in parent_node:
                values = self.get_node_children(parent_node, 'concentration')
                if 'error' in parent_node:  
                    errors = self.get_node_children(parent_node, 'error')
                if 'concentration_units' in parent_node:
                    units = self.get_node_children(parent_node, 'concentration_units')
                if 'strain' in parent_node: 
                    strains = self.get_node_children(parent_node, 'strain')
                if 'growth_status' in parent_node: 
                    statuses = self.get_node_children(parent_node, 'growth_status')
                if 'growth_media' in parent_node: 
                    medias = self.get_node_children(parent_node, 'growth_media')
                if 'temperature' in parent_node: 
                    temperatures = self.get_node_children(parent_node, 'temperature')
                if 'growth_system' in parent_node: 
                    systems = self.get_node_children(parent_node, 'growth_system')
                if 'reference' in parent_node: 
                    references = self.get_node_children(parent_node, 'reference')

                for i_conc in range(len(values)):
                    value = float(self.get_node_text(values[i_conc]))
                    error = float(self.get_node_text(errors[i_conc]) or 'nan')
                    unit = self.get_node_text(units[i_conc])
                    if unit == 'uM' or unit == '&#181;M':
                        pass
                    else:
                        raise ValueError('Unsupport units: {}'.format(unit))

                    if temperatures[i_conc]:
                        temperature, unit = self.get_node_text(temperatures[i_conc]).split(' ')
                        temperature = float(temperature)
                        if unit != 'oC':
                            raise ValueError('Unsupport units: {}'.format(unit))
                    else:
                        temperature = None

                    concentration = Concentration(
                        value=value,
                        error=error,
                        strain=self.get_node_text(strains[i_conc]) or None,
                        growth_status=self.get_node_text(statuses[i_conc]) or None,
                        media=self.get_node_text(medias[i_conc]) or None,
                        temperature=temperature,
                        growth_system=self.get_node_text(systems[i_conc]) or None,
                    )
                    db_session.add(concentration)

                    if 'pubmed_id' in references[i_conc]:
                        pmid_nodes = self.get_node_children(references[i_conc], 'pubmed_id')
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

            id = self.get_node_text(entry_details.get('het_id','default'))
            if id:
                compound.cross_references.append(self.get_or_create_object(Resource, namespace='ligandexpo', id=id))

            id = self.get_node_text(entry_details['hmdb_id'])
            if id:
                compound.cross_references.append(self.get_or_create_object(Resource, namespace='hmdb', id=id))

            id = self.get_node_text(entry_details['kegg_id'])
            if id:
                compound.cross_references.append(self.get_or_create_object(Resource, namespace='kegg.compound', id=id))

            id = self.get_node_text(entry_details.get('msds_url','default'))
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

            if self.commit_intermediate_results and (i_entry % 100 == 99):
                db_session.commit()

        if self.verbose:
            print('  done')

        # commit changes to database
        if self.verbose:
            print('Saving database ...')

        db_session.commit()

        if self.verbose:
            print('  done')

    def get_node_children(self, node, children_name):
        """ Get the children of an XML node

        Args:
            node (:obj:`jxmlease.cdatanode.XMLNode`): XML node
            children_name (:obj:`str`): tag names of the desired children

        Returns:
            :obj:`list` of :obj:`XMLNode`: list of child nodes
        """

        nodes = node[children_name]
        if isinstance(nodes, jxmlease.listnode.XMLListNode):
            return nodes
        return [nodes]

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
