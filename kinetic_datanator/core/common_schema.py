# -*- coding: utf-8 -*-

"""
This code is a common schema for all the kinetic_datanator modules

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-07-31
:Copyright: 2017, Karr Lab
:License: MIT
"""

##TODO: Create the definition of tables for the format and then include with the query function in sql

from sqlalchemy import Column, BigInteger, Integer, Float, String, Text, ForeignKey, Boolean, Table, create_engine, Numeric
from sqlalchemy.orm import relationship, backref, sessionmaker
from kinetic_datanator.core import data_source
from kinetic_datanator.data_source import corum, pax, jaspar, array_express, ecmdb, sabio_rk
import sqlalchemy.ext.declarative
from six import BytesIO
import six

Base = sqlalchemy.ext.declarative.declarative_base()

class Taxon(Base):
    """ Represents a Taxon

    Attributes:
        ncbi_id (:obj:`int`): `NCBI taxonomy <https://www.ncbi.nlm.nih.gov/taxonomy>`_ id
        species_name (:obj:`str`): Name of Species

    """
    __tablename__ = 'taxon'

    ncbi_id = Column(Integer, primary_key = True)
    species_name = Column(String(255), unique = True)


class CellLine(Base):
    """ Represents a CellLine


    """
    __tablename__ = 'celline'

    id = Column(Integer, primary_key = True)
    name = Column(String)

class Protein(Base):
    """


    """

    __tablename__ = 'protein'

    id = Column(Integer, primary_key = True)
    name = Column(String)
    protein_id = Column(String)

    ##Insert association table for subunits

class Subunit(Base):
    """

    """
    __tablename__ = 'subunit'

    id = Column(Integer, primary_key = True)
    uniprot_id = Column(String)
    entrez_id = Column(Integer)

class Class(Base):
    """

    """
    __tablename__ = 'class'

    id = Column(Integer, primary_key = True)
    name = Column(String)

class Family(Base):
    """

    """
    __tablename__ = 'family'

    id = Column(Integer, primary_key = True)
    name = Column(String)

class Resource(Base):
    """ Represents a Resource


    """

    __tablename__ = 'resource'
    _id = Column(Integer, primary_key = True)
    namespace = Column(String)
    id = Column(String, unique = True)

class Synonyms(Base):
    """

    """
    __tablename__ = 'synonym'
    id = Column(Integer, primary_key = True)
    name = Column(String)

class Compartment(Base):
    """

    """
    __tablename__ = 'compartment'
    id = Column(Integer, primary_key = True)
    name = Column(String, unique = True, index = True)

class Parameter(Base):
    """

    """
    __tablename__ = 'parameter'
    id = Column(Integer, primary_key = True)
    value = Column(Numeric)
    unit = Column(String)

class Reaction(Base):
    """

    """
    __tablename__ = 'reaction'
    id = Column(Integer, primary_key = True)
    reactants = Column(String)
    products = Column(String)
    enzymes = Column(String) ##Link to protien category if protein

class Structure(Base):
    """

    """
    __tablename__ = 'structure'
    id = Column(Integer, primary_key = True)
    inchi = Column(String) #Figure out this inchi buisness

class Conditions(Base):
    """

    """
    __tablename__ = 'conditions'
    id = Column(Integer, primary_key = True)
    PH = Column(Float)
    temperature = Column(Float)
    ##ETC Growth media

class Method(Base):
    """

    """
    __tablename__ = 'method'
    id = Column(Integer, primary_key = True)
    name = Column(String)

class Entry(Base):
    """

    """
    __tablename__= 'entry'
    id = Column(Integer, primary_key = True)

class Database(Base):
    """
    """
    __tablename__ = 'database'
    id = Column(Integer, primary_key = True)
    name = Column(String, unique = True)


class CommonSchema(data_source.CachedDataSource):
    """
    A Local SQLlite copy of the aggregation of data_source modules

    NOTES: For now.. Keeping the SQL local, out of cache
    """
    base_model = Base

    def load_content(self):

        session = self.session

        """ WORKING """
        # When truly implementing the db will be load from Karr Cache: self.cache_dirname

        corumdb = corum.Corum(cache_dirname = self.cache_dirname , clear_content = True, load_content=False, download_backup=False, max_entries = 5)
        # corumdb = corum.Corum(name = 'corum', clear_content = True, load_content=False, download_backup=False, max_entries = 5)
        corumdb.load_content()
        corum_session = corumdb.session

        self.corum_taxon = corum_session.query(corum.Taxon).all()
        self.corum_observation = corum_session.query(corum.Observation).all()
        self.corum_complex = corum_session.query(corum.Complex).all()
        self.corum_subunit = corum_session.query(corum.Subunit).all()

        paxdb = pax.Pax(cache_dirname = self.cache_dirname, clear_content = True,  load_content=False, download_backup=False, max_entries = 5)
        # paxdb = pax.Pax(name = 'pax', clear_content = True,  load_content=False, download_backup=False, max_entries = 5)
        paxdb.load_content()
        pax_session = paxdb.session

        self.pax_taxon = pax_session.query(pax.Taxon).all()
        self.pax_dataset = pax_session.query(pax.Dataset).all()
        self.pax_protein = pax_session.query(pax.Protein).all()
        self.pax_observation = pax_session.query(pax.Observation).all()

        jaspardb = jaspar.Jaspar(cache_dirname = self.cache_dirname, clear_content = True, load_content=False, download_backup=False, max_entries = 5)
        # jaspardb = jaspar.Jaspar(name = 'jaspar', clear_content = True, load_content=False, download_backup=False, max_entries = 5)
        jaspardb.load_content()
        jaspar_session = jaspardb.session

        self.jaspar_transcriptionfactor = jaspar_session.query(jaspar.TranscriptionFactor).all()
        self.jaspar_matrix = jaspar_session.query(jaspar.Matrix).all()
        self.jaspar_matrixposition = jaspar_session.query(jaspar.MatrixPosition).all()
        self.jaspar_collection = jaspar_session.query(jaspar.Collection).all()
        self.jaspar_type = jaspar_session.query(jaspar.Type).all()
        self.jaspar_family = jaspar_session.query(jaspar.Family).all()
        self.jaspar_class = jaspar_session.query(jaspar.Class).all()
        self.jaspar_species = jaspar_session.query(jaspar.Species).all()
        self.jaspar_resource = jaspar_session.query(jaspar.Resource).all()
        self.jaspar_subunit = jaspar_session.query(jaspar.Subunit).all()

        ecmDB = ecmdb.Ecmdb(cache_dirname = self.cache_dirname, clear_content = True, load_content=False, download_backup=False, max_entries = 5)
        # ecmDB = ecmdb.Ecmdb(name = 'ecmdb', clear_content = True, load_content=False, download_backup=False, max_entries = 20)
        ecmDB.load_content()
        ecmdb_session = ecmDB.session

        self.ecmdb_synonym = ecmdb_session.query(ecmdb.Synonym).all()
        self.ecmdb_compartment = ecmdb_session.query(ecmdb.Compartment).all()
        self.ecmdb_concentration = ecmdb_session.query(ecmdb.Concentration).all()
        self.ecmdb_resource = ecmdb_session.query(ecmdb.Resource).all()
        self.ecmdb_compound = ecmdb_session.query(ecmdb.Compound).all()

        sabiodb = sabio_rk.SabioRk(cache_dirname = self.cache_dirname, clear_content = True, load_content=False, download_backup=False, max_entries = 5)
        # sabiodb = sabio_rk.SabioRk(name = 'sabio', clear_content = True, load_content=False, download_backup=False, max_entries = 5)
        sabiodb.load_content()
        sabio_session = sabiodb.session

        self.sabio_synonym = sabio_session.query(sabio_rk.Synonym).all()
        self.sabio_resource = sabio_session.query(sabio_rk.Resource).all()
        self.sabio_entry = sabio_session.query(sabio_rk.Entry).all()
        self.sabio_reactionparticipant = sabio_session.query(sabio_rk.ReactionParticipant).all()
        self.sabio_parameter = sabio_session.query(sabio_rk.Parameter).all()
        self.sabio_kineticlaw = sabio_session.query(sabio_rk.KineticLaw).all()
        self.sabio_compoundstructure = sabio_session.query(sabio_rk.CompoundStructure).all()


        """ NOT WORKING """
        #TODO: AE needs 1. normalization to schema(max_entries implementation/working functions) 2. working tests

        # arrayexpressdb = array_express.ArrayExpress(name = 'array_express', load_content=False, download_backup=False, max_entries = 5)
        # arrayexpressdb.load_experiments_from_text()
        # array_session = arrayexpressdb.session

        self.add_content_to_db()

    def add_content_to_db(self):

        # Taxon
        temp_taxon = self.jaspar_species + self.pax_taxon + self.corum_taxon + self.sabio_kineticlaw
        for items in temp_taxon:
            taxon = self.get_or_create_object(Taxon, ncbi_id = items.ncbi_id)
            if hasattr(items, 'species_name'):
                taxon.species_name = items.species_name

        # Resource
        temp_resource = self.corum_observation + self.pax_dataset + self.jaspar_resource + self.ecmdb_resource + self.sabio_resource
        for items in temp_resource:
            if hasattr(items, 'pubmed_id'):
                resource = self.get_or_create_object(Resource, namespace = 'pubmed', id = items.pubmed_id)
            elif hasattr(items, 'publication'):
                resource = self.get_or_create_object(Resource, namespace = 'url', id = items.publication)
            else:
                resource = self.get_or_create_object(Resource, namespace = items.namespace, id = items.id)

        # Compartment
        temp_compartment = self.ecmdb_compartment
        for items in temp_compartment:
            compartment = self.get_or_create_object(Compartment, name = items.name)


        if self.verbose:
            print('Comitting')
        self.session.commit()
