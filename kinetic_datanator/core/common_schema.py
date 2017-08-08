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


class Entry(Base):
    """
    Represents an Entry within the Common-Schema


    """
    __tablename__ = 'entry'

    id = Column(Integer, primary_key = True)
    name = Column(String(255), index = True)

    taxon_id = Column(Integer, ForeignKey('taxon.ncbi_id'), index = True)
    taxon = relationship('Taxon', backref = 'entry')

    resource_id = Column(Integer, ForeignKey('resource._id'), index = True)
    resource = relationship('Resource', backref = 'entry')

    cellline_id = Column(Integer, ForeignKey('cellline.id'), index = True)
    cellline = relationship('CellLine', backref = 'entry')

    method_id = Column(Integer, ForeignKey('method.id'), index = True)
    method = relationship('Method', backref = 'entry')


class Taxon(Base):
    """ Represents a Taxon

    Attributes:
        ncbi_id (:obj:`int`): `NCBI taxonomy <https://www.ncbi.nlm.nih.gov/taxonomy>`_ id
        species_name (:obj:`str`): Name of Species

    """
    __tablename__ = 'taxon'

    ncbi_id = Column(Integer, primary_key = True)
    species_name = Column(String(255), unique = True)

class Resource(Base):
    """ Represents a Resource

    Attributes:
        namespace (:obj:`str`): Type of Resource
        id (:obj:`str`): Resource Identifier

    """
    __tablename__ = 'resource'
    _id = Column(Integer, primary_key = True)
    namespace = Column(String)
    id = Column(String, unique = True)

class CellCompartment(Base):
    """

    """
    __tablename__ = 'cellcompartment'
    id = Column(Integer, primary_key = True)
    name = Column(String, unique = True, index = True)


class PaxDataSet(Base):
    """

    """

    __tablename__ = 'pax_dataset'
    id  = Column(Integer, primary_key = True)
    file_name = Column(String, unique = True)
    score = Column(Float)
    weight = Column(Integer)
    coverage = Column(Integer)

    entry_id = Column(Integer, ForeignKey('entry.id'), index = True)
    entry = relationship('Entry', backref = 'pax_dataset')


class PaxAbundanceData(Base):
    """

    """
    __tablename__ = 'pax_abundance_data'

    id  = Column(Integer, primary_key = True)
    abundance = Column(Float)

    dataset_id = Column(Integer, ForeignKey('pax_dataset.id'), index=True)
    pax_dataset = relationship('PaxDataSet', backref='pax_abundance_data')

    subunit_id = Column(Integer, ForeignKey('subunit.id'), index = True)
    subunit = relationship('Subunit', backref = 'pax_abundance_data' )

class Subunit(Base):
    """

    """
    __tablename__ = 'subunit'

    id = Column(Integer, primary_key = True)
    uniprot_id = Column(String(255))
    entrez_id = Column(Integer)
    ensembl_id = Column(String(255), unique = True)
    name = Column(String)
    gene_name = Column(String(255))
    gene_syn  = Column(String(255))
    sequence = Column(String(255), unique = True)
    coefficient = Column(Integer)
    sequence = Column(String(255))
    molecular_weight = Column(Float)

    # Possibly change to subunits
    proteincomplex_id = Column(Integer, ForeignKey('proteincomplex.id'), index = True)
    proteincomplex = relationship('ProteinComplex', backref = 'subunit')


class CellLine(Base):
    """ Represents a CellLine


    """
    __tablename__ = 'cellline'

    id = Column(Integer, primary_key = True)
    name = Column(String)

class ProteinComplex(Base):
    """


    """

    __tablename__ = 'proteincomplex'

    id = Column(Integer, primary_key = True)
    name = Column(String(255))
    go_id = Column(String(255))
    go_dsc = Column(String(255))
    funcat_id = Column(String(255))
    funcat_dsc = Column(String(255))
    su_cmt = Column(String(255))
    complex_cmt = Column(String(255))
    disease_cmt = Column(String(255))

    entry_id = Column(Integer, ForeignKey('entry.id'), index = True)
    entry = relationship('Entry', backref = 'proteincomplex')

class Method(Base):
    """

    """
    __tablename__ = 'method'
    id = Column(Integer, primary_key = True)
    name = Column(String(255))
    comment = Column(String(255))

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


class Synonyms(Base):
    """

    """
    __tablename__ = 'synonym'
    id = Column(Integer, primary_key = True)
    name = Column(String)


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



class CommonSchema(data_source.CachedDataSource):
    """
    A Local SQLlite copy of the aggregation of data_source modules

    NOTES: For now.. Keeping the SQL local, out of cache
    """
    base_model = Base

    def load_content(self):

        session = self.session

        corumdb = corum.Corum(cache_dirname = self.cache_dirname , clear_content = True, load_content=False, download_backup=False, max_entries = 5)
        # corumdb = corum.Corum(name = 'corum', clear_content = True, load_content=False, download_backup=False, max_entries = 5)
        corumdb.load_content()
        self.corum_session = corumdb.session

        paxdb = pax.Pax(cache_dirname = self.cache_dirname, clear_content = True,  load_content=False, download_backup=False, max_entries = 5)
        # paxdb = pax.Pax(name = 'pax', clear_content = True,  load_content=False, download_backup=False, max_entries = 5)
        paxdb.load_content()
        self.pax_session = paxdb.session

        jaspardb = jaspar.Jaspar(cache_dirname = self.cache_dirname, clear_content = True, load_content=False, download_backup=False, max_entries = 5)
        # jaspardb = jaspar.Jaspar(name = 'jaspar', clear_content = True, load_content=False, download_backup=False, max_entries = 5)
        jaspardb.load_content()
        self.jaspar_session = jaspardb.session


        # ecmDB = ecmdb.Ecmdb(cache_dirname = self.cache_dirname, clear_content = True, load_content=False, download_backup=False, max_entries = 5)
        # ecmDB = ecmdb.Ecmdb(name = 'ecmdb', clear_content = True, load_content=False, download_backup=False, max_entries = 20)
        # ecmDB.load_content()
        # self.ecmdb_session = ecmDB.session
        #
        # self.ecmdb_synonym = ecmdb_session.query(ecmdb.Synonym).all()
        # self.ecmdb_compartment = ecmdb_session.query(ecmdb.Compartment).all()
        # self.ecmdb_concentration = ecmdb_session.query(ecmdb.Concentration).all()
        # self.ecmdb_resource = ecmdb_session.query(ecmdb.Resource).all()
        # self.ecmdb_compound = ecmdb_session.query(ecmdb.Compound).all()

        # sabiodb = sabio_rk.SabioRk(cache_dirname = self.cache_dirname, clear_content = True, load_content=False, download_backup=False, max_entries = 5)
        # sabiodb = sabio_rk.SabioRk(name = 'sabio', clear_content = True, load_content=False, download_backup=False, max_entries = 20)
        # sabiodb.load_content()
        # self.sabio_session = sabiodb.session

        # self.sabio_synonym = sabio_session.query(sabio_rk.Synonym).all()
        # self.sabio_resource = sabio_session.query(sabio_rk.Resource).all()
        # self.sabio_entry = sabio_session.query(sabio_rk.Entry).all()
        # self.sabio_reactionparticipant = sabio_session.query(sabio_rk.ReactionParticipant).all()
        # self.sabio_parameter = sabio_session.query(sabio_rk.Parameter).all()
        # self.sabio_kineticlaw = sabio_session.query(sabio_rk.KineticLaw).all()
        # self.sabio_compoundstructure = sabio_session.query(sabio_rk.CompoundStructure).all()
        # self.sabio_enzyme = sabio_session.query(sabio_rk.Enzyme).all()
        # self.sabio_enzymesubunit = sabio_session.query(sabio_rk.EnzymeSubunit).all()
        # self.sabio_compartment = sabio_session.query(sabio_rk.Compartment).all()


        """ NOT WORKING """
        #TODO: AE needs 1. normalization to schema(max_entries implementation/working functions) 2. working tests

        # arrayexpressdb = array_express.ArrayExpress(name = 'array_express', load_content=False, download_backup=False, max_entries = 5)
        # arrayexpressdb.load_experiments_from_text()
        # array_session = arrayexpressdb.session

        self.add_content_to_db()



    def add_content_to_db(self):

        # Pax
        pax_dataset = self.pax_session.query(pax.Dataset).all()
        pax_protein = self.pax_session.query(pax.Protein).all()
        pax_observation = self.pax_session.query(pax.Observation).all()

        for item in pax_dataset:
            dataset = self.get_or_create_object(PaxDataSet, file_name = item.file_name,
                score = item.score, weight = item.weight, coverage = item.coverage)
            dataset.entry = self.get_or_create_object(Entry, name = 'protein abundance')
            dataset.entry.taxon = self.get_or_create_object(Taxon, ncbi_id = item.taxon_ncbi_id,
                species_name = self.pax_session.query(pax.Taxon).get(item.taxon_ncbi_id).species_name)
            dataset.entry.resource = self.get_or_create_object(Resource, namespace = 'url', id = item.publication)
        for item in pax_protein:
            subunit = self.get_or_create_object(Subunit, ensembl_id = item.string_id)
        for item in pax_observation:
            obs = self.get_or_create_object(PaxAbundanceData, abundance = item.abundance,
                pax_dataset = self.session.query(PaxDataSet).get(item.dataset_id))
            ensembl_id = self.pax_session.query(pax.Protein).get(item.protein_id).string_id
            obs.subunit = self.session.query(Subunit).filter_by(ensembl_id = ensembl_id).first()

        # Corum
        corum_complex = self.corum_session.query(corum.Complex).all()
        corum_subunit = self.corum_session.query(corum.Subunit).all()
        corum_observation = self.corum_session.query(corum.Observation).all()

        for item in corum_complex:
            pcomplex = self.get_or_create_object(ProteinComplex, name = item.complex_name, go_id = item.go_id,
                go_dsc = item.go_dsc, funcat_id = item.funcat_id, funcat_dsc = item.funcat_dsc, su_cmt = item.su_cmt,
                complex_cmt = item.complex_cmt, disease_cmt = item.disease_cmt)
            obs = self.corum_session.query(corum.Observation).get(item.observation_id)
            #FIXME: Don't know why this needs to be entered in the function to display
            #everything and why pcomplex.entry.taxon = etc. is not working
            entry = self.get_or_create_object(Entry, name = 'protein complex',
                    taxon = self.get_or_create_object(Taxon, ncbi_id = obs.taxon_ncbi_id,
                    species_name = self.corum_session.query(corum.Taxon).get(obs.taxon_ncbi_id).swissprot_id),
                    resource = self.get_or_create_object(Resource, namespace = 'pubmed', id = obs.pubmed_id),
                    cellline = self.get_or_create_object(CellLine, name = obs.cell_line),
                    method = self.get_or_create_object(Method, name = 'purification', comment = obs.pur_method))
            pcomplex.entry = entry

        for item in corum_subunit:
            subunit = self.get_or_create_object(Subunit, uniprot_id = item.su_uniprot,
                entrez_id = item.su_entrezs, name = item.protein_name, gene_name=item.gene_name,
                gene_syn = item.gene_syn)
            subunit.proteincomplex = self.session.query(ProteinComplex).get(item.complex_id)

        #TODO: ADD Jaspar
        jaspar_transcriptionfactor = self.jaspar_session.query(jaspar.TranscriptionFactor).all()
        jaspar_matrix = self.jaspar_session.query(jaspar.Matrix).all()
        jaspar_matrixposition = self.jaspar_session.query(jaspar.MatrixPosition).all()
        jaspar_collection = self.jaspar_session.query(jaspar.Collection).all()
        jaspar_type = self.jaspar_session.query(jaspar.Type).all()
        jaspar_family = self.jaspar_session.query(jaspar.Family).all()
        jaspar_class = self.jaspar_session.query(jaspar.Class).all()
        jaspar_species = self.jaspar_session.query(jaspar.Species).all()
        jaspar_resource = self.jaspar_session.query(jaspar.Resource).all()
        jaspar_subunit = self.jaspar_session.query(jaspar.Subunit).all()
        jaspar_transcriptionfactor_subunit = self.jaspar_session.query(jaspar.transcription_factor_subunit).all()




        if self.verbose:
            print('Comitting')
        self.session.commit()
