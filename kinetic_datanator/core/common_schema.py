# -*- coding: utf-8 -*-

"""
This code is a common schema for all the kinetic_datanator modules

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-07-31
:Copyright: 2017, Karr Lab
:License: MIT
"""

from sqlalchemy import Column, BigInteger, Integer, Float, String, Text, ForeignKey, Boolean, Table,  Numeric, or_
from sqlalchemy.orm import relationship, backref, sessionmaker
from kinetic_datanator.core import data_source
from kinetic_datanator.data_source import corum, pax, jaspar, jaspar, array_express, ecmdb, sabio_rk, intact, uniprot
import sqlalchemy.ext.declarative
from six import BytesIO
import six
from bioservices import UniProt
from ete3 import NCBITaxa
import pandas as pd
import numpy
import os
import time
import re
# from sqlalchemy import MetaData
# from sqlalchemy_schemadisplay import create_schema_graph


Base = sqlalchemy.ext.declarative.declarative_base()


class Observation(Base):
    """
    Represents an Observation of a Physical Entity or Property in the Common Schema

    Attributes:
        id (:obj:`int`): Common Schema Observation Identifier
        _metadata_id (:obj:`int` of :obj:`Metadata`): Related Metadata ID

    """

    __tablename__ = 'observation'

    id = Column(Integer, primary_key = True)

    _metadata_id = Column(Integer, ForeignKey('_metadata.id'))
    _metadata = relationship('Metadata', backref = 'observation')

    __mapper_args__ = {'polymorphic_identity' : 'observation'}


_metadata_taxon = Table(
    '_metadata_taxon', Base.metadata,
    Column('_metadata_id', Integer, ForeignKey('_metadata.id'), index=True),
    Column('taxon_id', Integer, ForeignKey('taxon.ncbi_id'), index=True),
)
    # :obj:`Table`: Metadata:Taxon many-to-many association table

_metadata_method = Table(
    '_metadata_method', Base.metadata,
    Column('_metadata_id', Integer, ForeignKey('_metadata.id'), index=True),
    Column('method_id', Integer, ForeignKey('method.id'), index=True),
)
    # :obj:`Table`: Metadata:Method many-to-many association table

_metadata_resource = Table(
    '_metadata_resource', Base.metadata,
    Column('_metadata_id', Integer, ForeignKey('_metadata.id'), index=True),
    Column('resource_id', Integer, ForeignKey('resource.id'), index=True),
)
    # :obj:`Table`: Metadata:Resource many-to-many association table

_metadata_cell_line = Table(
    '_metadata_cell_line', Base.metadata,
    Column('_metadata_id', Integer, ForeignKey('_metadata.id'), index=True),
    Column('cell_line_id', Integer, ForeignKey('cell_line.id'), index=True),
)
    # :obj:`Table`: Metadata:CellLine many-to-many association table

_metadata_synonym = Table(
    '_metadata_synonym', Base.metadata,
    Column('_metadata_id', Integer, ForeignKey('_metadata.id'), index=True),
    Column('synonym_id', Integer, ForeignKey('synonym.id'), index=True),
)
    # :obj:`Table`: Metadata:Synonym many-to-many association table

_metadata_conditions = Table(
    '_metadata_conditions', Base.metadata,
    Column('_metadata_id', Integer, ForeignKey('_metadata.id'), index=True),
    Column('conditions_id', Integer, ForeignKey('conditions.id'), index=True),
)
    # :obj:`Table`: Metadata:Conditions many-to-many association table

_metadata_compartment = Table(
    '_metadata_compartment', Base.metadata,
    Column('_metadata_id', Integer, ForeignKey('_metadata.id'), index=True),
    Column('compartment_id', Integer, ForeignKey('cell_compartment.id'), index=True),
)
    # :obj:`Table`: Metadata:Conditions many-to-many association table

subunit_interaction = Table(
    'subunit_interaction', Base.metadata,
    Column('protein_subunit_id', Integer, ForeignKey('protein_subunit.subunit_id'), index = True),
    Column('interaction_id', Integer, ForeignKey('protein_interactions.interaction_id'), index = True)
)



class Metadata(Base):
    """
    Table representing Metadata identifiers for entities and properties

    Attributes:
        name (:obj:`str`): Name of the entity or property

    """
    __tablename__ = '_metadata'

    id = Column(Integer, primary_key = True)
    name = Column(String(255), unique = True)

    taxon = relationship('Taxon', secondary = _metadata_taxon, backref ='_metadata')
    method = relationship('Method', secondary = _metadata_method, backref ='_metadata')
    resource = relationship('Resource', secondary = _metadata_resource, backref ='_metadata')
    cell_line = relationship('CellLine', secondary = _metadata_cell_line, backref ='_metadata')
    synonym = relationship('Synonym', secondary = _metadata_synonym, backref ='_metadata')
    conditions = relationship('Conditions', secondary = _metadata_conditions, backref = '_metadata')
    cell_compartment = relationship('CellCompartment', secondary = _metadata_compartment, backref = '_metadata')

class Method(Base):
    """
    Represents the method of collection for a given entity or Property

    Attributes:
        name (:obj:`str`): Name of the Method
        comments (:obj:`str`): Comments on the method

    """

    id = Column(Integer, primary_key = True)
    name = Column(String(255))
    comments = Column(String(255))

    __tablename__ = 'method'

class Taxon(Base):
    """
    Represents the species of a given physical entity or property

    Attributes:
        ncbi_id (:obj:`int`): NCBI id of the species
        name (:obj:`str`): Name of the species

    """

    ncbi_id = Column(Integer, primary_key = True)
    name = Column(String(255))

    __tablename__ = 'taxon'

class Synonym(Base):
    """
    Represents a synonym of a given physical entity or property

    Attributes:
        name (:obj:`str`): Name of the Synonym

    """

    id = Column(Integer, primary_key = True)
    name = Column(String(255))

    __tablename__ = 'synonym'

class Resource(Base):
    """
    Represents a resource of a given physical entity or property

    Attributes:
        namespace (:obj:`str`): Name of the classifier of the resource (Ex. Pubmed)
        _id (:obj:`str`): Identifier of the resource

    """
    id = Column(Integer, primary_key = True)
    namespace = Column(String(255))
    _id =  Column(String(255))

    __tablename__ = 'resource'

class CellLine(Base):
    """
    Represents a cell line of a given physical entity or property

    Attributes:
        name (:obj:`str`): Name of the Cell Line

    """
    id = Column(Integer, primary_key = True)
    name = Column(String(255))

    __tablename__ = 'cell_line'

class Conditions(Base):
    """
    Represents the conditions of a given physical entity or property

    Attributes:
        growth_status (:obj:`str`): Type of growth status
        media (:obj:`str`): Media composition
        temperature (:obj:`float`): Temperature of the sample (C)
        ph (:obj:`float`): pH of the sample
        growth_system (:obj:`str`): Type of growth system

    """

    id = Column(Integer, primary_key = True)
    growth_status = Column(String(255))
    media = Column(String(255))
    temperature = Column(Float)
    ph = Column(Float)
    growth_system = Column(String(255))

    __tablename__ = 'conditions'



class CellCompartment(Base):
    """
    Represents a cell compartment of a given physical entity or property

    Ties especially to the reacitons because this is where the reactions occur

    Attributes:
        name (:obj:`str`): Name of the Cell Compartment

    """
    __tablename__ = 'cell_compartment'

    id = Column(Integer, primary_key = True)
    name = Column(String(255), unique = True, index=True)

class PhysicalEntity(Observation):
    """
    Represents a Physical Entity in the Common Schema

    Attributes:
        observation_id (:obj:`int`): Common Schema Observation Identifier
        type (:obj:`str`): Type of Physical Entity (Ex. Compound)
        name (:obj:`str`): Name of the Physical Entity (Ex. Water )
    """

    __tablename__ = 'physical_entity'
    __mapper_args__ = {'polymorphic_identity': 'physical_entity'}

    observation_id = Column(Integer, ForeignKey('observation.id'), primary_key = True)
    type = Column(String(255))
    name = Column(String(255))

class ProteinSubunit(PhysicalEntity):
    """
    Represents a Protein Subunit - An instance of Physical Entity

    Attributes:
        subunit_id (:obj:`int`): Common Schema Observation Identifier
        subunit_name (:obj:`str`): Name of the Protein Subunit
        uniprot_id  (:obj:`str`): Uniprot ID for the Subunit
        entrez_id (:obj:`int`): Entrez ID for the Subunit
        gene_name (:obj:`str`): Name of Gene which subunit is derived
        gene_syn (:obj:`str`): Synonym of Gene
        class_name (:obj:`str`): Class Name of Subunit
        family_name (:obj:`str`): Family Name of Subunit
        coefficient (:obj:`str`): Number of units required for complex
        sequence (:obj:`str`): Sequence of Subunit
        molecular_weight (:obj:`float`): Molecular weight of subunit
    """

    __tablename__ = 'protein_subunit'
    __mapper_args__ = {'polymorphic_identity': 'protein_subunit'}

    subunit_id = Column(Integer, ForeignKey('physical_entity.observation_id'), primary_key = True)
    subunit_name = Column(String(255))
    uniprot_id = Column(String(255))
    entrez_id = Column(Integer)
    gene_name = Column(String(255))
    gene_syn  = Column(String(255))
    class_name = Column(String(255))
    family_name = Column(String(255))
    coefficient = Column(Integer)
    canonical_sequence = Column(String(255))
    mass = Column(Integer)
    length = Column(Integer)
    molecular_weight = Column(Float)
    pax_load = Column(Integer)

    interaction = relationship('ProteinInteractions', secondary = subunit_interaction , backref = 'protein_subunit')

    proteincomplex_id = Column(Integer, ForeignKey('protein_complex.complex_id'))
    proteincomplex = relationship('ProteinComplex', backref = 'protein_subunit', foreign_keys=[proteincomplex_id])

class ProteinComplex(PhysicalEntity):
    """
    Represents a Protein Complex - An instance of Physical Entity

    Attributes:
        complex_id (:obj:`int`): Common Schema Observation Identifier
        complex_name (:obj:`str`): Name of the Protein Complex
        go_id (:obj:`str`): GO functional annotation
        go_dsc (:obj:`str`): Description of the annotation
        funcat_id (:obj:`str`): FUNCAT functional annotation
        funcat_dsc (:obj:`str`): Description of the annotation
        su_cmt (:obj:`str`): Subunit comments
        complex_cmt (:obj:`str`): Compex comments
        disease_cmt (:obj:`str`): Disease comments
        class_name (:obj:`str`): Class Name of Subunit
        family_name (:obj:`str`): Family Name of Subunit
        molecular_weight (:obj:`float`): Molecular weight of subunit
    """
    __tablename__ = 'protein_complex'
    __mapper_args__ = {'polymorphic_identity': 'protein_complex'}

    complex_id = Column(Integer, ForeignKey('physical_entity.observation_id'), primary_key = True)
    complex_name = Column(String(255))
    go_id = Column(String(255))
    go_dsc = Column(String(255))
    funcat_id = Column(String(255))
    funcat_dsc = Column(String(255))
    su_cmt = Column(String(255))
    complex_cmt = Column(String(255))
    disease_cmt = Column(String(255))
    class_name = Column(String(255))
    family_name = Column(String(255))
    molecular_weight = Column(Float)


class Compound(PhysicalEntity):
    """
    Represents a Compound - An instance of Physical Entity

    Attributes:
        compound_id (:obj:`int`): Common Schema Observation Identifier
        compound_name (:obj:`str`): Name of the Compound
        description (:obj:`str`):
        comment = Column(String(255))
        _is_name_ambiguous = Column(Boolean)

    """
    __tablename__ = 'compound'
    __mapper_args__ = {'polymorphic_identity': 'compound'}

    compound_id = Column(Integer, ForeignKey('physical_entity.observation_id'), primary_key = True)
    compound_name = Column(String(255))
    description = Column(String(255))
    comment = Column(String(255))
    _is_name_ambiguous = Column(Boolean)

    structure_id = Column(Integer, ForeignKey('structure.struct_id'))
    structure = relationship('Structure', backref = 'compound')


class PhysicalProperty(Observation):
    """
    Represents a Physical Property in the Common Schema

    Attributes:
        observation_id (:obj:`int`): Common Schema Observation Identifier
        type (:obj:`str`): Type of Physical Property (Ex. Concentration)
        name (:obj:`str`): Name of the Physical Property
    """
    observation_id = Column(Integer, ForeignKey('observation.id'), primary_key = True)
    type = Column(String(255))
    name = Column(String(255))

    __tablename__ = 'physical_property'
    __mapper_args__ = {'polymorphic_identity': 'physical_property'}

class Structure(PhysicalProperty):
    """
    Represents a structure of a compound

    Attributes:
        _value_smiles (:obj:`str`): Smiles format for compound representation
        _value_inchi (:obj:`str`): Inchi format for compound representation
        _structure_formula_connectivity (:obj:`str`): Connectivity of compound

    """

    __tablename__ = 'structure'
    __mapper_args__ = {'polymorphic_identity': 'structure'}

    struct_id = Column(Integer, ForeignKey('physical_property.observation_id'), primary_key = True)
    _value_smiles= Column(String(255))
    _value_inchi = Column(String(255))
    _structure_formula_connectivity = Column(String(255))

class Concentration(PhysicalProperty):
    """
    Represents the concentration of an entity

    Attributes:
        value (:obj:`float`): concentration of a tagged compound
        error (:obj:`float`): uncertainty of corresponding concentration value
    """


    __tablename__ = 'concentration'
    __mapper_args__ = {'polymorphic_identity': 'concentration'}

    concentration_id = Column(Integer, ForeignKey('physical_property.observation_id'), primary_key = True)

    compound_id = Column(Integer, ForeignKey('compound.compound_id'))
    compound = relationship('Compound', backref = 'concentration')

    value = Column(Float)
    error = Column(Float)

class KineticLaw(PhysicalProperty):
    """
    Represents the concentration of an entity

    Attributes:
        enzyme_id (:obj:`int`): ID of enzyme driving the kinetic law
            enzyme_type (:obj:`str`): Enzyme classification (Ex. Modifier-Catalyst)
            tissue (:obj:`str`): Tissue from which kinetic law stems from
            mechanism (:obj:`str`): Rate kinetics of Kinetic Law
            equation (:obj:`str`): Equation of the rate kinetics
    """

    __tablename__ = 'kinetic_law'
    __mapper_args__ = {'polymorphic_identity': 'kinetic_law'}

    kineticlaw_id = Column(Integer, ForeignKey('physical_property.observation_id'), primary_key = True)

    enzyme_id = Column(Integer, ForeignKey('protein_complex.complex_id'), index=True)
    enzyme = relationship(ProteinComplex, backref = 'kinetic_law')

    enzyme_type = Column(String(255))
    tissue = Column(String(255))
    mechanism = Column(String(255))
    equation = Column(String(255))

class Reaction(Base):
    """
    Represents a reaction

    Attributes:
        compound_id (:obj:`int`): ID of the corresponding compound
        coefficient (:obj:`float`): Stoichiometric coefficient
        _is_reactant (:obj:`bool`): Indicates of corresponding compound is a reactant
        _is_product (:obj:`bool`): Indicates of corresponding compound is a product
        _is_modifier (:obj:`bool`): Indicates of corresponding compound is a modifier
        rxn_type (:obj:`str`): Classifer of reaction

    """

    __tablename__ = 'reaction'

    reaction_id = Column(Integer, primary_key = True)
    compound_id = Column(Integer, ForeignKey('compound.compound_id'))
    compound = relationship(Compound, backref = 'reaction' )
    compartment_id = Column(Integer, ForeignKey('cell_compartment.id'))
    compartment = relationship(CellCompartment, backref = 'reaction')
    coefficient = Column(Float)
    _is_reactant = Column(Boolean)
    _is_product = Column(Boolean)
    _is_modifier = Column(Boolean)
    rxn_type = Column(String(255))


    kinetic_law_id = Column(Integer, ForeignKey('kinetic_law.kineticlaw_id'))


class AbundanceDataSet(PhysicalProperty):
    """
    Represents a dataset for protein abundance

    Attributes:
        file_name (:obj:`str`): Name of data set which data stems from
        score (:obj:`float`): Quality of the experimental analysis (PAXdb Measure)
        weight (:obj:`int`):
        coverage (:obj:`int`): Percentage of genome coverage
    """

    __tablename__ = 'abundance_dataset'
    __mapper_args__ = {'polymorphic_identity': 'abundance_dataset'}

    dataset_id = Column(Integer, ForeignKey('physical_property.observation_id'), primary_key = True)
    file_name = Column(String, unique = True)
    score = Column(Float)
    weight = Column(Integer)
    coverage = Column(Integer)



class DNABindingDataset(PhysicalProperty):
    """
    Represents a dataset for Transcription Factor Binding

    Attributes:
        version (:obj:`int`): Represents the version of binding matrix
        complex_id (:obj:`int`): Relation ID for transcription factor complex
        subunit_id (:obj:`int`):  Relation ID for transcription factor subunit
    """
    __tablename__ = 'dna_binding_dataset'
    __mapper_args__ = {'polymorphic_identity': 'dna_binding_dataset'}

    dataset_id = Column(Integer, ForeignKey('physical_property.observation_id'), primary_key = True)
    version = Column(Integer)

    complex_id = Column(Integer, ForeignKey('protein_complex.complex_id'))
    tf = relationship('ProteinComplex', backref = 'dna_binding_dataset')

    subunit_id = Column(Integer, ForeignKey('protein_subunit.subunit_id'))
    subunit = relationship('ProteinSubunit', backref = 'dna_binding_dataset')

class Parameter(Base):
    """
    Represents a parameter for a given kinetic law and compound

    Attributes:
        kinetic_law_id (:obj:`int`): corresponding kinetic law to the parameter
        sabio_type (:obj:`int`): sabio identifier for type of parameter
        compound_id (:obj:`int`): corresponding compound for the parameter
        value (:obj:`float`): Value of parameter
        error (:obj:`float`): Uncertainty of the value
        units (:obj:`str`): units of the parameter
        observed_name (:obj:`str`): observed name of parameter
        observed_sabio_type (:obj:`int`): observed sabio type
        observed_value (:obj:`float`): observed value of parameter
        observed_error (:obj:`float`): observed error of parameter
        observed_units (:obj:`str`): observed units of parameter

    """

    __tablename__ = 'parameter'

    parameter_id = Column(Integer, primary_key = True)

    kinetic_law_id = Column(Integer, ForeignKey('kinetic_law.kineticlaw_id') )
    kinetic_law = relationship(KineticLaw, backref = 'parameter')

    sabio_type = Column(Integer, index=True)
    compound_id = Column(Integer, ForeignKey('compound.compound_id'), index=True)
    compound = relationship(Compound, backref= 'parameter')

    value = Column(Float)
    error = Column(Float)
    units = Column(String(255), index=True)

    observed_name = Column(String(255))
    observed_sabio_type = Column(Integer)
    observed_value = Column(Float)
    observed_error = Column(Float)
    observed_units = Column(String(255))

class AbundanceData(Base):
    """
    Represents protein abundance data from the Pax DB database

    Attributes:
        abundance (:obj:`float`): Represents protein abundance from given observation in ppm
        dataset_id  (:obj:`int`): Represents the dataset from which the abundance stems from
        subunit_id  (:obj:`int`): Represents the protein frmo which the abundance stems from
    """

    __tablename__ = 'abundance_data'

    abundance_id  = Column(Integer, primary_key = True)
    abundance = Column(Float)

    dataset_id  = Column(Integer, ForeignKey('abundance_dataset.dataset_id'))
    dataset = relationship('AbundanceDataSet', backref = 'abundance_data', foreign_keys=[dataset_id])

    subunit_id = Column(Integer, ForeignKey('protein_subunit.subunit_id'), index=True)
    subunit = relationship('ProteinSubunit', backref = 'pax_abundance_data' )

    pax_load = Column(Integer)
    uniprot_id = Column(Integer)

class DNABindingData(Base):
    """
    Represents Matrix binding profile for a protein transcription factor

    Attributes:
        position (:obj:`int`): Position in the sequence
        frequency_a (:obj:`int`): Frequency of A
        frequency_c (:obj:`int`): Frequency of C
        frequency_g (:obj:`int`): Frequency of G
        frequency_t (:obj:`int`): Frequency of T
        jaspar_id (:obj:`int`): ID of Jaspar Matrix (used for bulk insert mapping)
        dataset_id  (:obj:`int`): Represents the dataset from which the data stems from
    """
    __tablename__ = 'dna_bidning_data'

    position_id = Column(Integer, primary_key = True)
    position = Column(Integer, index=True)
    frequency_a = Column(Integer)
    frequency_c = Column(Integer)
    frequency_g = Column(Integer)
    frequency_t = Column(Integer)
    jaspar_id = Column(Integer)

    dataset_id  = Column(Integer, ForeignKey('dna_binding_dataset.dataset_id'))
    dataset = relationship('DNABindingDataset', backref = 'dna_binding_data', foreign_keys=[dataset_id])

class ProteinInteractions(PhysicalProperty):
    """
    Represents a protein-protein interaction

    Attributes:
        participant_a (:obj:`str`): Participant A in the interaction
        participant_b (:obj:`str`): Participant B in the interaction
        interaction (:obj:`str`): Interaction ID
        site_a (:obj:`str`): Binding Site of Participant A
        site_b (:obj:`str`): Binding Site of Participant B
        stoich_a (:obj:`str`): Stoichiometry of Participant A
        stoich_b (:obj:`str`): Stoichiometry of Participant B

    """
    __tablename__ = 'protein_interactions'
    __mapper_args__ = {'polymorphic_identity': 'protein_interactions'}

    interaction_id = Column(Integer, ForeignKey('physical_property.observation_id'), primary_key = True)
    participant_a = Column(String(255))
    participant_b = Column(String(255))
    interaction = Column(String(255))
    site_a = Column(String(255))
    site_b = Column(String(255))
    stoich_a = Column(String(255))
    stoich_b = Column(String(255))
    interaction_type = Column(String(255))
    publication = Column(String(255))


class Progress(Base):
    """
    Represents amount loaded of large DBs (Ex. Pax and Sabio)
    Attributes:
        database_name (:obj:`str`): Name of observed databse
        amount_loaded (:obj:`int`): Amount of entries loaded in Common Schema

    """
    __tablename__ = 'progress'

    database_name = Column(String, primary_key = True)
    amount_loaded = Column(Integer)


class CommonSchema(data_source.HttpDataSource):
    """
    A Local SQLlite copy of the aggregation of data_source modules

    """
    base_model = Base

    def __init__(self, name=None, cache_dirname=None, clear_content=False, load_content=False, max_entries=float('inf'),
                 commit_intermediate_results=False, download_backup=True, verbose=False, load_entire_small_DBs = False,
                  clear_requests_cache=False, download_request_backup=False):

        """
        Args:
            name (:obj:`str`, optional): name
            cache_dirname (:obj:`str`, optional): directory to store the local copy of the data source and the HTTP requests cache
            clear_content (:obj:`bool`, optional): if :obj:`True`, clear the content of the sqlite local copy of the data source
            load_content (:obj:`bool`, optional): if :obj:`True`, load the content of the local sqlite database from the external source
            max_entries (:obj:`float`, optional): maximum number of entries to save locally
            commit_intermediate_results (:obj:`bool`, optional): if :obj:`True`, commit the changes throughout the loading
                process. This is particularly helpful for restarting this method when webservices go offline.
            download_backup (:obj:`bool`, optional): if :obj:`True`, load the local copy of the data source from the Karr Lab server
            verbose (:obj:`bool`, optional): if :obj:`True`, print status information to the standard output
        """

        super(CommonSchema, self).__init__(name=name, cache_dirname=cache_dirname, clear_content=clear_content,
                                      load_content=False, max_entries=max_entries,
                                      commit_intermediate_results=commit_intermediate_results,
                                      download_backup=download_backup, verbose=verbose,
                                      clear_requests_cache=clear_requests_cache, download_request_backup=download_request_backup)

        self.load_entire_small_DBs = load_entire_small_DBs

        if download_backup and load_content:
            self.pax_loaded = self.session.query(Progress).filter_by(database_name = 'Pax').first().amount_loaded
            self.sabio_loaded = self.session.query(Progress).filter_by(database_name = 'Sabio').first().amount_loaded
            self.intact_loaded = self.session.query(Progress).filter_by(database_name = 'IntAct').first().amount_loaded
            self.load_small_db_switch = False
            self.session.query(Progress).delete()
            self.load_content()
        elif load_content:
            self.pax_loaded = 0
            self.sabio_loaded = 0
            self.intact_loaded = 0
            self.load_small_db_switch = True
            self.load_content()

    def load_content(self):
        ## Initiate Observation and direct Subclasses
        observation = Observation()
        observation.physical_entity = PhysicalEntity()
        self.entity = observation.physical_entity
        observation.physical_property = PhysicalProperty()
        self.property = observation.physical_property

        # Add all DBs
        self.add_intact()
        if self.verbose:
            print('IntAct Done')
        self.add_paxdb()
        if self.verbose:
            print('Pax Done')
        if self.load_small_db_switch:
            self.add_corumdb()
            if self.verbose:
                print('Corum Done')
            self.add_jaspardb()
            if self.verbose:
                print('Jaspar Done')
            self.add_ecmdb()
            if self.verbose:
                print('ECMDB Done')
        self.add_sabiodb()
        if self.verbose:
            print('Sabio Done')


        ## Add missing subunit information
        # self.add_uniprot()
        # if self.verbose:
        #     print('Uniprot Done')

        ## Add missing Taxon information
        self.fill_missing_ncbi_names()

    def create_schema_png(self):
        if switch:
            # create the pydot graph object by autoloading all tables via a bound metadata object
            graph = create_schema_graph(metadata=MetaData(self.engine),
               show_datatypes=False, # The image would get nasty big if we'd show the datatypes
               show_indexes=False, # ditto for indexes
               rankdir='TB', # From left to right (instead of top to bottom)
               concentrate=False # Don't try to join the relation lines together
            )
            graph.write_png(os.getcwd())

    # def fill_missing_subunit_info(self):
    #     t0 = time.time()
    #     u = UniProt(verbose = False)
    #
    #
    #     subunits = self.session.query(ProteinSubunit).filter_by(entrez_id = None)
    #     while(subunits.count() != 0):
    #         initial_count = subunits.count()
    #         chunk = 1000
    #         uni = []
    #         for items in subunits.limit(chunk).all():
    #             uni.append(str(items.uniprot_id))
    #         entrez_dict = u.mapping(fr = 'ACC', to = 'P_ENTREZGENEID', query = uni)
    #         for protein in subunits.limit(chunk).all():
    #             if protein.entrez_id == None and protein.uniprot_id in entrez_dict.keys():
    #                 protein.entrez_id = int(entrez_dict[protein.uniprot_id][0])
    #         if subunits.count() == initial_count:
    #             break
    #
    #     self.session.commit()
    #
    #     subunits = self.session.query(ProteinSubunit).filter_by(gene_name = None)
    #     while(subunits.count() != 0):
    #         initial_count = subunits.count()
    #         chunk = 200
    #         uni = []
    #         for items in subunits.limit(chunk).all():
    #             uni.append(str(items.uniprot_id))
    #         df = u.get_df(uni, nChunk = 200)
    #         df.set_index('Entry', inplace = True)
    #         for protein in subunits.limit(chunk).all():
    #             if protein.uniprot_id in df.index:
    #                 protein.subunit_name = str(df.loc[protein.uniprot_id,'Protein names'])
    #                 protein.gene_name = str(df.loc[protein.uniprot_id,'Gene names']).replace('[', '').replace(']','')
    #         if subunits.count() == initial_count:
    #             break
    #
    #     self.session.commit()
    #
    #
    #     subunits = self.session.query(ProteinSubunit).filter_by(canonical_sequence = None)
    #     while(subunits.count() != 0):
    #         initial_count = subunits.count()
    #         chunk = 1000
    #         uni = []
    #         for items in subunits.limit(chunk).all():
    #             uni.append(str(items.uniprot_id))
    #         df = u.get_df(uni, nChunk = 200)
    #         df.set_index('Entry', inplace = True)
    #         for protein in subunits.limit(chunk).all():
    #             if protein.uniprot_id in df.index:
    #                 protein.canonical_sequence = str(df.loc[protein.uniprot_id,'Sequence'])
    #                 protein.mass = str(df.loc[protein.uniprot_id,'Mass'])
    #                 if type(df.loc[protein.uniprot_id,'Length']) == numpy.int64:
    #                     protein.length = int(df.loc[protein.uniprot_id,'Length'])
    #                 else:
    #                     protein.length = int(df.loc[protein.uniprot_id,'Length'].iloc[0])
    #         if subunits.count() == initial_count:
    #             break

        if self.verbose:
            print('Total time taken for Uniprot fillings: ' + str(time.time()-t0) + ' secs')

        if self.verbose:
            print('Comitting')
        self.session.commit()

    def fill_missing_ncbi_names(self):
        t0 = time.time()

        ncbi = NCBITaxa()

        ncbi_ids = []
        species = self.session.query(Taxon).all()
        for items in species:
            ncbi_ids.append(items.ncbi_id)

        species_dict = ncbi.get_taxid_translator(ncbi_ids)

        for tax in species:
            if tax.ncbi_id in species_dict.keys():
                tax.name = species_dict[tax.ncbi_id]

        if self.verbose:
            print('Total time taken for NCBI fillings: ' + str(time.time()-t0) + ' secs')

        if self.verbose:
            print('Comitting')
        self.session.commit()

    def add_paxdb(self):
        """
        Adds Pax Database from Pax.sqlite file in Karr Server

        """

        t0 = time.time()
        paxdb = pax.Pax(cache_dirname = self.cache_dirname, verbose = self.verbose)
        pax_ses = paxdb.session
        u = UniProt(verbose = False)

        _entity = self.entity
        _property = self.property


        if self.max_entries == float('inf'):
            pax_dataset = pax_ses.query(pax.Dataset).all()
        else:
            pax_dataset = pax_ses.query(pax.Dataset).filter(pax.Dataset.id.in_\
                (range(self.pax_loaded+1, self.pax_loaded+1 + int(self.max_entries/5))))

        for dataset in pax_dataset:
            metadata = self.get_or_create_object(Metadata, name = dataset.file_name)
            metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = dataset.taxon_ncbi_id))
            metadata.resource.append(self.get_or_create_object(Resource, namespace = 'url', _id = dataset.publication))
            _property.abundance_dataset = self.get_or_create_object(AbundanceDataSet, type = 'Protein Abundance Dataset',
                name = dataset.file_name, file_name = dataset.file_name, score = dataset.score, weight = dataset.weight,
                coverage= dataset.coverage, _metadata = metadata)
            abundance = pax_ses.query(pax.Observation).filter_by(dataset_id = dataset.id).all()
            uni = [str(d.protein.uniprot_id) for d in abundance]

            self.session.bulk_insert_mappings(AbundanceData,
                [
                    dict(abundance = data.abundance, pax_load = data.dataset_id, \
                        uniprot_id = data.protein.uniprot_id)\
                        for data in abundance
                ])

            self.session.bulk_insert_mappings(ProteinSubunit,
                [
                    dict(uniprot_id = data.protein.uniprot_id, type = 'Protein Subunit',
                    pax_load = data.dataset_id) for data in abundance
                ])

            self.session.commit()

            for subunit in self.session.query(ProteinSubunit).filter_by(pax_load = dataset.id).all():
                subunit._metadata = metadata

            for rows in self.session.query(AbundanceData).filter_by(pax_load = dataset.id).all():
                rows.dataset = _property.abundance_dataset
                rows.subunit = self.session.query(ProteinSubunit).filter_by(pax_load = dataset.id).filter_by(uniprot_id = rows.uniprot_id).first()

        self.get_or_create_object(Progress, database_name = 'Pax', amount_loaded = self.pax_loaded + (self.max_entries/5))

        if self.verbose:
            print('Total time taken for Pax: ' + str(time.time()-t0) + ' secs')

        if self.verbose:
            print('Comitting')
        self.session.commit()

    def add_corumdb(self):
        """
        Adds Corum Database from Corum.sqlite file in Karr Server

        """
        t0 = time.time()
        corumdb = corum.Corum(cache_dirname = self.cache_dirname, verbose = self.verbose)
        corum_ses = corumdb.session

        _entity = self.entity
        _property = self.property

        corum_complex = corum_ses.query(corum.Complex).all()
        corum_subunit = corum_ses.query(corum.Subunit).all()

        max_entries = self.max_entries

        if self.load_entire_small_DBs:
            max_entries = float('inf')

        entries = 0
        for complx in corum_complex:
            if entries < max_entries:
                entry = complx.observation
                metadata = self.get_or_create_object(Metadata, name = complx.complex_name)
                metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = entry.taxon_ncbi_id))
                metadata.resource.append(self.get_or_create_object(Resource, namespace = 'pubmed', _id = str(entry.pubmed_id)))
                metadata.method.append(self.get_or_create_object(Method, name = 'purification', comments = entry.pur_method))
                metadata.cell_line.append(self.get_or_create_object(CellLine, name = entry.cell_line))
                _entity.protein_complex = self.get_or_create_object(ProteinComplex,
                    type = 'Protein Complex', name = complx.complex_name, complex_name = complx.complex_name, go_id = complx.go_id,
                    go_dsc = complx.go_dsc, funcat_id = complx.funcat_id, funcat_dsc = complx.funcat_dsc, su_cmt = complx.su_cmt,
                    complex_cmt = complx.complex_cmt, disease_cmt = complx.disease_cmt,  _metadata = metadata)
                entries += 1

        entries = 0
        for subunit in corum_subunit:
            if entries < max_entries:
                complx = self.session.query(ProteinComplex).filter_by(complex_name = subunit.complex.complex_name).first()
                entry = self.session.query(Observation).get(complx.complex_id)
                _entity.protein_subunit = self.get_or_create_object(ProteinSubunit,
                    type = 'Protein Subunit', uniprot_id = subunit.su_uniprot,
                    entrez_id = subunit.su_entrezs, name = subunit.protein_name, subunit_name = subunit.protein_name, gene_name=subunit.gene_name,
                    gene_syn = subunit.gene_syn, proteincomplex = complx, _metadata = self.session.query(Metadata).get(entry._metadata_id))
                entries += 1
        if self.verbose:
            print('Total time taken for Corum: ' + str(time.time()-t0) + ' secs')

        if self.verbose:
            print('Comitting')
        self.session.commit()

    def add_jaspardb(self):
        t0 = time.time()
        jaspardb = jaspar.Jaspar(cache_dirname = self.cache_dirname,verbose = self.verbose)

        jasp_ses = jaspardb.session

        _entity = self.entity
        _property = self.property

        def list_to_string(list_):
            if list_:
                ans = ''
                for word in list_:
                    ans += word
                return ans
            else:
                return ''

        matrix = jasp_ses.query(jaspar.Matrix).all()

        max_entries = self.max_entries

        if self.load_entire_small_DBs:
            max_entries = float('inf')

        entries = 0
        for entry in matrix:
            if entries <= max_entries:
                annotations = jasp_ses.query(jaspar.Annotation).filter_by(ID = entry.ID)
                class_ = [c.VAL for c in annotations.filter_by(TAG = 'class').all()]
                family_ = [f.VAL for f  in annotations.filter_by(TAG = 'family').all()]
                pubmed = [p.VAL for p in annotations.filter_by(TAG = 'medline').all()]
                type_ = [t.VAL for t in annotations.filter_by(TAG = 'type').all()]
                species = [s.TAX_ID for s in jasp_ses.query(jaspar.Species).filter_by(ID = entry.ID).all()]
                protein = [p.ACC for p in jasp_ses.query(jaspar.Protein).filter_by(ID = entry.ID).all()]

                metadata = self.get_or_create_object(Metadata, name = entry.NAME + ' Binding Motif')
                metadata.method.append(self.get_or_create_object(Method, name = list_to_string(type_)))
                metadata.resource = [self.get_or_create_object(Resource, namespace = 'pubmed', _id = ref) for ref in pubmed]
                metadata.taxon = [self.get_or_create_object(Taxon, ncbi_id = int(tax)) for tax in species if tax != '-' ]
                if '::' in entry.NAME:
                    _entity.protein_complex = self.get_or_create_object(ProteinComplex, type = 'Transcription Factor Complex',
                    name = entry.NAME, complex_name = entry.NAME, complex_cmt = 'transcription factor', class_name = list_to_string(class_),
                    family_name = list_to_string(family_), _metadata = metadata)
                    _property.dna_binding_dataset = self.get_or_create_object(DNABindingDataset, type = 'DNA Binding Dataset',
                    name = entry.NAME, version = entry.VERSION, tf= _entity.protein_complex, _metadata = metadata)
                else:
                    prot = protein[0] if protein else None
                    _entity.protein_subunit = self.get_or_create_object(ProteinSubunit, uniprot_id = prot,
                    type = 'Transcription Factor Subunit', name = entry.NAME, subunit_name = entry.NAME, gene_name = entry.NAME,
                    class_name = list_to_string(class_), family_name = list_to_string(family_), _metadata = metadata)
                    _property.dna_binding_dataset = self.get_or_create_object(DNABindingDataset, type = 'DNA Binding Dataset',
                    name = entry.NAME, version = entry.VERSION, subunit = _entity.protein_subunit, _metadata = metadata)
                subquery = jasp_ses.query(jaspar.Data).filter_by(ID = entry.ID)
                for position in range(1,1+max(set([c.col for c in subquery.all()]))):
                    freq = subquery.filter_by(col = position)
                    self.get_or_create_object(DNABindingData, position = position, frequency_a = freq.filter_by(row = 'A').first().val,
                    frequency_c = freq.filter_by(row = 'C').first().val, frequency_t = freq.filter_by(row = 'T').first().val,
                    frequency_g = freq.filter_by(row = 'G').first().val, jaspar_id = entry.ID, dataset = _property.dna_binding_dataset)
            entries += 1

        if self.verbose:
            print('Total time taken for Jaspar: ' + str(time.time()-t0) + ' secs')

        if self.verbose:
            print('Comitting')
        self.session.commit()

    #def add_jaspardb(self):
        # """ DEPRECATED"""
        # t0 = time.time()
        # jaspardb = jaspar.Jaspar(cache_dirname = self.cache_dirname, clear_content = False,
        #     load_content= False, download_backup= True, verbose = self.verbose)
        # jasp_ses = jaspardb.session
        #
        # _entity = self.entity
        # _property = self.property
        #
        # jaspar_matrix = jasp_ses.query(jaspar.Matrix).all()
        # position = jasp_ses.query(jaspar.MatrixPosition).all()
        #
        # self.session.bulk_insert_mappings(DNABindingData,
        #     [
        #         dict(position = pos.position,
        #         frequency_a = pos.frequency_a, frequency_c = pos.frequency_c, frequency_g = pos.frequency_g,
        #         frequency_t = pos.frequency_t, jaspar_id = pos.matrix_id) for pos in position
        #     ])
        #
        # max_entries = self.max_entries
        #
        # if self.load_entire_small_DBs:
        #     max_entries = float('inf')
        #
        # entries = 0
        # for item in jaspar_matrix:
        #     if entries < max_entries:
        #         tf = item.transcription_factor
        #         if item.type_id:
        #             type_name = item.type.name
        #         if item.transcription_factor.classes:
        #             class_name = item.transcription_factor.classes[0].name
        #         if item.transcription_factor.families:
        #             fam_name = item.transcription_factor.families[0].name
        #         metadata = self.get_or_create_object(Metadata, name = tf.name)
        #         for speice in item.transcription_factor.species:
        #             metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = speice.ncbi_id))
        #         metadata.method.append(self.get_or_create_object(Method, name = type_name))
        #         for docs in item.references:
        #             metadata.resource.append(self.get_or_create_object(Resource, namespace = 'pubmed', _id = docs.pubmed_id))
        #         if '::' in tf.name:
        #             su = tf.name.split('::')
        #             dialoge = 'Subunits are: '
        #             for items in su:
        #                 dialoge += (items + ' ')
        #             _entity.protein_complex = self.get_or_create_object(ProteinComplex, type = 'Transcription Factor Complex', name = tf.name,
        #                 complex_name = tf.name, su_cmt = dialoge, complex_cmt = 'transcription factor', class_name = class_name, family_name = fam_name, _metadata = metadata)
        #             _property.dna_binding_dataset = self.get_or_create_object(DNABindingDataset, type = 'DNA Binding Dataset', name = tf.name,
        #                 version = item.version, tf = _entity.protein_complex, _metadata = metadata)
        #             for row in self.session.query(DNABindingData).filter_by(jaspar_id = item.id).all():
        #                 row.dataset = _property.dna_binding_dataset
        #         else:
        #             if tf.subunits:
        #                 uniprot_id = tf.subunits[0].uniprot_id
        #             _entity.protein_subunit = self.get_or_create_object(ProteinSubunit, uniprot_id = uniprot_id, type = 'Transcription Factor Subunit',
        #                 name = tf.name, subunit_name = tf.name, gene_name = tf.name,
        #                 class_name = class_name, family_name = fam_name, _metadata = metadata)
        #             _property.dna_binding_dataset = self.get_or_create_object(DNABindingDataset, type = 'DNA Binding Dataset', name = tf.name,
        #                 version = item.version, subunit = _entity.protein_subunit, _metadata = metadata)
        #             for row in self.session.query(DNABindingData).filter_by(jaspar_id = item.id).all():
        #                 row.dataset = _property.dna_binding_dataset
        #         entries += 1
        #
        # if self.verbose:
        #     print('Total time taken for Jaspar: ' + str(time.time()-t0) + ' secs')
        #
        # if self.verbose:
        #     print('Comitting')
        # self.session.commit()

    def add_ecmdb(self):
        """
        Adds ECMDB Database from Ecmdb.sqlite file in Karr Server

        """
        t0 = time.time()
        ecmDB = ecmdb.Ecmdb(cache_dirname = self.cache_dirname, verbose = self.verbose)
        ecm_ses = ecmDB.session

        _entity = self.entity
        _property = self.property

        ecmdb_compound = ecm_ses.query(ecmdb.Compound).all()

        max_entries = self.max_entries

        if self.load_entire_small_DBs:
            max_entries = float('inf')

        entries = 0
        for compound in ecmdb_compound:
            if entries < max_entries:
                concentration = compound.concentrations
                ref = compound.cross_references
                compart = compound.compartments
                syn = compound.synonyms
                metadata = self.get_or_create_object(Metadata, name = compound.name)
                metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = 562, name = 'E.Coli'))
                metadata.resource = [self.get_or_create_object(Resource, namespace = docs.namespace, _id = docs.id) for docs in ref]
                metadata.cell_compartment = [self.get_or_create_object(CellCompartment, name = areas.name) for areas in compart]
                metadata.synonym = [self.get_or_create_object(Synonym, name = syns.name) for syns in syn]
                _property.structure = self.get_or_create_object(Structure, type = 'Structure', name = compound.name,
                    _value_inchi = compound.structure,
                   _structure_formula_connectivity = compound._structure_formula_connectivity, _metadata = metadata)
                _entity.compound = self.get_or_create_object(Compound, type = 'Compound', name = compound.name,
                    compound_name = compound.name, description = compound.description, comment = compound.comment, structure = _property.structure,
                    _metadata = metadata)
                index = 0 if concentration else -1
                if index == -1: continue
                for rows in concentration:
                    new_metadata = self.get_or_create_object(Metadata, name = compound.name+ ' Concentration ' + str(index))
                    new_metadata.taxon = metadata.taxon
                    new_metadata.cell_compartment = metadata.cell_compartment
                    new_metadata.synonym = metadata.synonym
                    new_metadata.cell_line.append(self.get_or_create_object(CellLine, name = rows.strain))
                    new_metadata.conditions.append(self.get_or_create_object(Conditions, growth_status = rows.growth_status,
                        media = rows.media, temperature = rows.temperature, growth_system = rows.growth_system))
                    new_metadata.resource = [self.get_or_create_object(Resource, namespace = docs.namespace, _id = docs.id) for docs in rows.references]
                    _property.concentration = self.get_or_create_object(Concentration, type = 'Concentration', name = compound.name+ ' Concentration '+str(index),
                        value = rows.value, error = rows.error, _metadata = new_metadata, compound = _entity.compound)
                    index += 1
                entries += 1

        if self.verbose:
            print('Total time taken for ECMDB: ' + str(time.time()-t0) + ' secs')

        if self.verbose:
            print('Comitting')
        self.session.commit()

    def add_sabiodb(self):
        """
        Adds Sabio Database from sabio.sqlite file in Karr Server

        """
        t0 = time.time()
        sabiodb = sabio_rk.SabioRk(cache_dirname = self.cache_dirname, verbose = self.verbose)
        sabio_ses = sabiodb.session

        _entity = self.entity
        _property = self.property

        if self.max_entries == float('inf'):
            sabio_entry = sabio_ses.query(sabio_rk.Entry)
        else:
            sabio_entry = sabio_ses.query(sabio_rk.Entry).filter(sabio_rk.Entry._id.in_\
                (range(self.sabio_loaded+1, self.sabio_loaded+1 + int(self.max_entries*50))))

        counter = 1
        for item in sabio_entry:
            metadata = self.get_or_create_object(Metadata, name = 'Kinetic Law ' + str(item.id)) if item._type == 'kinetic_law' else self.get_or_create_object(Metadata, name = item.name)
            metadata.synonym = [self.get_or_create_object(Synonym, name = synonyms.name) for synonyms in item.synonyms]
            metadata.taxon = [self.get_or_create_object(Taxon, ncbi_id = docs.id) for docs in item.cross_references if docs.namespace == 'taxonomy']
            uniprot = [docs.id for docs in item.cross_references  if docs.namespace == 'uniprot']
            metadata.resource = [self.get_or_create_object(Resource, namespace = docs.namespace, _id = docs.id) for docs in item.cross_references]
            compartment = self.get_or_create_object(CellCompartment, name = item.name) if item._type == 'compartment' else None
            if compartment: continue

            if item._type == 'compound':
                structure = item.structures
                for struct in structure:
                    _property.structure = self.get_or_create_object(Structure, type = 'Structure', name = item.name,
                        _value_smiles = struct.value, _value_inchi = struct._value_inchi,
                        _structure_formula_connectivity = struct._value_inchi_formula_connectivity, _metadata = metadata)\
                        if struct.format == 'smiles' else None
                _entity.compound = self.get_or_create_object(Compound, type = 'Compound',
                    name = item.name, compound_name = item.name,
                    _is_name_ambiguous = sabio_ses.query(sabio_rk.Compound).get(item._id)._is_name_ambiguous,
                    structure = _property.structure, _metadata = metadata)
                continue

            elif item._type == 'enzyme':
                complx = sabio_ses.query(sabio_rk.Enzyme).get(item._id)
                _entity.protein_complex = self.get_or_create_object(ProteinComplex, type = 'Enzyme' ,
                    name = item.name , complex_name = item.name,
                    molecular_weight = item.molecular_weight, funcat_dsc = 'Enzyme', _metadata = metadata)
                continue

            elif item._type == 'enzyme_subunit':
                result = self.session.query(ProteinComplex).filter_by(complex_name = item.enzyme.name).first()
                _entity.protein_subunit = self.get_or_create_object(ProteinSubunit, type = 'Enzyme Subunit',
                    name = item.name, subunit_name = item.name,
                    uniprot_id = uniprot, coefficient = item.coefficient,
                    molecular_weight = item.molecular_weight,
                    proteincomplex = result, _metadata = metadata)
                continue

            elif item._type == 'kinetic_law':

                catalyst = self.session.query(ProteinComplex).filter_by(complex_name = item.enzyme.name).first() if item.enzyme_id else None
                metadata.resource = [self.get_or_create_object(Resource, namespace = resource.namespace, _id = resource.id) for resource in item.references]
                metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = item.taxon))
                metadata.cell_line.append(self.get_or_create_object(CellLine, name = item.taxon_variant))
                metadata.conditions.append(self.get_or_create_object(Conditions, temperature = item.temperature, ph = item.ph, media = item.media))
                _property.kinetic_law = self.get_or_create_object(KineticLaw, type = 'Kinetic Law', enzyme = catalyst,
                    enzyme_type = item.enzyme_type, tissue = item.tissue, mechanism = item.mechanism, equation = item.equation, _metadata = metadata)

                def common_schema_compound(sabio_object):
                    compound_name = sabio_object.name
                    return self.session.query(Compound).filter_by(compound_name = compound_name).first()

                def common_schema_compartment(sabio_object):
                    if sabio_object:
                        compartment_name = sabio_object.name
                        return self.session.query(CellCompartment).filter_by(name = compartment_name).first()
                    else: return None

                reactants = [Reaction(compound = common_schema_compound(r.compound),
                    compartment = common_schema_compartment(r.compartment), _is_reactant = 1, rxn_type = r.type,
                     kinetic_law_id = _property.kinetic_law.kineticlaw_id) for r in item.reactants if item.reactants]

                products = [Reaction(compound = common_schema_compound(p.compound),
                    compartment = common_schema_compartment(p.compartment), _is_product = 1, rxn_type = p.type,
                    kinetic_law_id = _property.kinetic_law.kineticlaw_id) for p in item.products if item.products]

                modifier = [Reaction(compound = common_schema_compound(m.compound),
                    compartment = common_schema_compartment(m.compartment), _is_modifier = 1, rxn_type = m.type,
                    kinetic_law_id = _property.kinetic_law.kineticlaw_id) for m in item.modifiers if item.products]

                for param in item.parameters:
                    parameter = Parameter(sabio_type = param.type, value = param.value, error = param.error,
                        units = param.units, observed_name = param.observed_name, kinetic_law = _property.kinetic_law,
                        observed_sabio_type = param.observed_type, observed_value = param.observed_value,
                        observed_error = param.observed_error, observed_units = param.observed_units)
                    parameter.compound = common_schema_compound(param.compound) if param.compound else None
                continue


        self.get_or_create_object(Progress, database_name = 'Sabio', amount_loaded = self.sabio_loaded+ (self.max_entries*50))

        if self.verbose:
            print('Total time taken for Sabio: ' + str(time.time()-t0) + ' secs')

        if self.verbose:
            print('Comitting')
        self.session.commit()

    def add_intact(self):
        t0 = time.time()
        intactdb = intact.IntAct(cache_dirname = self.cache_dirname)

        intact_ses = intactdb.session
        _entity = self.entity
        _property = self.property

        if self.max_entries == float('inf'):
            interactiondb = intactdb.session.query(intact.ProteinInteractions).all()
        else:
            interactiondb = intactdb.session.query(intact.ProteinInteractions).filter(intact.ProteinInteractions.index.in_\
                (range(self.max_entries))).all()

        self.session.bulk_insert_mappings(ProteinInteractions,
            [{'name' : e.interactor_a + "+" + e.interactor_b, 'type' : 'Protein Interaction',
            'participant_a' : e.interactor_a, 'participant_b' : e.interactor_b, 'publication' : e.publications,
            'interaction' : e.interaction, 'site_a' : e.feature_a, 'site_b' : e.feature_b,
            'stoich_a' : e.stoich_a, 'stoich_b' : e.stoich_b, 'interaction_type': e.interaction_type} for e in interactiondb]
            )

        index = self.intact_loaded
        for row in self.session.query(ProteinInteractions).all():
            metadata = self.get_or_create_object(Metadata, name = 'Protein Interaction ' + str(index))
            metadata.resource.append(self.get_or_create_object(Resource, namespace = 'pubmed', _id = re.split(':||', row.publication)[1]))
            row._metadata = metadata
            if 'uniprotkb:' in  row.participant_a:
                row.protein_subunit.append(self.get_or_create_object(ProteinSubunit,\
                uniprot_id = str(row.participant_a.replace('uniprotkb:', ''))))
            if 'uniprotkb:' in  row.participant_b:
                row.protein_subunit.append(self.get_or_create_object(ProteinSubunit,\
                uniprot_id = str(row.participant_b.replace('uniprotkb:', ''))))
            index += 1

        self.get_or_create_object(Progress, database_name = 'IntAct', amount_loaded = self.intact_loaded + index)

        if self.verbose:
            print('Total time taken for IntAct: ' + str(time.time()-t0) + ' secs')

        if self.verbose:
            print('Comitting')
        self.session.commit()

    # def add_uniprot(self):
    #     t0 = time.time()
    #
    #     unidb = uniprot.Uniprot(cache_dirname = self.cache_dirname)
    #     unidb_ses = unidb.session
    #     _entity = self.entity
    #     _property = self.property
    #
    #     com_unis = self.session.query(ProteinSubunit).all()
    #
    #     for subunit in com_unis:
    #         info = unidb_ses.query(uniprot.UniprotData).filter_by(uniprot_id = subunit.uniprot_id).first()
    #         if info:
    #             subunit.subunit_name = info.entry_name if !(subunit.subunit_name)
    #             subunit.entrez_id = info.entrez_id if not subunit.entrez_id else None
    #             # subunit.gene_name = info.gene_name
    #             subunit.canonical_sequence = info.canonical_sequence if not subunit.canonical_sequence else None
    #             subunit.length = info.length if not subunit.length else None
    #             subunit.mass = info.mass if not subunit.mass else None
    #
    #     if self.verbose:
    #         print('Total time taken for Uniprot: ' + str(time.time()-t0) + ' secs')
    #
    #     if self.verbose:
    #         print('Comitting')
    #     self.session.commit()
