# -*- coding: utf-8 -*-

"""
This code is a common schema for all the kinetic_datanator modules

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-07-31
:Copyright: 2017, Karr Lab
:License: MIT
"""

from sqlalchemy import Column, BigInteger, Integer, Float, String, Text, ForeignKey, Boolean, Table, create_engine, Numeric, or_
from sqlalchemy.orm import relationship, backref, sessionmaker
from kinetic_datanator.core import data_source
from kinetic_datanator.data_source import corum, pax, jaspar, array_express, ecmdb, sabio_rk
import sqlalchemy.ext.declarative
from six import BytesIO
import six
from bioservices import UniProt
from ete3 import NCBITaxa
import pandas as pd
import os
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

    Attributes:
        name (:obj:`str`): Name of the Cell Compartment

    """
    __tablename__ = 'cell_compartment'

    id = Column(Integer, primary_key = True)
    name = Column(String(255), unique = True, index = True)

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
    compound_name = Column(String(255), unique = True)
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

    enzyme_id = Column(Integer, ForeignKey('protein_complex.complex_id'), index = True)
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

    complex_id = Column(Integer, ForeignKey('protein_complex.complex_id') )
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

    subunit_id = Column(Integer, ForeignKey('protein_subunit.subunit_id'), index = True)
    subunit = relationship('ProteinSubunit', backref = 'pax_abundance_data' )

class DNABindingData(Base):
    """
    Represents Matrix binding profile for a protein transcription factor

    Attributes:
        position (:obj:`int`): Position in the sequence
        frequency_a (:obj:`int`): Frequency of A
        frequency_c (:obj:`int`): Frequency of C
        frequency_g (:obj:`int`): Frequency of G
        frequency_t (:obj:`int`): Frequency of T
        dataset_id  (:obj:`int`): Represents the dataset from which the data stems from
    """
    __tablename__ = 'dna_bidning_data'

    position_id = Column(Integer, primary_key = True)
    position = Column(Integer, index=True)
    frequency_a = Column(Integer)
    frequency_c = Column(Integer)
    frequency_g = Column(Integer)
    frequency_t = Column(Integer)

    dataset_id  = Column(Integer, ForeignKey('dna_binding_dataset.dataset_id'))
    dataset = relationship('DNABindingDataset', backref = 'dna_binding_data', foreign_keys=[dataset_id])


class CommonSchema(data_source.HttpDataSource):
    """
    A Local SQLlite copy of the aggregation of data_source modules

    """
    base_model = Base

    def load_content(self):
        ## Initiate Observation and direct Subclasses
        observation = Observation()
        observation.physical_entity = PhysicalEntity()
        self.entity = observation.physical_entity
        observation.physical_property = PhysicalProperty()
        self.property = observation.physical_property

        ## Switches
        self.clear = False
        self.load = False
        self.download = True

        ## Add all DBs
        self.add_paxdb()
        if self.verbose:
            print('Pax Done')
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
        self.fill_missing_subunit_info()

        ## Add missing Taxon information
        self.fill_missing_ncbi_names()

        if self.verbose:
            print('Comitting')
        self.session.commit()

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

    def fill_missing_subunit_info(self):
        u = UniProt(verbose = False)

        subunits = self.session.query(ProteinSubunit).all()
        uni = []
        for items in subunits:
            uni.append(str(items.uniprot_id))
        entrez_dict = u.mapping(fr = 'ACC', to = 'P_ENTREZGENEID', query = uni)
        for protein in subunits:
            # IF statement created to account for issues in UniProt Entrez ID fetching
            if protein.entrez_id == None and protein.uniprot_id in entrez_dict.keys():
                protein.entrez_id = int(entrez_dict[protein.uniprot_id][0])
            subunit_data = pd.read_csv(six.StringIO(u.search(protein.uniprot_id, \
                columns = 'id, entry name, protein names, genes, sequence, length, mass', limit = 1)), \
                sep = '\t')
            protein.subunit_name = str(subunit_data['Protein names'].iloc[0])
            protein.gene_name = str(subunit_data['Gene names'].iloc[0])
            protein.canonical_sequence = str(subunit_data['Sequence'].iloc[0])
            protein.mass = int(subunit_data['Mass'].iloc[0].replace(',',''))
            protein.length = int(subunit_data['Length'].iloc[0])

    def fill_missing_ncbi_names(self):
        ncbi = NCBITaxa()

        ncbi_ids = []
        species = self.session.query(Taxon).all()
        for items in species:
            ncbi_ids.append(items.ncbi_id)

        species_dict = ncbi.get_taxid_translator(ncbi_ids)

        for tax in species:
            tax.name = species_dict[tax.ncbi_id]

    def add_paxdb(self):
        paxdb = pax.Pax(cache_dirname = self.cache_dirname, clear_content = self.clear,
            load_content= self.load, download_backup= self.download, max_entries = self.max_entries/5, verbose = self.verbose)
        pax_ses = paxdb.session

        _entity = self.entity
        _property = self.property

        pax_dataset = pax_ses.query(pax.Dataset).all()
        entries = 0
        for item in pax_dataset:
            if entries < (self.max_entries/5):
                metadata = self.get_or_create_object(Metadata, name = item.file_name)
                metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = item.taxon_ncbi_id))
                metadata.resource.append(self.get_or_create_object(Resource, namespace = 'url', _id = item.publication))
                _property.abundance_dataset = self.get_or_create_object(AbundanceDataSet, type = 'Protein Abundance Dataset',
                    name = item.file_name, file_name = item.file_name, score = item.score, weight = item.weight, coverage= item.coverage, _metadata = metadata)
                abundance = pax_ses.query(pax.Observation).filter_by(dataset_id = item.id).all()
                for data in abundance:
                    uniprot_id = pax_ses.query(pax.Protein).get(data.protein_id).uniprot_id
                    _entity.protein_subunit = self.get_or_create_object(ProteinSubunit,
                        type = 'Protein Subunit', uniprot_id = uniprot_id, _metadata = metadata)
                    abundance_data = self.get_or_create_object(AbundanceData,
                        abundance = data.abundance, dataset = _property.abundance_dataset,
                        subunit = _entity.protein_subunit)
                entries += 1

    def add_corumdb(self):
        corumdb = corum.Corum(cache_dirname = self.cache_dirname, clear_content = self.clear,
            load_content= self.load, download_backup= self.download, max_entries = self.max_entries, verbose = self.verbose)
        corum_ses = corumdb.session

        _entity = self.entity
        _property = self.property

        corum_complex = corum_ses.query(corum.Complex).all()
        corum_subunit = corum_ses.query(corum.Subunit).all()

        entries = 0
        for row in corum_complex:
            if entries < self.max_entries:
                entry = corum_ses.query(corum.Observation).get(row.observation_id)
                metadata = self.get_or_create_object(Metadata, name = row.complex_name)
                metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = entry.taxon_ncbi_id))
                metadata.resource.append(self.get_or_create_object(Resource, namespace = 'pubmed', _id = str(entry.pubmed_id)))
                metadata.method.append(self.get_or_create_object(Method, name = 'purification', comments = entry.pur_method))
                metadata.cell_line.append(self.get_or_create_object(CellLine, name = entry.cell_line))
                _entity.protein_complex = self.get_or_create_object(ProteinComplex,
                    type = 'Protein Complex', name = row.complex_name, complex_name = row.complex_name, go_id = row.go_id,
                    go_dsc = row.go_dsc, funcat_id = row.funcat_id, funcat_dsc = row.funcat_dsc, su_cmt = row.su_cmt,
                    complex_cmt = row.complex_cmt, disease_cmt = row.disease_cmt,  _metadata = metadata)
                entries += 1

        entries = 0
        for row in corum_subunit:
            if entries < self.max_entries:
                corum_cmplx = corum_ses.query(corum.Complex).get(row.complex_id)
                complex_ = self.session.query(ProteinComplex).filter_by(complex_name = corum_cmplx.complex_name).first()
                entry = self.session.query(Observation).get(complex_.complex_id)
                _entity.protein_subunit = self.get_or_create_object(ProteinSubunit,
                    type = 'Protein Subunit', uniprot_id = row.su_uniprot,
                    entrez_id = row.su_entrezs, name = row.protein_name, subunit_name = row.protein_name, gene_name=row.gene_name,
                    gene_syn = row.gene_syn, proteincomplex = complex_, _metadata = self.session.query(Metadata).get(entry._metadata_id))
                entries += 1

    def add_jaspardb(self):
        jaspardb = jaspar.Jaspar(cache_dirname = self.cache_dirname, clear_content = self.clear,
            load_content= self.load, download_backup= self.download, verbose = self.verbose)
        jasp_ses = jaspardb.session

        _entity = self.entity
        _property = self.property

        jaspar_matrix = jasp_ses.query(jaspar.Matrix).all()

        entries = 0
        for item in jaspar_matrix:
            if entries < self.max_entries:
                temp_tf = jasp_ses.query(jaspar.TranscriptionFactor).filter_by(id = item.transcription_factor_id).first()
                tax = jasp_ses.query(jaspar.transcription_factor_species).filter_by(transcription_factor__id = temp_tf._id).all()
                ref = jasp_ses.query(jaspar.matrix_resource).filter_by(matrix_id = item.id).all()
                type_name = jasp_ses.query(jaspar.Type).get(item.type_id).name
                temp_class = jasp_ses.query(jaspar.transcription_factor_class).filter_by(transcription_factor__id = temp_tf._id).first()
                if temp_class is not None:
                    class_name = jasp_ses.query(jaspar.Class).get(temp_class.class_id).name
                temp_fam = jasp_ses.query(jaspar.transcription_factor_family).filter_by(transcription_factor__id = temp_tf._id).first()
                if temp_fam is not None:
                    fam_name = jasp_ses.query(jaspar.Family).get(temp_fam.family_id).name
                for speice in tax:
                    ncbi = jasp_ses.query(jaspar.Species).get(speice.species_id).ncbi_id
                    for docs in ref:
                        doc = docs.resource_pubmed_id
                        if '::' in temp_tf.name:
                            su = temp_tf.name.split('::')
                            dialoge = 'Subunits are: '
                            for items in su:
                                dialoge += (items + ' ')
                            metadata = self.get_or_create_object(Metadata, name = temp_tf.name)
                            metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = ncbi))
                            metadata.method.append(self.get_or_create_object(Method, name = type_name))
                            metadata.resource.append(self.get_or_create_object(Resource, namespace = 'pubmed', _id = doc))
                            _entity.protein_complex = self.get_or_create_object(ProteinComplex, type = 'Transcription Factor Complex', name = temp_tf.name,
                                complex_name = temp_tf.name, su_cmt = dialoge, complex_cmt = 'transcription factor', class_name = class_name, family_name = fam_name, _metadata = metadata)
                            _property.dna_binding_dataset = self.get_or_create_object(DNABindingDataset, type = 'DNA Binding Dataset', name = temp_tf.name,
                                version = item.version, tf = _entity.protein_complex, _metadata = metadata)
                            for pos in jaspar_matrixposition:
                                data = self.get_or_create_object(DNABindingData, position = pos.position,
                                frequency_a = pos.frequency_a, frequency_c = pos.frequency_c, frequency_g = pos.frequency_g,
                                frequency_t = pos.frequency_t, dataset = _property.dna_binding_dataset)
                        else:
                            sub_id = jasp_ses.query(jaspar.transcription_factor_subunit).filter_by(transcription_factor__id = temp_tf._id).first().subunit_id
                            uniprot_id = jasp_ses.query(jaspar.Subunit).get(sub_id).uniprot_id
                            metadata = self.get_or_create_object(Metadata, name = temp_tf.name)
                            metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = ncbi))
                            metadata.method.append(self.get_or_create_object(Method, name = type_name))
                            metadata.resource.append(self.get_or_create_object(Resource, namespace = 'pubmed', _id = doc))
                            _entity.protein_subunit = self.get_or_create_object(ProteinSubunit, uniprot_id = uniprot_id, type = 'Transcription Factor Subunit',
                                name = temp_tf.name, subunit_name = temp_tf.name, gene_name = temp_tf.name,
                                class_name = class_name, family_name = fam_name, _metadata = metadata)
                            _property.dna_binding_dataset = self.get_or_create_object(DNABindingDataset, type = 'DNA Binding Dataset', name = temp_tf.name,
                                version = item.version, subunit = _entity.protein_subunit, _metadata = metadata)
                            jaspar_matrixposition = jasp_ses.query(jaspar.MatrixPosition).filter_by(matrix_id = item.id).all()
                            for pos in jaspar_matrixposition:
                                data = self.get_or_create_object(DNABindingData, position = pos.position,
                                frequency_a = pos.frequency_a, frequency_c = pos.frequency_c, frequency_g = pos.frequency_g,
                                frequency_t = pos.frequency_t, dataset = _property.dna_binding_dataset)
                entries += 1

    def add_ecmdb(self):
        ecmDB = ecmdb.Ecmdb(cache_dirname = self.cache_dirname, clear_content = self.clear,
            load_content= self.load, download_backup= self.download, max_entries = self.max_entries, verbose = self.verbose)
        ecm_ses = ecmDB.session

        _entity = self.entity
        _property = self.property

        ecmdb_compound = ecm_ses.query(ecmdb.Compound).all()

        entries = 0
        for item in ecmdb_compound:
            if entries < self.max_entries:
                concentration = ecm_ses.query(ecmdb.Concentration).filter_by(compound_id = item._id).all()
                ref = ecm_ses.query(ecmdb.compound_resource).filter_by(compound__id = item._id).all()
                compart = ecm_ses.query(ecmdb.compound_compartment).filter_by(compound__id = item._id).all()
                syn = ecm_ses.query(ecmdb.compound_synonym).filter_by(compound__id = item._id).all()
                metadata = self.get_or_create_object(Metadata, name = item.name)
                metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = 562, name = 'E.Coli'))
                for docs in ref:
                    resource = ecm_ses.query(ecmdb.Resource).get(docs.resource__id)
                    metadata.resource.append(self.get_or_create_object(Resource, namespace = resource.namespace, _id = resource.id))
                for areas in compart:
                    compartment = ecm_ses.query(ecmdb.Compartment).get(areas.compartment__id)
                    metadata.cell_compartment.append(self.get_or_create_object(CellCompartment, name = compartment.name))
                for syns in syn:
                    synonym = ecm_ses.query(ecmdb.Synonym).get(syns.synonym__id)
                    metadata.synonym.append(self.get_or_create_object(Synonym, name = synonym.name))
                _property.structure = self.get_or_create_object(Structure, type = 'Structure', name = item.name,
                    _value_inchi = item.structure,
                   _structure_formula_connectivity = item._structure_formula_connectivity, _metadata = metadata)
                _entity.compound = self.get_or_create_object(Compound, type = 'Compound', name = item.name,
                    compound_name = item.name, description = item.description, comment = item.comment, structure = _property.structure,
                    _metadata = metadata)
                if concentration:
                    index = 0
                    for rows in concentration:
                        metadata.cell_line.append(self.get_or_create_object(CellLine, name = rows.strain))
                        metadata.conditions.append(self.get_or_create_object(Conditions, growth_status = rows.growth_status,
                            media = rows.media, temperature = rows.temperature, growth_system = rows.growth_system))
                        refs = ecm_ses.query(ecmdb.concentration_resource).filter_by(concentration__id = rows._id).all()
                        for docs in refs:
                            resource = ecm_ses.query(ecmdb.Resource).get(docs.resource__id)
                            metadata.resource.append(self.get_or_create_object(Resource, namespace = resource.namespace, _id = resource.id))
                        _property.concentration = self.get_or_create_object(Concentration, type = 'Concentration', name = item.name+ ' Concentration '+str(index),
                            value = rows.value, error = rows.error, _metadata = metadata, compound = _entity.compound)
                        index += 1
                entries += 1

    def add_sabiodb(self):
        sabiodb = sabio_rk.SabioRk(cache_dirname = self.cache_dirname, clear_content = self.clear,
            load_content= self.load, download_backup= self.download, max_entries = self.max_entries*5, verbose = self.verbose)
        sabio_ses = sabiodb.session

        _entity = self.entity
        _property = self.property

        sabio_entry = sabio_ses.query(sabio_rk.Entry).all()
        entries = 0
        counter = 1
        for item in sabio_entry:
            if entries < (self.max_entries*5):
                if item.name:
                    metadata = self.get_or_create_object(Metadata, name = item.name)
                elif item._type == 'kinetic_law':
                    metadata = self.get_or_create_object(Metadata, name = 'Kinetic Law ' + str(counter))
                    counter += 1
                syn = sabio_ses.query(sabio_rk.entry_synonym).filter_by(entry__id = item._id).all()
                res = sabio_ses.query(sabio_rk.entry_resource).filter_by(entry__id = item._id).all()
                for synonyms in syn:
                    syns = sabio_ses.query(sabio_rk.Synonym).get(synonyms.synonym__id)
                    metadata.synonym.append(self.get_or_create_object(Synonym, name = syns.name))
                for docs in res:
                    resource = sabio_ses.query(sabio_rk.Resource).get(docs.resource__id)
                    if resource.namespace == 'taxonomy':
                        metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = resource.id))
                    elif resource.namespace == 'uniprot':
                        uniprot = resource.id
                    else:
                        metadata.resource.append(self.get_or_create_object(Resource, namespace = resource.namespace,
                            _id = resource.id))
                if item._type == 'compound':
                    structure = sabio_ses.query(sabio_rk.compound_compound_structure).filter_by(compound__id = item._id).all()
                    if structure:
                        for shape in structure:
                            struct = sabio_ses.query(sabio_rk.CompoundStructure).get(shape.compound_structure__id)
                            if struct.format == 'smiles':
                                _property.structure = self.get_or_create_object(Structure, type = 'Structure', name = item.name,
                                _value_smiles = struct.value, _value_inchi = struct._value_inchi,
                                _structure_formula_connectivity = struct._value_inchi_formula_connectivity, _metadata = metadata)
                                break
                    else:
                        _property.structure = None
                    _entity.compound = self.get_or_create_object(Compound, type = 'Compound',
                        name = item.name, compound_name = item.name,
                        _is_name_ambiguous = sabio_ses.query(sabio_rk.Compound).get(item._id)._is_name_ambiguous,
                        structure = _property.structure, _metadata = metadata)
                elif item._type == 'enzyme':
                    complx = sabio_ses.query(sabio_rk.Enzyme).get(item._id)
                    _entity.protein_complex = self.get_or_create_object(ProteinComplex, type = 'Enzyme' ,
                        name = item.name , complex_name = item.name,
                        molecular_weight = complx.molecular_weight, funcat_dsc = 'Enzyme', _metadata = metadata)
                elif item._type == 'enzyme_subunit':
                    subunit = sabio_ses.query(sabio_rk.EnzymeSubunit).get(item._id)
                    complx = sabio_ses.query(sabio_rk.Entry).get(subunit.enzyme_id)
                    result = self.session.query(ProteinComplex).filter_by(complex_name = complx.name).first()
                    _entity.protein_subunit = self.get_or_create_object(ProteinSubunit, type = 'Enzyme Subunit',
                        name = item.name, subunit_name = item.name,
                        uniprot_id = uniprot, coefficient = subunit.coefficient,
                        molecular_weight = subunit.molecular_weight,
                        proteincomplex = result, _metadata = metadata)
                elif item._type == 'kinetic_law':
                    res = sabio_ses.query(sabio_rk.kinetic_law_resource).filter_by(kinetic_law__id = item._id).all()
                    law = sabio_ses.query(sabio_rk.KineticLaw).get(item._id)
                    entry = sabio_ses.query(sabio_rk.Entry).get(law.enzyme_id)
                    result = self.session.query(ProteinComplex).filter_by(complex_name = entry.name).first()
                    for docs in res:
                        resource = sabio_ses.query(sabio_rk.Resource).get(docs.resource__id)
                        metadata.resource.append(self.get_or_create_object(Resource, namespace = resource.namespace, _id = resource.id))
                    metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = law.taxon))
                    metadata.cell_line.append(self.get_or_create_object(CellLine, name = law.taxon_variant))
                    metadata.conditions.append(self.get_or_create_object(Conditions, temperature = law.temperature, ph = law.ph, media = law.media))
                    _property.kinetic_law = self.get_or_create_object(KineticLaw, type = 'Kinetic Law', enzyme = result,
                        enzyme_type = law.enzyme_type, tissue = law.tissue, mechanism = law.mechanism, equation = law.equation, _metadata = metadata)
                    rxn = sabio_ses.query(sabio_rk.ReactionParticipant).filter(or_(sabio_rk.ReactionParticipant.reactant_kinetic_law_id == law._id, sabio_rk.ReactionParticipant.product_kinetic_law_id == law._id, sabio_rk.ReactionParticipant.modifier_kinetic_law_id == law._id)).all()
                    for row in rxn:
                        compound_name = sabio_ses.query(sabio_rk.Entry).get(row.compound_id).name
                        _compound = self.session.query(Compound).filter_by(compound_name = compound_name).first()
                        if row.reactant_kinetic_law_id:
                            reaction = Reaction(compound = _compound,
                                coefficient = row.coefficient, _is_reactant = 1, rxn_type = row.type,
                                kinetic_law_id = _property.kinetic_law.kineticlaw_id)
                        elif row.product_kinetic_law_id:
                            reaction = Reaction( compound = _compound,
                                coefficient = row.coefficient, _is_product = 1, rxn_type = row.type,
                                kinetic_law_id = _property.kinetic_law.kineticlaw_id)
                        elif row.modifier_kinetic_law_id:
                            reaction = Reaction(compound = _compound,
                                coefficient = row.coefficient, _is_modifier = 1, rxn_type = row.type,
                                kinetic_law_id = _property.kinetic_law.kineticlaw_id)
                    total_param = sabio_ses.query(sabio_rk.Parameter).filter_by(kinetic_law_id = law._id).all()
                    for param in total_param:
                        if param.compound_id:
                            compound_name = sabio_ses.query(sabio_rk.Entry).get(param.compound_id).name
                            _compound = self.session.query(Compound).filter_by(compound_name = compound_name).first()
                            parameter = Parameter(sabio_type = param.type, value = param.value, error = param.error,
                                units = param.units, observed_name = param.observed_name, kinetic_law = _property.kinetic_law,
                                observed_sabio_type = param.observed_type, observed_value = param.observed_value, compound = _compound,
                                observed_error = param.observed_error, observed_units = param.observed_units)
                        else:
                            parameter = Parameter(sabio_type = param.type, value = param.value, error = param.error,
                                units = param.units, observed_name = param.observed_name, kinetic_law = _property.kinetic_law,
                                observed_sabio_type = param.observed_type, observed_value = param.observed_value,
                                observed_error = param.observed_error, observed_units = param.observed_units)
                entries += 1


    # def add_arrayexpressdb(self):
    #     arrayexpressdb = array_express.ArrayExpress(name = 'array_express', load_content=False, download_backup=False)
    #     arrayexpressdb.load_content
    #     array_session = arrayexpressdb.session
