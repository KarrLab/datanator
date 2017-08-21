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

Base = sqlalchemy.ext.declarative.declarative_base()

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
    sequence = Column(String(255))
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

    compound_id (:obj:`int`): Common Schema Observation Identifier
    compound_name (:obj:`str`): Name of the Compound
    description
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

    """
    observation_id = Column(Integer, ForeignKey('observation.id'), primary_key = True)
    type = Column(String(255))
    name = Column(String(255))

    __tablename__ = 'physical_property'
    __mapper_args__ = {'polymorphic_identity': 'physical_property'}

class Structure(PhysicalProperty):
    """


    """

    __tablename__ = 'structure'
    __mapper_args__ = {'polymorphic_identity': 'structure'}

    struct_id = Column(Integer, ForeignKey('physical_property.observation_id'), primary_key = True)
    _value_smiles= Column(String(255))
    _value_inchi = Column(String(255))
    _structure_formula_connectivity = Column(String(255))

# class Reaction(PhysicalProperty):
#     __tablename__ = 'reaction'
#     __mapper_args__ = {'polymorphic_identity': 'reaction'}

class Concentration(PhysicalProperty):
    """

    """

    __tablename__ = 'concentration'
    __mapper_args__ = {'polymorphic_identity': 'concentration'}

    concentration_id = Column(Integer, ForeignKey('physical_property.observation_id'), primary_key = True)
    value = Column(Float)
    error = Column(Float)

class KineticLaw(PhysicalProperty):
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
    __tablename__ = 'reaction'
    # __mapper_args__ = {'polymorphic_identity': 'reaction'}

    reaction_id = Column(Integer, primary_key = True)
    compound_id = Column(Integer, ForeignKey('compound.compound_id'))
    compound = relationship(Compound, backref = 'reaction' )
    coefficient = Column(Float)
    _is_reactant = Column(Boolean)
    _is_product = Column(Boolean)
    _is_modifier = Column(Boolean)
    rxn_type = Column(String(255))

    kinetic_law_id = Column(Integer, ForeignKey('kinetic_law.kineticlaw_id'))

class Parameter(Base):
    """

    """

    __tablename__ = 'parameter'
    # __mapper_args__ = {'polymorphic_identity': 'parameter'}

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



class AbundanceDataSet(PhysicalProperty):
    """

    """

    __tablename__ = 'abundance_dataset'
    __mapper_args__ = {'polymorphic_identity': 'abundance_dataset'}

    dataset_id = Column(Integer, ForeignKey('physical_property.observation_id'), primary_key = True)
    file_name = Column(String, unique = True)
    score = Column(Float)
    weight = Column(Integer)
    coverage = Column(Integer)

class AbundanceData(Base):
    """

    """
    __tablename__ = 'abundance_data'
    # __mapper_args__ = {'polymorphic_identity': 'abundance_data'}

    abundance_id  = Column(Integer, primary_key = True)
    abundance = Column(Float)

    dataset_id  = Column(Integer, ForeignKey('abundance_dataset.dataset_id'))
    dataset = relationship('AbundanceDataSet', backref = 'abundance_data', foreign_keys=[dataset_id])

    subunit_id = Column(Integer, ForeignKey('protein_subunit.subunit_id'), index = True)
    subunit = relationship('ProteinSubunit', backref = 'pax_abundance_data' )

class DNABindingDataset(PhysicalProperty):
    """


    """
    __tablename__ = 'dna_binding_dataset'
    __mapper_args__ = {'polymorphic_identity': 'dna_binding_dataset'}

    dataset_id = Column(Integer, ForeignKey('physical_property.observation_id'), primary_key = True)
    version = Column(Integer)

    complex_id = Column(Integer, ForeignKey('protein_complex.complex_id') )
    tf = relationship('ProteinComplex', backref = 'dna_binding_dataset')

    subunit_id = Column(Integer, ForeignKey('protein_subunit.subunit_id'))
    subunit = relationship('ProteinSubunit', backref = 'dna_binding_dataset')

class DNABindingData(Base):
    """

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



class Metadata(Base):
    """

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

    """

    id = Column(Integer, primary_key = True)
    name = Column(String(255))
    comments = Column(String(255))

    __tablename__ = 'method'

class Taxon(Base):
    """

    """
    ncbi_id = Column(Integer, primary_key = True)
    name = Column(String(255))

    __tablename__ = 'taxon'

class Synonym(Base):
    """

    """

    id = Column(Integer, primary_key = True)
    name = Column(String(255))

    __tablename__ = 'synonym'

class Resource(Base):
    """


    """
    id = Column(Integer, primary_key = True)
    namespace = Column(String(255))
    _id =  Column(String(255))

    __tablename__ = 'resource'

class CellLine(Base):
    """

    """
    id = Column(Integer, primary_key = True)
    name = Column(String(255))

    __tablename__ = 'cell_line'

class Conditions(Base):
    """

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

    """
    __tablename__ = 'cell_compartment'

    id = Column(Integer, primary_key = True)
    name = Column(String(255), unique = True, index = True)


class CommonSchema(data_source.CachedDataSource):
    """
    A Local SQLlite copy of the aggregation of data_source modules

    NOTES: For now.. Keeping the SQL local, out of cache
    """
    base_model = Base

    def load_content(self):
        observation = Observation()
        observation.physical_entity = PhysicalEntity()
        self.obs_pe = observation.physical_entity
        observation.physical_property = PhysicalProperty()
        self.obs_pp = observation.physical_property

        ## Add all DBs
        self.add_all()

        if self.verbose:
            print('Comitting')
        self.session.commit()

    def add_all(self):
        self.add_paxdb()
        self.add_corumdb()
        self.add_jaspardb()
        self.add_ecmdb()
        self.add_sabiodb()

    def add_paxdb(self):
        if self._is_local:
            paxdb = pax.Pax(name = 'pax', clear_content = True,  load_content=False, download_backup=False, max_entries = self.max_entries)
        else:
            paxdb = pax.Pax(cache_dirname = self.cache_dirname, clear_content = True,  load_content=False, download_backup=False, max_entries = self.max_entries)
        paxdb.load_content()
        self.pax_session = paxdb.session

        obs_pe = self.obs_pe
        obs_pp = self.obs_pp

        pax_dataset = self.pax_session.query(pax.Dataset).all()

        for item in pax_dataset:
            metadata = self.get_or_create_object(Metadata, name = item.file_name)
            metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = item.taxon_ncbi_id,
                name = self.pax_session.query(pax.Taxon).get(item.taxon_ncbi_id).species_name))
            metadata.resource.append(self.get_or_create_object(Resource, namespace = 'url', _id = item.publication))
            obs_pp.abundance_dataset = self.get_or_create_object(AbundanceDataSet, type = 'Protein Abundance Dataset',
                name = item.file_name, file_name = item.file_name, score = item.score, weight = item.weight, coverage= item.coverage, _metadata = metadata )
            abundance = self.pax_session.query(pax.Observation).filter_by(dataset_id = item.id).all()
            for data in abundance:
                uniprot_id = self.pax_session.query(pax.Protein).get(data.protein_id).uniprot_id
                obs_pe.protein_subunit = self.get_or_create_object(ProteinSubunit,
                    type = 'Protein Subunit', name = uniprot_id, uniprot_id = uniprot_id,
                    _metadata = metadata)
                abundance_data = (self.get_or_create_object(AbundanceData,
                    abundance = data.abundance, dataset = obs_pp.abundance_dataset,
                    subunit = obs_pe.protein_subunit))

    def add_corumdb(self):
        if self._is_local:
            corumdb = corum.Corum(name = 'corum', clear_content = True, load_content=False, download_backup=False, max_entries = self.max_entries)
        else:
            corumdb = corum.Corum(cache_dirname = self.cache_dirname , clear_content = True, load_content=False, download_backup=False, max_entries = self.max_entries)
        corumdb.load_content()
        self.corum_session = corumdb.session

        obs_pe = self.obs_pe
        obs_pp = self.obs_pp

        corum_complex = self.corum_session.query(corum.Complex).all()
        corum_subunit = self.corum_session.query(corum.Subunit).all()

        for row in corum_complex:
            entry = self.corum_session.query(corum.Observation).get(row.observation_id)
            metadata = self.get_or_create_object(Metadata, name = row.complex_name)
            metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = entry.taxon_ncbi_id,
                name = self.corum_session.query(corum.Taxon).get(entry.taxon_ncbi_id).swissprot_id))
            metadata.resource.append(self.get_or_create_object(Resource, namespace = 'pubmed', _id = str(entry.pubmed_id)))
            metadata.method.append(self.get_or_create_object(Method, name = 'purification', comments = entry.pur_method))
            metadata.cell_line.append(self.get_or_create_object(CellLine, name = entry.cell_line))
            obs_pe.protein_complex = self.get_or_create_object(ProteinComplex,
                type = 'Protein Complex', name = row.complex_name, complex_name = row.complex_name, go_id = row.go_id,
                go_dsc = row.go_dsc, funcat_id = row.funcat_id, funcat_dsc = row.funcat_dsc, su_cmt = row.su_cmt,
                complex_cmt = row.complex_cmt, disease_cmt = row.disease_cmt,  _metadata = metadata)

        for row in corum_subunit:
            corum_cmplx = self.corum_session.query(corum.Complex).get(row.complex_id)
            complex_ = self.session.query(ProteinComplex).filter_by(complex_name = corum_cmplx.complex_name).first()
            entry = self.session.query(Observation).get(complex_.complex_id)
            obs_pe.protein_subunit = self.get_or_create_object(ProteinSubunit,
                type = 'Protein Subunit', uniprot_id = row.su_uniprot,
                entrez_id = row.su_entrezs, name = row.protein_name, subunit_name = row.protein_name, gene_name=row.gene_name,
                gene_syn = row.gene_syn, proteincomplex = complex_, _metadata = self.session.query(Metadata).get(entry._metadata_id))

    def add_jaspardb(self):
        if self._is_local:
            jaspardb = jaspar.Jaspar(name = 'jaspar', clear_content = True, load_content=False, download_backup=False, max_entries = self.max_entries)
        else:
            jaspardb = jaspar.Jaspar(cache_dirname = self.cache_dirname, clear_content = True, load_content=False, download_backup=False, max_entries = self.max_entries)
        jaspardb.load_content()
        self.jaspar_session = jaspardb.session

        obs_pe = self.obs_pe
        obs_pp = self.obs_pp

        jaspar_matrix = self.jaspar_session.query(jaspar.Matrix).all()
        jaspar_matrixposition = self.jaspar_session.query(jaspar.MatrixPosition).all()

        for item in jaspar_matrix:
            temp_tf = self.jaspar_session.query(jaspar.TranscriptionFactor).filter_by(id = item.transcription_factor_id).first()
            tax = self.jaspar_session.query(jaspar.transcription_factor_species).filter_by(transcription_factor__id = temp_tf._id).all()
            ref = self.jaspar_session.query(jaspar.matrix_resource).filter_by(matrix_id = item.id).all()
            type_name = self.jaspar_session.query(jaspar.Type).get(item.type_id).name
            temp_class = self.jaspar_session.query(jaspar.transcription_factor_class).filter_by(transcription_factor__id = temp_tf._id).first()
            if temp_class is not None:
                class_name = self.jaspar_session.query(jaspar.Class).get(temp_class.class_id).name
            temp_fam = self.jaspar_session.query(jaspar.transcription_factor_family).filter_by(transcription_factor__id = temp_tf._id).first()
            if temp_fam is not None:
                fam_name = self.jaspar_session.query(jaspar.Family).get(temp_fam.family_id).name
            for speice in tax:
                ncbi = self.jaspar_session.query(jaspar.Species).get(speice.species_id).ncbi_id
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
                        obs_pe.protein_complex = self.get_or_create_object(ProteinComplex, type = 'Transcription Factor Complex', name = temp_tf.name,
                            complex_name = temp_tf.name, su_cmt = dialoge, complex_cmt = 'transcription factor', class_name = class_name, family_name = fam_name, _metadata = metadata)
                        obs_pp.dna_binding_dataset = self.get_or_create_object(DNABindingDataset, type = 'DNA Binding Dataset', name = temp_tf.name,
                            version = item.version, tf = obs_pe.protein_complex, _metadata = metadata)
                        for pos in jaspar_matrixposition:
                            data = self.get_or_create_object(DNABindingData, position = pos.position,
                            frequency_a = pos.frequency_a, frequency_c = pos.frequency_c, frequency_g = pos.frequency_g,
                            frequency_t = pos.frequency_t, dataset = obs_pp.dna_binding_dataset)
                    else:
                        sub_id = self.jaspar_session.query(jaspar.transcription_factor_subunit).filter_by(transcription_factor__id = temp_tf._id).first().subunit_id
                        uniprot_id = self.jaspar_session.query(jaspar.Subunit).get(sub_id).uniprot_id
                        metadata = self.get_or_create_object(Metadata, name = temp_tf.name)
                        metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = ncbi))
                        metadata.method.append(self.get_or_create_object(Method, name = type_name))
                        metadata.resource.append(self.get_or_create_object(Resource, namespace = 'pubmed', _id = doc))
                        obs_pe.protein_subunit = self.get_or_create_object(ProteinSubunit, uniprot_id = uniprot_id, type = 'Transcription Factor Subunit',
                            name = temp_tf.name, subunit_name = temp_tf.name, gene_name = temp_tf.name,
                            class_name = class_name, family_name = fam_name, _metadata = metadata)
                        obs_pp.dna_binding_dataset = self.get_or_create_object(DNABindingDataset, type = 'DNA Binding Dataset', name = temp_tf.name,
                            version = item.version, subunit = obs_pe.protein_subunit, _metadata = metadata)
                        for pos in jaspar_matrixposition:
                            data = self.get_or_create_object(DNABindingData, position = pos.position,
                            frequency_a = pos.frequency_a, frequency_c = pos.frequency_c, frequency_g = pos.frequency_g,
                            frequency_t = pos.frequency_t, dataset = obs_pp.dna_binding_dataset)

    def add_ecmdb(self):
        if self._is_local:
            ecmDB = ecmdb.Ecmdb(name = 'ecmdb', clear_content = True, load_content=False, download_backup=False, max_entries = self.max_entries)
        else:
            ecmDB = ecmdb.Ecmdb(cache_dirname = self.cache_dirname, clear_content = True, load_content=False, download_backup=False, max_entries = self.max_entries)
        ecmDB.load_content()
        self.ecmdb_session = ecmDB.session

        obs_pe = self.obs_pe
        obs_pp = self.obs_pp

        ##TODO: Figure out a connection between concentration and compound
        ecmdb_compound = self.ecmdb_session.query(ecmdb.Compound).all()

        for item in ecmdb_compound:
            concentration = self.ecmdb_session.query(ecmdb.Concentration).filter_by(compound_id = item._id).all()
            ref = self.ecmdb_session.query(ecmdb.compound_resource).filter_by(compound__id = item._id).all()
            compart = self.ecmdb_session.query(ecmdb.compound_compartment).filter_by(compound__id = item._id).all()
            syn = self.ecmdb_session.query(ecmdb.compound_synonym).filter_by(compound__id = item._id).all()
            metadata = self.get_or_create_object(Metadata, name = item.name)
            metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = 562, name = 'E.Coli'))
            for docs in ref:
                resource = self.ecmdb_session.query(ecmdb.Resource).get(docs.resource__id)
                metadata.resource.append(self.get_or_create_object(Resource, namespace = resource.namespace, _id = resource.id))
            for areas in compart:
                compartment = self.ecmdb_session.query(ecmdb.Compartment).get(areas.compartment__id)
                metadata.cell_compartment.append(self.get_or_create_object(CellCompartment, name = compartment.name))
            for syns in syn:
                synonym = self.ecmdb_session.query(ecmdb.Synonym).get(syns.synonym__id)
                metadata.synonym.append(self.get_or_create_object(Synonym, name = synonym.name))
            if concentration:
                for rows in concentration:
                    metadata.cell_line.append(self.get_or_create_object(CellLine, name = rows.strain))
                    metadata.conditions.append(self.get_or_create_object(Conditions, growth_status = rows.growth_status,
                        media = rows.media, temperature = rows.temperature, growth_system = rows.growth_system))
                    refs = self.ecmdb_session.query(ecmdb.concentration_resource).filter_by(concentration__id = rows._id).all()
                    for docs in refs:
                        resource = self.ecmdb_session.query(ecmdb.Resource).get(docs.resource__id)
                        metadata.resource.append(self.get_or_create_object(Resource, namespace = resource.namespace, _id = resource.id))
                    obs_pp.concentration = self.get_or_create_object(Concentration, type = 'Concentration', name = item.name,
                        value = rows.value, error = rows.error, _metadata = metadata)
            obs_pp.structure = self.get_or_create_object(Structure, type = 'Structure', name = item.name,
                _value_inchi = item.structure,
               _structure_formula_connectivity = item._structure_formula_connectivity, _metadata = metadata)
            obs_pe.compound = self.get_or_create_object(Compound, type = 'Compound', name = item.name,
                compound_name = item.name, description = item.description, comment = item.comment, structure = obs_pp.structure,
                _metadata = metadata)

    def add_sabiodb(self):
        if self._is_local:
            sabiodb = sabio_rk.SabioRk(name = 'sabio', clear_content = True, load_content=False, download_backup=False, max_entries = self.max_entries)
        else:
            sabiodb = sabio_rk.SabioRk(cache_dirname = self.cache_dirname, clear_content = True, load_content=False, download_backup=False, max_entries = self.max_entries)
        sabiodb.load_content()
        self.sabio_session = sabiodb.session

        obs_pe = self.obs_pe
        obs_pp = self.obs_pp

        sabio_entry = self.sabio_session.query(sabio_rk.Entry).all()
        counter = 1
        for item in sabio_entry:
            if item.name:
                metadata = self.get_or_create_object(Metadata, name = item.name)
            elif item._type == 'kinetic_law':
                metadata = self.get_or_create_object(Metadata, name = 'Kinetic Law ' + str(counter))
                counter += 1
            syn = self.sabio_session.query(sabio_rk.entry_synonym).filter_by(entry__id = item._id).all()
            res = self.sabio_session.query(sabio_rk.entry_resource).filter_by(entry__id = item._id).all()
            for synonyms in syn:
                syns = self.sabio_session.query(sabio_rk.Synonym).get(synonyms.synonym__id)
                metadata.synonym.append(self.get_or_create_object(Synonym, name = syns.name))
            for docs in res:
                resource = self.sabio_session.query(sabio_rk.Resource).get(docs.resource__id)
                if resource.namespace == 'taxonomy':
                    metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = resource.id))
                elif resource.namespace == 'uniprot':
                    uniprot = resource.id
                else:
                    metadata.resource.append(self.get_or_create_object(Resource, namespace = resource.namespace,
                        _id = resource.id))
            if item._type == 'compound':
                structure = self.sabio_session.query(sabio_rk.compound_compound_structure).filter_by(compound__id = item._id).all()
                if structure:
                    for shape in structure:
                        struct = self.sabio_session.query(sabio_rk.CompoundStructure).get(shape.compound_structure__id)
                        if struct.format == 'smiles':
                            obs_pp.structure = self.get_or_create_object(Structure, type = 'Structure', name = item.name,
                            _value_smiles = struct.value, _value_inchi = struct._value_inchi,
                            _structure_formula_connectivity = struct._value_inchi_formula_connectivity, _metadata = metadata)
                            break
                else:
                    obs_pp.structure = None
                obs_pe.compound = self.get_or_create_object(Compound, type = 'Compound',
                    name = item.name, compound_name = item.name,
                    _is_name_ambiguous = self.sabio_session.query(sabio_rk.Compound).get(item._id)._is_name_ambiguous,
                    structure = obs_pp.structure, _metadata = metadata)
            elif item._type == 'enzyme':
                complx = self.sabio_session.query(sabio_rk.Enzyme).get(item._id)
                obs_pe.protein_complex = self.get_or_create_object(ProteinComplex, type = 'Enzyme' ,
                    name = item.name , complex_name = item.name,
                    molecular_weight = complx.molecular_weight, funcat_dsc = 'Enzyme', _metadata = metadata)
            elif item._type == 'enzyme_subunit':
                subunit = self.sabio_session.query(sabio_rk.EnzymeSubunit).get(item._id)
                complx = self.sabio_session.query(sabio_rk.Entry).get(subunit.enzyme_id)
                result = self.session.query(ProteinComplex).filter_by(complex_name = complx.name).first()
                obs_pe.protein_subunit = self.get_or_create_object(ProteinSubunit, type = 'Enzyme Subunit',
                    name = item.name, subunit_name = item.name,
                    uniprot_id = uniprot, coefficient = subunit.coefficient,
                    sequence = subunit.coefficient, molecular_weight = subunit.molecular_weight,
                    proteincomplex = result, _metadata = metadata)
            elif item._type == 'kinetic_law':
                res = self.sabio_session.query(sabio_rk.kinetic_law_resource).filter_by(kinetic_law__id = item._id).all()
                law = self.sabio_session.query(sabio_rk.KineticLaw).get(item._id)
                entry = self.sabio_session.query(sabio_rk.Entry).get(law.enzyme_id)
                result = self.session.query(ProteinComplex).filter_by(complex_name = entry.name).first()
                for docs in res:
                    resource = self.sabio_session.query(sabio_rk.Resource).get(docs.resource__id)
                    metadata.resource.append(self.get_or_create_object(Resource, namespace = resource.namespace, _id = resource.id))
                metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = law.taxon))
                metadata.cell_line.append(self.get_or_create_object(CellLine, name = law.taxon_variant))
                metadata.conditions.append(self.get_or_create_object(Conditions, temperature = law.temperature, ph = law.ph, media = law.media))
                obs_pp.kinetic_law = self.get_or_create_object(KineticLaw, type = 'Kinetic Law', enzyme = result,
                    enzyme_type = law.enzyme_type, tissue = law.tissue, mechanism = law.mechanism, equation = law.equation, _metadata = metadata)
                rxn = self.sabio_session.query(sabio_rk.ReactionParticipant).filter(or_(sabio_rk.ReactionParticipant.reactant_kinetic_law_id == law._id, sabio_rk.ReactionParticipant.product_kinetic_law_id == law._id, sabio_rk.ReactionParticipant.modifier_kinetic_law_id == law._id)).all()
                for row in rxn:
                    compound_name = self.sabio_session.query(sabio_rk.Entry).get(row.compound_id).name
                    _compound = self.session.query(Compound).filter_by(compound_name = compound_name).first()
                    if row.reactant_kinetic_law_id:
                        reaction = Reaction(compound = _compound,
                            coefficient = row.coefficient, _is_reactant = 1, rxn_type = row.type,
                            kinetic_law_id = obs_pp.kinetic_law.kineticlaw_id)
                    elif row.product_kinetic_law_id:
                        reaction = Reaction( compound = _compound,
                            coefficient = row.coefficient, _is_product = 1, rxn_type = row.type,
                            kinetic_law_id = obs_pp.kinetic_law.kineticlaw_id)
                    elif row.modifier_kinetic_law_id:
                        reaction = Reaction(compound = _compound,
                            coefficient = row.coefficient, _is_modifier = 1, rxn_type = row.type,
                            kinetic_law_id = obs_pp.kinetic_law.kineticlaw_id)
                total_param = self.sabio_session.query(sabio_rk.Parameter).filter_by(kinetic_law_id = law._id).all()
                for param in total_param:
                    if param.compound_id:
                        compound_name = self.sabio_session.query(sabio_rk.Entry).get(param.compound_id).name
                        _compound = self.session.query(Compound).filter_by(compound_name = compound_name).first()
                        parameter = Parameter(sabio_type = param.type, value = param.value, error = param.error,
                            units = param.units, observed_name = param.observed_name, kinetic_law = obs_pp.kinetic_law,
                            observed_sabio_type = param.observed_type, observed_value = param.observed_value, compound = _compound,
                            observed_error = param.observed_error, observed_units = param.observed_units)
                    else:
                        parameter = Parameter(sabio_type = param.type, value = param.value, error = param.error,
                            units = param.units, observed_name = param.observed_name, kinetic_law = obs_pp.kinetic_law,
                            observed_sabio_type = param.observed_type, observed_value = param.observed_value,
                            observed_error = param.observed_error, observed_units = param.observed_units)

    # def add_arrayexpressdb(self):
        # arrayexpressdb = array_express.ArrayExpress(name = 'array_express', load_content=False, download_backup=False, max_entries = 5)
        # arrayexpressdb.load_experiments_from_text()
        # array_session = arrayexpressdb.session
