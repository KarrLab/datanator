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
    Column('resource_id', Integer, ForeignKey('resource._id'), index=True),
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

    """

    __tablename__ = 'observation'

    id = Column(Integer, primary_key = True)

    _metadata_id = Column(Integer, ForeignKey('_metadata.id'))
    _metadata = relationship('Metadata', backref = 'observation')


class PhysicalEntity(Observation):
    """

    """

    __tablename__ = 'physical_entity'
    __mapper_args__ = {'polymorphic_identity': 'physical_entity'}

    observation_id = Column(Integer, ForeignKey('observation.id'), primary_key = True)
    type = Column(String(255))
    name = Column(String(255))

class ProteinSubunit(PhysicalEntity):
    """

    """
    __tablename__ = 'protein_subunit'
    __mapper_args__ = {'polymorphic_identity': 'protein_subunit'}

    subunit_id = Column(Integer, ForeignKey('physical_entity.observation_id'), primary_key = True)
    subunit_name = Column(String(255))
    uniprot_id = Column(String(255))
    entrez_id = Column(Integer)
    ensembl_id = Column(String(255), unique = True)
    gene_name = Column(String(255))
    gene_syn  = Column(String(255))
    class_name = Column(String(255))
    family_name = Column(String(255))
    sequence = Column(String(255), unique = True)
    coefficient = Column(Integer)
    sequence = Column(String(255))
    molecular_weight = Column(Float)

    proteincomplex_id = Column(Integer, ForeignKey('protein_complex.complex_id'))
    proteincomplex = relationship('ProteinComplex', backref = 'protein_subunit', foreign_keys=[proteincomplex_id])

class ProteinComplex(PhysicalEntity):
    """

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


class Compound(PhysicalEntity):
    """

    """
    __tablename__ = 'compound'
    __mapper_args__ = {'polymorphic_identity': 'compound'}

    compound_id = Column(Integer, ForeignKey('physical_entity.observation_id'), primary_key = True)
    compound_name = Column(String(255), unique = True)
    description = Column(String(255))
    comment = Column(String(255))

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
    structure = Column(String(255))
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

        session = self.session

        paxdb = pax.Pax(cache_dirname = self.cache_dirname, clear_content = True,  load_content=False, download_backup=False, max_entries = self.max_entries)
        # paxdb = pax.Pax(name = 'pax', clear_content = True,  load_content=False, download_backup=False, max_entries = self.max_entries)
        paxdb.load_content()
        self.pax_session = paxdb.session

        corumdb = corum.Corum(cache_dirname = self.cache_dirname , clear_content = True, load_content=False, download_backup=False, max_entries = self.max_entries)
        # corumdb = corum.Corum(name = 'corum', clear_content = True, load_content=False, download_backup=False, max_entries = self.max_entries)
        corumdb.load_content()
        self.corum_session = corumdb.session

        jaspardb = jaspar.Jaspar(cache_dirname = self.cache_dirname, clear_content = True, load_content=False, download_backup=False, max_entries = self.max_entries)
        # jaspardb = jaspar.Jaspar(name = 'jaspar', clear_content = True, load_content=False, download_backup=False, max_entries = self.max_entries)
        jaspardb.load_content()
        self.jaspar_session = jaspardb.session

        ecmDB = ecmdb.Ecmdb(cache_dirname = self.cache_dirname, clear_content = True, load_content=False, download_backup=False, max_entries = self.max_entries)
        # ecmDB = ecmdb.Ecmdb(name = 'ecmdb', clear_content = True, load_content=False, download_backup=False, max_entries = self.max_entries)
        ecmDB.load_content()
        self.ecmdb_session = ecmDB.session

        # # # sabiodb = sabio_rk.SabioRk(cache_dirname = self.cache_dirname, clear_content = True, load_content=False, download_backup=False, max_entries = self.max_entries)
        # sabiodb = sabio_rk.SabioRk(name = 'sabio', clear_content = True, load_content=False, download_backup=False, max_entries = 10)
        # sabiodb.load_content()
        # self.sabio_session = sabiodb.session

        """ NOT WORKING """
        #TODO: AE needs 1. normalization to schema(max_entries implementation/working functions) 2. working tests

        # arrayexpressdb = array_express.ArrayExpress(name = 'array_express', load_content=False, download_backup=False, max_entries = 5)
        # arrayexpressdb.load_experiments_from_text()
        # array_session = arrayexpressdb.session

        self.add_content_to_db()

    def add_content_to_db(self):
        observation = Observation()
        observation.physical_entity = PhysicalEntity()
        observation.physical_property = PhysicalProperty()

        # # Pax --- DONE
        pax_dataset = self.pax_session.query(pax.Dataset).all()

        for item in pax_dataset:
            metadata = self.get_or_create_object(Metadata, name = item.file_name)
            metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = item.taxon_ncbi_id,
                name = self.pax_session.query(pax.Taxon).get(item.taxon_ncbi_id).species_name))
            metadata.resource.append(self.get_or_create_object(Resource, namespace = 'url', _id = item.publication))
            observation.physical_property.abundance_dataset = self.get_or_create_object(AbundanceDataSet, type = 'Protein Abundance Dataset',
                name = item.file_name, file_name = item.file_name, score = item.score, weight = item.weight, coverage= item.coverage, _metadata = metadata )
            abundance = self.pax_session.query(pax.Observation).filter_by(dataset_id = item.id).all()
            for data in abundance:
                ## TODO: May have to do some IF exists in this. Also give names to the physical entities
                uniprot_id = self.pax_session.query(pax.Protein).get(data.protein_id).uniprot_id
                observation.physical_entity.protein_subunit = self.get_or_create_object(ProteinSubunit,
                    type = 'Protein Subunit', name = uniprot_id, uniprot_id = uniprot_id,
                    _metadata = metadata)
                abundance_data = (self.get_or_create_object(AbundanceData,
                    abundance = data.abundance, dataset = observation.physical_property.abundance_dataset,
                    subunit = observation.physical_entity.protein_subunit))

        # # Corum --- DONE
        corum_complex = self.corum_session.query(corum.Complex).all()
        corum_subunit = self.corum_session.query(corum.Subunit).all()

        for item in corum_complex:
            obs = self.corum_session.query(corum.Observation).get(item.observation_id)
            metadata = self.get_or_create_object(Metadata, name = item.complex_name )
            metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = obs.taxon_ncbi_id,
                name = self.corum_session.query(corum.Taxon).get(obs.taxon_ncbi_id).swissprot_id))
            metadata.resource.append(self.get_or_create_object(Resource, namespace = 'pubmed', _id = obs.pubmed_id))
            metadata.method.append(self.get_or_create_object(Method, name = 'purification', comments = obs.pur_method))
            metadata.cell_line.append(self.get_or_create_object(CellLine, name = obs.cell_line))
            observation.physical_entity.protein_complex = self.get_or_create_object(ProteinComplex,
                type = 'Protein Complex', name = item.complex_name, complex_name = item.complex_name, go_id = item.go_id,
                go_dsc = item.go_dsc, funcat_id = item.funcat_id, funcat_dsc = item.funcat_dsc, su_cmt = item.su_cmt,
                complex_cmt = item.complex_cmt, disease_cmt = item.disease_cmt,  _metadata = metadata)

        for item in corum_subunit:
            cmplx = self.corum_session.query(corum.Complex).get(item.complex_id)
            result = self.session.query(ProteinComplex).filter_by(complex_name = cmplx.complex_name).first()
            obs = self.session.query(Observation).get(result.complex_id)
            observation.physical_entity.protein_subunit = self.get_or_create_object(ProteinSubunit,
                type = 'Protein Subunit', uniprot_id = item.su_uniprot,
                entrez_id = item.su_entrezs, name = item.protein_name, subunit_name = item.protein_name, gene_name=item.gene_name,
                gene_syn = item.gene_syn, proteincomplex = result, _metadata = self.session.query(Metadata).get(obs._metadata_id))


        # Jaspar -- DONE
        jaspar_matrix = self.jaspar_session.query(jaspar.Matrix).all()
        jaspar_matrixposition = self.jaspar_session.query(jaspar.MatrixPosition).all()
        jaspar_transcription_factor = self.jaspar_session.query(jaspar.TranscriptionFactor).all()

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
                        observation.physical_entity.protein_complex = self.get_or_create_object(ProteinComplex, type = 'Transcription Factor Complex', name = temp_tf.name,
                            complex_name = temp_tf.name, su_cmt = dialoge, complex_cmt = 'transcription factor', class_name = class_name, family_name = fam_name, _metadata = metadata)
                        observation.physical_property.dna_binding_dataset = self.get_or_create_object(DNABindingDataset, type = 'DNA Binding Dataset', name = temp_tf.name,
                            version = item.version, tf = observation.physical_entity.protein_complex, _metadata = metadata)
                        for pos in jaspar_matrixposition:
                            data = self.get_or_create_object(DNABindingData, position = pos.position,
                            frequency_a = pos.frequency_a, frequency_c = pos.frequency_c, frequency_g = pos.frequency_g,
                            frequency_t = pos.frequency_t, dataset = observation.physical_property.dna_binding_dataset)
                    else:
                        sub_id = self.jaspar_session.query(jaspar.transcription_factor_subunit).filter_by(transcription_factor__id = temp_tf._id).first().subunit_id
                        uniprot_id = self.jaspar_session.query(jaspar.Subunit).get(sub_id).uniprot_id
                        metadata = self.get_or_create_object(Metadata, name = temp_tf.name)
                        metadata.taxon.append(self.get_or_create_object(Taxon, ncbi_id = ncbi))
                        metadata.method.append(self.get_or_create_object(Method, name = type_name))
                        metadata.resource.append(self.get_or_create_object(Resource, namespace = 'pubmed', _id = doc))
                        observation.physical_entity.protein_subunit = self.get_or_create_object(ProteinSubunit, uniprot_id = uniprot_id, type = 'Transcription Factor Subunit',
                            name = temp_tf.name, subunit_name = temp_tf.name, gene_name = temp_tf.name,
                            class_name = class_name, family_name = fam_name, _metadata = metadata)
                        observation.physical_property.dna_binding_dataset = self.get_or_create_object(DNABindingDataset, type = 'DNA Binding Dataset', name = temp_tf.name,
                            version = item.version, subunit = observation.physical_entity.protein_subunit, _metadata = metadata)
                        for pos in jaspar_matrixposition:
                            data = self.get_or_create_object(DNABindingData, position = pos.position,
                            frequency_a = pos.frequency_a, frequency_c = pos.frequency_c, frequency_g = pos.frequency_g,
                            frequency_t = pos.frequency_t, dataset = observation.physical_property.dna_binding_dataset)



        #ECMDB -- DONE
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
                    observation.physical_property.concentration = self.get_or_create_object(Concentration, type = 'Concentration', name = item.name,
                        value = rows.value, error = rows.error, _metadata = metadata)
            observation.physical_property.structure = self.get_or_create_object(Structure, type = 'Structure', name = item.name,
            structure = item.structure, _structure_formula_connectivity = item._structure_formula_connectivity, _metadata = metadata)
            observation.physical_entity.compound = self.get_or_create_object(Compound, type = 'Compound', name = item.name,
                compound_name = item.name, description = item.description, comment = item.comment, structure = observation.physical_property.structure,
                _metadata = metadata)

        # # SabioRk
        # entry = self.sabio_session.query(sabio_rk.Entry).all()
        #
        # for item in sabio_entry:
        #     if item._type == 'compound':
        #         structure = self.sabio_session.query(sabio_rk.compound_compound_structure).filter_by(item._id).all()
        #         for shape in structure:
        #             struct = self.sabio_session.query(sabio_rk.CompoundStructure).get(shape.compound_structure__id)
        #             compound = self.get_or_create_object(Compound, name = item.name,
        #                 _is_name_ambiguous = self.sabio_session.query(sabio_rk.Compound).get(item._id)._is_name_ambiguous)
        #             compound.structure = self.get_or_create_object(CompoundStructure, value = struct.value,
        #                 format = struct.format, _structure_formula_connectivity = struct._value_inchi_formula_connectivity)
        #             compound.biomolecule = self.get_or_create_object(BioMolecule, type = item._type, name = item.name)


        if self.verbose:
            print('Comitting')
        self.session.commit()
