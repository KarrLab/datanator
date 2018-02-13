from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate
from kinetic_datanator.config import config


app = Flask(__name__)
app.config.from_object(config.Config)
db = SQLAlchemy(app)
migrate = Migrate(app, db)


class Observation(db.Model):
    """
    Represents an Observation of a Physical Entity or Property in the Common Schema
    Attributes:
        id (:obj:`int`): Common Schema Observation Identifier
        _metadata_id (:obj:`int` of :obj:`Metadata`): Related Metadata ID
    """

    __tablename__ = 'observation'

    id = db.Column(db.Integer, primary_key=True)

    _metadata_id = db.Column(db.Integer, db.ForeignKey('_metadata.id'))
    _metadata = db.relationship('Metadata', backref='observation')

    __mapper_args__ = {'polymorphic_identity': 'observation'}


_metadata_taxon = db.Table(
    '_metadata_taxon', db.Model.metadata,
    db.Column('_metadata_id', db.Integer,
              db.ForeignKey('_metadata.id'), index=True),
    db.Column('taxon_id', db.Integer, db.ForeignKey(
        'taxon.ncbi_id'), index=True),
)
# :obj:`db.Table`: Metadata:Taxon many-to-many association table

_metadata_method = db.Table(
    '_metadata_method', db.Model.metadata,
    db.Column('_metadata_id', db.Integer,
              db.ForeignKey('_metadata.id'), index=True),
    db.Column('method_id', db.Integer, db.ForeignKey('method.id'), index=True),
)
# :obj:`db.Table`: Metadata:Method many-to-many association table

_metadata_resource = db.Table(
    '_metadata_resource', db.Model.metadata,
    db.Column('_metadata_id', db.Integer,
              db.ForeignKey('_metadata.id'), index=True),
    db.Column('resource_id', db.Integer,
              db.ForeignKey('resource.id'), index=True),
)
# :obj:`db.Table`: Metadata:Resource many-to-many association table

_metadata_cell_line = db.Table(
    '_metadata_cell_line', db.Model.metadata,
    db.Column('_metadata_id', db.Integer,
              db.ForeignKey('_metadata.id'), index=True),
    db.Column('cell_line_id', db.Integer,
              db.ForeignKey('cell_line.id'), index=True),
)
# :obj:`db.Table`: Metadata:CellLine many-to-many association table

_metadata_synonym = db.Table(
    '_metadata_synonym', db.Model.metadata,
    db.Column('_metadata_id', db.Integer,
              db.ForeignKey('_metadata.id'), index=True),
    db.Column('synonym_id', db.Integer,
              db.ForeignKey('synonym.id'), index=True),
)
# :obj:`db.Table`: Metadata:Synonym many-to-many association table

_metadata_conditions = db.Table(
    '_metadata_conditions', db.Model.metadata,
    db.Column('_metadata_id', db.Integer,
              db.ForeignKey('_metadata.id'), index=True),
    db.Column('conditions_id', db.Integer,
              db.ForeignKey('conditions.id'), index=True),
)
# :obj:`db.Table`: Metadata:Conditions many-to-many association table

_metadata_compartment = db.Table(
    '_metadata_compartment', db.Model.metadata,
    db.Column('_metadata_id', db.Integer,
              db.ForeignKey('_metadata.id'), index=True),
    db.Column('compartment_id', db.Integer, db.ForeignKey(
        'cell_compartment.id'), index=True),
)
# :obj:`db.Table`: Metadata:Conditions many-to-many association table

_metadata_characteristic = db.Table(
    '_metadata_characteristic', db.Model.metadata,
    db.Column('_metadata_id', db.Integer,
              db.ForeignKey('_metadata.id'), index=True),
    db.Column('characteristic_id', db.Integer, db.ForeignKey(
        'characteristic.id'), index=True),
)
# :obj:`db.Table`: Metadata:Conditions many-to-many association table

_metadata_variable = db.Table(
    '_metadata_variable', db.Model.metadata,
    db.Column('_metadata_id', db.Integer,
              db.ForeignKey('_metadata.id'), index=True),
    db.Column('variable_id', db.Integer, db.ForeignKey(
        'variable.id'), index=True),
)
# :obj:`db.Table`: Metadata:Conditions many-to-many association table

subunit_interaction = db.Table(
    'subunit_interaction', db.Model.metadata,
    db.Column('protein_subunit_id', db.Integer, db.ForeignKey(
        'protein_subunit.subunit_id'), index=True),
    db.Column('interaction_id', db.Integer, db.ForeignKey(
        'protein_interactions.interaction_id'), index=True)
)

rnaseqdataset_rnaseqexperiment = db.Table(
    'rnaseqdataset_rnaseqexperiment', db.Model.metadata,
    db.Column('experiment_id', db.Integer, db.ForeignKey(
        'rna_seq_experiment.experiment_id'), index=True),
    db.Column('sample_id', db.Integer, db.ForeignKey(
        'rna_seq_dataset.sample_id'), index=True)
)


class Metadata(db.Model):
    """
    db.Table representing Metadata identifiers for entities and properties
    Attributes:
        name (:obj:`str`): Name of the entity or property
    """
    __tablename__ = '_metadata'

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(255), unique=True)

    description = db.Column(db.Text())
    taxon = db.relationship(
        'Taxon', secondary=_metadata_taxon, backref='_metadata')
    method = db.relationship(
        'Method', secondary=_metadata_method, backref='_metadata')
    resource = db.relationship(
        'Resource', secondary=_metadata_resource, backref='_metadata')
    cell_line = db.relationship(
        'CellLine', secondary=_metadata_cell_line, backref='_metadata')
    synonym = db.relationship(
        'Synonym', secondary=_metadata_synonym, backref='_metadata')
    conditions = db.relationship(
        'Conditions', secondary=_metadata_conditions, backref='_metadata')
    cell_compartment = db.relationship(
        'CellCompartment', secondary=_metadata_compartment, backref='_metadata')
    characteristic = db.relationship(
        'Characteristic', secondary=_metadata_characteristic, backref='_metadata')
    variable = db.relationship(
        'Variable', secondary=_metadata_variable, backref='_metadata')


class Method(db.Model):
    """
    Represents the method of collection for a given entity or Property
    Attributes:
        name (:obj:`str`): Name of the Method
        comments (:obj:`str`): Comments on the method
    """
    __tablename__ = 'method'
    __searchable__ = ['name']

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(255))
    comments = db.Column(db.String(255))
    performer = db.Column(db.String(255))
    hardware = db.Column(db.String(255))
    software = db.Column(db.String(255))

class Characteristic(db.Model):
    """
    Represents the method of collection for a given entity or Property
    Attributes:
        name (:obj:`str`): Name of the Method
        comments (:obj:`str`): Comments on the method
    """
    __tablename__ = 'characteristic'
    #__searchable__ = ['name']

    id = db.Column(db.Integer, primary_key=True)
    category = db.Column(db.String(255))
    value = db.Column(db.String(255))

class Variable(db.Model):
    """
    Represents the method of collection for a given entity or Property
    Attributes:
        name (:obj:`str`): Name of the variable
        comments (:obj:`str`): Comments on the method
    """
    __tablename__ = 'variable'
    #__searchable__ = ['name']

    id = db.Column(db.Integer, primary_key=True)
    category = db.Column(db.String(255))
    value = db.Column(db.String(255))
    units = db.Column(db.String(255))


class Taxon(db.Model):
    """
    Represents the species of a given physical entity or property
    Attributes:
        ncbi_id (:obj:`int`): NCBI id of the species
        name (:obj:`str`): Name of the species
    """
    __tablename__ = 'taxon'
    __searchable__ = ['ncbi_id', 'name']

    ncbi_id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(255))


class Synonym(db.Model):
    """
    Represents a synonym of a given physical entity or property
    Attributes:
        name (:obj:`str`): Name of the Synonym
    """

    __searchable__ = ['name']

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(255))

    __tablename__ = 'synonym'


class Resource(db.Model):
    """
    Represents a resource of a given physical entity or property
    Attributes:
        namespace (:obj:`str`): Name of the classifier of the resource (Ex. Pubmed)
        _id (:obj:`str`): Identifier of the resource
        release_date(:obj:`str`): The date that resource released the data
    """
    id = db.Column(db.Integer, primary_key=True)
    namespace = db.Column(db.String(255))
    _id = db.Column(db.String(255))
    release_date = db.Column(db.String(255))

    __tablename__ = 'resource'


class CellLine(db.Model):
    """
    Represents a cell line of a given physical entity or property
    Attributes:
        name (:obj:`str`): Name of the Cell Line
    """
    __tablename__ = 'cell_line'
    __searchable__ = ['name']

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(255))


class Conditions(db.Model):
    """
    Represents the conditions of a given physical entity or property
    Attributes:
        growth_status (:obj:`str`): Type of growth status
        media (:obj:`str`): Media composition
        temperature (:obj:`float`): Temperature of the sample (C)
        ph (:obj:`float`): pH of the sample
        growth_system (:obj:`str`): Type of growth system
    """

    id = db.Column(db.Integer, primary_key=True)
    growth_status = db.Column(db.String(255))
    media = db.Column(db.String(255))
    temperature = db.Column(db.Float)
    ph = db.Column(db.Float)
    growth_system = db.Column(db.String(255))

    __tablename__ = 'conditions'


class CellCompartment(db.Model):
    """
    Represents a cell compartment of a given physical entity or property
    Ties especially to the reacitons because this is where the reactions occur
    Attributes:
        name (:obj:`str`): Name of the Cell Compartment
    """
    __tablename__ = 'cell_compartment'
    __searchable__ = ['name']

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(255), unique=True, index=True)


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

    observation_id = db.Column(db.Integer, db.ForeignKey(
        'observation.id'), primary_key=True)
    type = db.Column(db.String(255))
    name = db.Column(db.String(255))


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
    __searchable__ = ['subunit_name', 'uniprot_id', 'gene_name', 'canonical_sequence']
    __mapper_args__ = {'polymorphic_identity': 'protein_subunit'}

    subunit_id = db.Column(db.Integer, db.ForeignKey(
        'physical_entity.observation_id'), primary_key=True)
    subunit_name = db.Column(db.String(255))
    uniprot_id = db.Column(db.String(255))
    entrez_id = db.Column(db.Integer)
    ec_number = db.Column(db.String(255))
    gene_name = db.Column(db.String(255))
    gene_syn = db.Column(db.String(255))
    class_name = db.Column(db.String(255))
    family_name = db.Column(db.String(255))
    coefficient = db.Column(db.Integer)
    canonical_sequence = db.Column(db.String(255))
    mass = db.Column(db.Integer)
    length = db.Column(db.Integer)
    molecular_weight = db.Column(db.Float)
    pax_load = db.Column(db.Integer)

    interaction = db.relationship(
        'ProteinInteractions', secondary=subunit_interaction, backref='protein_subunit')

    proteincomplex_id = db.Column(
        db.Integer, db.ForeignKey('protein_complex.complex_id'))
    proteincomplex = db.relationship(
        'ProteinComplex', backref='protein_subunit', foreign_keys=[proteincomplex_id])

    def __name__(self):
        return 'ProteinSubunit'



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
    __searchable__ = ['complex_name', 'go_id', 'funcat_id', 'su_cmt']
    __mapper_args__ = {'polymorphic_identity': 'protein_complex'}

    complex_id = db.Column(db.Integer, db.ForeignKey(
        'physical_entity.observation_id'), primary_key=True)
    complex_name = db.Column(db.String(255))
    go_id = db.Column(db.String(255))
    go_dsc = db.Column(db.String(255))
    funcat_id = db.Column(db.String(255))
    funcat_dsc = db.Column(db.String(255))
    su_cmt = db.Column(db.String(255))
    complex_cmt = db.Column(db.String(255))
    disease_cmt = db.Column(db.String(255))
    class_name = db.Column(db.String(255))
    family_name = db.Column(db.String(255))
    molecular_weight = db.Column(db.Float)

    def __name__(self):
        return 'ProteinComplex'


class Compound(PhysicalEntity):
    """
    Represents a Compound - An instance of Physical Entity
    Attributes:
        compound_id (:obj:`int`): Common Schema Observation Identifier
        compound_name (:obj:`str`): Name of the Compound
        description (:obj:`str`):
        comment = db.Column(db.String(255))
        _is_name_ambiguous = db.Column(db.Boolean)
    """
    __tablename__ = 'compound'
    __searchable__ = ['compound_name', 'description']
    __mapper_args__ = {'polymorphic_identity': 'compound'}

    compound_id = db.Column(db.Integer, db.ForeignKey(
        'physical_entity.observation_id'), primary_key=True)
    compound_name = db.Column(db.String(255))
    description = db.Column(db.String(255))
    comment = db.Column(db.String(255))
    _is_name_ambiguous = db.Column(db.Boolean)

    structure_id = db.Column(db.Integer, db.ForeignKey('structure.struct_id'))
    structure = db.relationship('Structure', backref='compound')

    def __name__(self):
        return 'Compound'


class PhysicalProperty(Observation):
    """
    Represents a Physical Property in the Common Schema
    Attributes:
        observation_id (:obj:`int`): Common Schema Observation Identifier
        type (:obj:`str`): Type of Physical Property (Ex. Concentration)
        name (:obj:`str`): Name of the Physical Property
    """
    observation_id = db.Column(db.Integer, db.ForeignKey(
        'observation.id'), primary_key=True)
    type = db.Column(db.String(255))
    name = db.Column(db.String(255))

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

    struct_id = db.Column(db.Integer, db.ForeignKey(
        'physical_property.observation_id'), primary_key=True)
    _value_smiles = db.Column(db.String(255))
    _value_inchi = db.Column(db.String(255))
    _structure_formula_connectivity = db.Column(db.String(255))


class Concentration(PhysicalProperty):
    """
    Represents the concentration of an entity
    Attributes:
        value (:obj:`float`): concentration of a tagged compound
        error (:obj:`float`): uncertainty of corresponding concentration value
    """

    __tablename__ = 'concentration'
    __mapper_args__ = {'polymorphic_identity': 'concentration'}

    concentration_id = db.Column(db.Integer, db.ForeignKey(
        'physical_property.observation_id'), primary_key=True)

    compound_id = db.Column(db.Integer, db.ForeignKey('compound.compound_id'))
    compound = db.relationship('Compound', backref='concentration')

    value = db.Column(db.Float)
    error = db.Column(db.Float)


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

    kinetic_law_id = db.Column(db.Integer, db.ForeignKey(
        'physical_property.observation_id'), primary_key=True)

    enzyme_id = db.Column(db.Integer, db.ForeignKey(
        'protein_complex.complex_id'), index=True)
    enzyme = db.relationship(ProteinComplex, backref='kinetic_law')

    enzyme_type = db.Column(db.String(255))
    tissue = db.Column(db.String(255))
    mechanism = db.Column(db.String(255))
    equation = db.Column(db.String(255))


class Reaction(db.Model):
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

    reaction_id = db.Column(db.Integer, primary_key=True)
    compound_id = db.Column(db.Integer, db.ForeignKey('compound.compound_id'))
    compound = db.relationship(Compound, backref='reaction')
    compartment_id = db.Column(
        db.Integer, db.ForeignKey('cell_compartment.id'))
    compartment = db.relationship(CellCompartment, backref='reaction')
    coefficient = db.Column(db.Float)
    _is_reactant = db.Column(db.Boolean)
    _is_product = db.Column(db.Boolean)
    _is_modifier = db.Column(db.Boolean)
    rxn_type = db.Column(db.String(255))

    kinetic_law_id = db.Column(
        db.Integer, db.ForeignKey('kinetic_law.kinetic_law_id'))
    kinetic_law = db.relationship(KineticLaw, backref='reaction')


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

    dataset_id = db.Column(db.Integer, db.ForeignKey(
        'physical_property.observation_id'), primary_key=True)
    file_name = db.Column(db.String, unique=True)
    score = db.Column(db.Float)
    weight = db.Column(db.Integer)
    coverage = db.Column(db.Integer)

class RNASeqDataSet(PhysicalProperty):
    __tablename__ = 'rna_seq_dataset'
    __mapper_args__ = {'polymorphic_identity': 'rna_seq_dataset'}

    sample_id = db.Column(db.Integer, db.ForeignKey(
        'physical_property.observation_id'), primary_key=True)
    experiment_accession_number = db.Column(db.String)
    sample_name = db.Column(db.String)
    assay = db.Column(db.String)
    ensembl_organism_strain = db.Column(db.String)
    read_type = db.Column(db.String)
    full_strain_specificity = db.Column(db.Boolean)

class RNASeqExperiment(PhysicalProperty):
    __tablename__ = 'rna_seq_experiment'
    __mapper_args__ = {'polymorphic_identity': 'rna_seq_experiment'}

    experiment_id = db.Column(db.Integer, db.ForeignKey(
        'physical_property.observation_id'), primary_key=True)
    samples = db.relationship(
        'RNASeqDataSet', secondary= rnaseqdataset_rnaseqexperiment, backref='experiment')
    accesion_number = db.Column(db.String)
    exp_name = db.Column(db.String)
    has_fastq_files = db.Column(db.Boolean)


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

    dataset_id = db.Column(db.Integer, db.ForeignKey(
        'physical_property.observation_id'), primary_key=True)
    version = db.Column(db.Integer)

    complex_id = db.Column(db.Integer, db.ForeignKey(
        'protein_complex.complex_id'))
    tf = db.relationship('ProteinComplex', backref='dna_binding_dataset')

    subunit_id = db.Column(db.Integer, db.ForeignKey(
        'protein_subunit.subunit_id'))
    subunit = db.relationship('ProteinSubunit', backref='dna_binding_dataset')


class Parameter(db.Model):
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

    parameter_id = db.Column(db.Integer, primary_key=True)

    kinetic_law_id = db.Column(
        db.Integer, db.ForeignKey('kinetic_law.kinetic_law_id'))
    kinetic_law = db.relationship(KineticLaw, backref='parameter')

    sabio_type = db.Column(db.Integer, index=True)
    compound_id = db.Column(db.Integer, db.ForeignKey(
        'compound.compound_id'), index=True)
    compound = db.relationship(Compound, backref='parameter')

    value = db.Column(db.Float)
    error = db.Column(db.Float)
    units = db.Column(db.String(255), index=True)

    observed_name = db.Column(db.String(255))
    observed_sabio_type = db.Column(db.Integer)
    observed_value = db.Column(db.Float)
    observed_error = db.Column(db.Float)
    observed_units = db.Column(db.String(255))


class AbundanceData(db.Model):
    """
    Represents protein abundance data from the Pax DB database
    Attributes:
        abundance (:obj:`float`): Represents protein abundance from given observation in ppm
        dataset_id  (:obj:`int`): Represents the dataset from which the abundance stems from
        subunit_id  (:obj:`int`): Represents the protein frmo which the abundance stems from
    """

    __tablename__ = 'abundance_data'

    abundance_id = db.Column(db.Integer, primary_key=True)
    abundance = db.Column(db.Float)

    dataset_id = db.Column(db.Integer, db.ForeignKey(
        'abundance_dataset.dataset_id'))
    dataset = db.relationship(
        'AbundanceDataSet', backref='abundance_data', foreign_keys=[dataset_id])

    subunit_id = db.Column(db.Integer, db.ForeignKey(
        'protein_subunit.subunit_id'), index=True)
    subunit = db.relationship('ProteinSubunit')

    pax_load = db.Column(db.Integer)
    uniprot_id = db.Column(db.Integer)


class DNABindingData(db.Model):
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

    position_id = db.Column(db.Integer, primary_key=True)
    position = db.Column(db.Integer, index=True)
    frequency_a = db.Column(db.Integer)
    frequency_c = db.Column(db.Integer)
    frequency_g = db.Column(db.Integer)
    frequency_t = db.Column(db.Integer)
    jaspar_id = db.Column(db.Integer)

    dataset_id = db.Column(db.Integer, db.ForeignKey(
        'dna_binding_dataset.dataset_id'))
    dataset = db.relationship(
        'DNABindingDataset', backref='dna_binding_data', foreign_keys=[dataset_id])


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
    __searchable__ = ['participant_a', 'participant_b']
    __mapper_args__ = {'polymorphic_identity': 'protein_interactions'}

    interaction_id = db.Column(db.Integer, db.ForeignKey(
        'physical_property.observation_id'), primary_key=True)
    participant_a = db.Column(db.String(255))
    participant_b = db.Column(db.String(255))
    interaction = db.Column(db.String(255))
    site_a = db.Column(db.String(255))
    site_b = db.Column(db.String(255))
    stoich_a = db.Column(db.String(255))
    stoich_b = db.Column(db.String(255))
    interaction_type = db.Column(db.String(255))
    publication = db.Column(db.String(255))


class Progress(db.Model):
    """
    Represents amount loaded of large DBs (Ex. Pax and Sabio)
    Attributes:
        database_name (:obj:`str`): Name of observed databse
        amount_loaded (:obj:`int`): Amount of entries loaded in Common Schema
    """
    __tablename__ = 'progress'

    database_name = db.Column(db.String, primary_key=True)
    amount_loaded = db.Column(db.Integer)
