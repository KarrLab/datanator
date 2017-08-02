"""
This module downloads the JASPAR database of transcription factor binding motifs (http://jaspar.genereg.net/)
via a seris of text files, parses them, and stores them in an SQLlite database.

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-08-01
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_source
from sqlalchemy import Column, BigInteger, Integer, String, Text, ForeignKey, Table, create_engine
from sqlalchemy.orm import relationship, backref, sessionmaker
import itertools
import operator
import requests
import six
import sqlalchemy.ext.declarative
import sqlalchemy.schema
import unicodecsv

""" - - - - - - - - - - - - - - - - - -  Table Definition Classes - - - - - - - - - - - - - - - - """

Base = sqlalchemy.ext.declarative.declarative_base()

transcription_factor_subunit = Table(
    'transcription_factor_subunit', Base.metadata,
    Column('transcription_factor__id', Integer, ForeignKey('transcription_factor._id'), index=True),
    Column('subunit_id', Integer, ForeignKey('subunit.id'), index=True),
)
# :obj:`Table`: TranscriptionFactor:Subunit many-to-many association table

transcription_factor_family = Table(
    'transcription_factor_family', Base.metadata,
    Column('transcription_factor__id', Integer, ForeignKey('transcription_factor._id'), index=True),
    Column('family_id', Integer, ForeignKey('family.id'), index=True),
)
# :obj:`Table`: TranscriptionFactor:Family many-to-many association table

transcription_factor_class = Table(
    'transcription_factor_class', Base.metadata,
    Column('transcription_factor__id', Integer, ForeignKey('transcription_factor._id'), index=True),
    Column('class_id', Integer, ForeignKey('class.id'), index=True),
)
# :obj:`Table`: TranscriptionFactor:Class many-to-many association table

matrix_resource = Table(
    'matrix_resource', Base.metadata,
    Column('matrix_id', Integer, ForeignKey('matrix.id'), index=True),
    Column('resource_id', Integer, ForeignKey('resource.id'), index=True),
)
# :obj:`Table`: Matrix:Resource many-to-many association table

transcription_factor_species = Table(
    'transcription_factor_species', Base.metadata,
    Column('transcription_factor__id', Integer, ForeignKey('transcription_factor._id'), index=True),
    Column('species_id', Integer, ForeignKey('species.id'), index=True),
)
# :obj:`Table`: TranscriptionFactor:Class many-to-many association table


class TranscriptionFactor(Base):
    """ Represents a transcription factor

    Attributes:
        subunits (:obj:`list` of :obj:`Subunit`): subunits
        species (:obj:`Species`): species
        classes (:obj:`list` of :obj:`Class`): classes
        families (:obj:`list` of :obj:`Family`): families
        collection (:obj:`Collection`): JASPAR collection that the transcription factor belongs to. Each
            collection is a group of transcription factor binding profiles that were published in a single
            paper.
    """
    __tablename__ = 'transcription_factor'

    _id = Column(Integer, primary_key=True)
    id = Column(String(255), unique=True)
    name = Column(String(255))

    subunits = relationship('Subunit', secondary=transcription_factor_subunit, backref='transcription_factors')

    species = relationship('Species', secondary=transcription_factor_species, backref='transcription_factors')

    classes = relationship('Class', secondary=transcription_factor_class, backref='transcription_factors')
    families = relationship('Family', secondary=transcription_factor_family, backref='transcription_factors')

    collection_id = Column(Integer, ForeignKey('collection.id'), index=True)
    collection = relationship('Collection', backref='transcription_factors')


class Matrix(Base):
    """ Represents a binding profile matrix

    Attributes:
        id (:obj:`int`): unique integer ID of the matrix assigned by JASPAR
        transcription_factor (:obj:`TranscriptionFactor`): transcription factor
        version (:obj:`int`): JASPAR version of the matrix
        positions (:obj:`list` of :obj:`MatrixPosition`): the A, C, G, and T frequencies at each position (1..L) within the matrix
        type_id (:obj:`Type`): method used to determine the binding profile matrix (e.g. SELEX, metamodel, phylogenetic)
        references (:obj:`list` of :obj:`Resource`): references to papers which published the binding profile matrix
    """
    __tablename__ = 'matrix'

    id = Column(Integer, primary_key=True)
    transcription_factor_id = Column(Integer, ForeignKey('transcription_factor.id'), index=True)
    transcription_factor = relationship('TranscriptionFactor', backref='matrices')
    version = Column(Integer)
    type_id = Column(Integer, ForeignKey('type.id'), index=True)
    type = relationship('Type', backref='matrices')
    references = relationship('Resource', secondary=matrix_resource, backref='matrices')

    sqlalchemy.schema.UniqueConstraint('transcription_factor_id', 'version')


class MatrixPosition(Base):
    """ Represents the A, C, G, and T frequencies at a position with a binding profile matrix

    Attributes:
        position (:obj:`int`): position within the matrix
        frequency_a (:obj:`int`): frequency of A at the position
        frequency_c (:obj:`int`): frequency of C at the position
        frequency_g (:obj:`int`): frequency of G at the position
        frequency_t (:obj:`int`): frequency of T at the position
        matrix_id (:obj:`int`): ID of the binding profile matrix that the position belongs to
        matrix (:obj:`Matrix`): binding profile matrix that the position belongs to
    """
    __tablename__ = 'matrix_position'

    id = Column(Integer, primary_key=True)
    position = Column(Integer, index=True)
    frequency_a = Column(Integer)
    frequency_c = Column(Integer)
    frequency_g = Column(Integer)
    frequency_t = Column(Integer)
    matrix_id = Column(Integer, ForeignKey('matrix.id'), index=True)
    matrix = relationship('Matrix', backref='positions')

    sqlalchemy.schema.UniqueConstraint('matrix_id', 'position')


class Collection(Base):
    """ Represents a JASPAR collection of transcription factors. Each collection is a group of
    transcription factor binding profiles that were published in a single paper.

    Attributes:
        id (:obj:`int`): unique numeric id assigned by SQLite
        name (:obj:`str`): name
        transcription_factors (:obj:`list` of :obj:`TranscriptionFactor`): transcription factors that
            belong to the collection
    """
    __tablename__ = 'collection'

    id = Column(Integer, primary_key=True)
    name = Column(String(255), unique=True)


class Type(Base):
    """ Represents a methodology used to construct a binding profile matrix such as 
    SELEX, metamodel, or phylogenetic.

    Attributes:
        id (:obj:`int`): unique numeric id assigned by SQLite
        name (:obj:`str`): name
        matrices (:obj:`list` of :obj:`Matrix`): matrices which were constructed using the type
    """
    __tablename__ = 'type'

    id = Column(Integer, primary_key=True)
    name = Column(String(255), unique=True)


class Family(Base):
    """ Represents a structural sub-class of transcription factors, based on the `TFClass system 
    <http://tfclass.bioinf.med.uni-goettingen.de/tfclass>`_

    Attributes:
        id (:obj:`int`): unique numeric id assigned by SQLite
        name (:obj:`str`): name
        transcription_factors (:obj:`list` of :obj:`TranscriptionFactor`): transcription factors that
            belong to the family
    """
    __tablename__ = 'family'

    id = Column(Integer, primary_key=True)
    name = Column(String(255), unique=True)


class Class(Base):
    """ Represents a structural class of transcription factors, based on the `TFClass system 
    <http://tfclass.bioinf.med.uni-goettingen.de/tfclass>`_

    Attributes:
        id (:obj:`int`): unique numeric id assigned by SQLite
        name (:obj:`str`): name
        transcription_factors (:obj:`list` of :obj:`TranscriptionFactor`): transcription factors that
            belong to the class
    """
    __tablename__ = 'class'

    id = Column(Integer, primary_key=True)
    name = Column(String(255), unique=True)


class Species(Base):
    """ Represents the species of a transcription factor

    Attributes:
        id (:obj:`int`): `NCBI taxonomy <https://www.ncbi.nlm.nih.gov/taxonomy>`_ id
        transcription_factors (:obj:`list` of :obj:`TranscriptionFactor`): transcription factors that belong
            to the species
    """
    __tablename__ = 'species'

    id = Column(Integer, primary_key=True)
    ncbi_id = Column(Integer, unique=True, index=True)


class Resource(Base):
    """ Represents a reference in which a binding profile was published

    Attributes:
        id (:obj:`int`): PubMed/Medline ID
        matrices (:obj:`list` of :obj:`Matrix`): matrices which were published in the reference
    """
    __tablename__ = 'resource'

    id = Column(BigInteger, primary_key=True)


class Subunit(Base):
    """ Represents a subunit of a transcription factor

    Attributes:
        id (:obj:`int`): unique numeric id assigned by SQLite
        uniprot_id (:obj:`str`): UniProt id
        transcription_factors (:obj:`list` of :obj:`TranscriptionFactor`): transcription factors that
            the subunit belongs to
    """
    __tablename__ = 'subunit'

    id = Column(Integer, primary_key=True)
    uniprot_id = Column(String(255), unique=True, index=True)


""" - - - - - - - Iterative Sorting and Type Casting Functions Used for Parsing - - - - - - - - - """


def type_cast_matrix_ids_to_ints(table):
    """ Type cast the first entry in each row of a table to an :obj:`int`

    Args:
        table (:obj:`list` of :obj:`list` of :obj:`str`): A list of lists whose first elements are strings

    Returns:
        :obj:`list` of :obj:`list`: A list of lists whose first elements are ints
    """
    for row in table:
        row[0] = int(row[0])
    return table


def type_cast_matrix_positions_and_frequencies(table):
    """ Type cast the matrix position and frequencies of matrix data into ints and floats

    Args:
        table (:obj:`list` of :obj:`list` of :obj:`str`): Table of matrix data (matrix id, base, position, frequency)

    Returns:
        :obj:`list` of :obj:`list`: table of matrix data with positions and frequencies converted to ints and floats
    """
    for row in table:
        row[2] = int(row[2])
        row[3] = float(row[3])
    return table


def group_by_matrix_ids(table):
    """ Group a table by matrix ids

    Args:
        table (:obj:`list` of :obj:`list`): List of JASPAR data, each of which represents a row in a JASPAR data file
            as a list with the first element being the matrix id of the data

    Returns:
        :obj:`dict`: JASPAR data, grouped by by matrix ids
    """
    sorted_table = sorted(table, key=operator.itemgetter(0))
    dictionary = {}
    for id, rows in itertools.groupby(sorted_table, key=operator.itemgetter(0)):
        dictionary[id] = [row[1:] for row in rows]
    return dictionary


def group_by_position(table):
    """ Group a table of binding profiles by position

    Args:
        table (:obj:`list` of :obj:`list` of binding profile positions): table of matrix binding positions

    Returns:
        :obj:`dict`: binding profiles, grouped by position
    """
    sorted_table = sorted(table, key=operator.itemgetter(1))
    dictionary = {}
    for position, rows in itertools.groupby(sorted_table, key=operator.itemgetter(1)):
        dictionary[position] = [[row[0], row[2]] for row in rows]
    return dictionary


""" - - - -  Collecting and Parsing Data From Website then Adding to Session DB - -  - - """


class Jaspar(data_source.HttpDataSource):
    """ A local SQLite copy of the `JASPAR <http://jaspar.genereg.net>`_ database of transcription factor binding profiles """
    base_model = Base
    ENDPOINT_DOMAINS = {
        'jaspar': 'http://jaspar.genereg.net/html/DOWNLOAD/database/',
    }

    def load_content(self):
        if self.verbose:
            print('Downloading and parsing JASPAR ...')
        self.download_parse_and_load_data()

        if self.verbose:
            print('Saving JASPAR ...')
        self.session.commit()

        if self.verbose:
            print('Done.')

    def download_parse_and_load_data(self):
        """ Downloads and parses all of the data from the JASPAR website and appends the data to an SQLlite session. """
        database_url = self.ENDPOINT_DOMAINS['jaspar']
        session = self.session
        req = self.requests_session

        # parse subunits (matrix id, UniProt id)
        response = req.get(database_url + 'MATRIX_PROTEIN.txt')
        response.raise_for_status()
        file = six.BytesIO(response.content)
        subunits_table = type_cast_matrix_ids_to_ints(list(unicodecsv.reader(file, delimiter='\t')))
        subunits = group_by_matrix_ids(subunits_table)
        if self.verbose:
            print('  Downloaded and parsed subunit data.')

        # parse species (matrix id, NCBI taxonomy ID)
        response = req.get(database_url + 'MATRIX_SPECIES.txt')
        response.raise_for_status()
        file = six.BytesIO(response.content)
        species_table = type_cast_matrix_ids_to_ints(list(unicodecsv.reader(file, delimiter='\t')))
        species = group_by_matrix_ids(species_table)
        if self.verbose:
            print('  Downloaded and parsed species data.')

        # parse annotations (classes, families, types, references) (matrix id, annotation type, annotated value)
        response = req.get(database_url + 'MATRIX_ANNOTATION.txt')
        response.raise_for_status()
        file = six.BytesIO(response.content)
        annotations_table = type_cast_matrix_ids_to_ints(list(unicodecsv.reader(file, delimiter='\t')))
        classes_table = []
        families_table = []
        references_table = []
        types_table = []
        for row in annotations_table:
            if row[1] == 'class':
                classes_table.append([row[0], row[2]])
            elif row[1] == 'family':
                families_table.append([row[0], row[2]])
            elif row[1] == 'medline':
                references_table.append([row[0], row[2]])
            elif row[1] == 'type':
                types_table.append([row[0], row[2]])
        classes = group_by_matrix_ids(classes_table)
        families = group_by_matrix_ids(families_table)
        references = group_by_matrix_ids(references_table)
        types = group_by_matrix_ids(types_table)
        if self.verbose:
            print('  Downloaded and parsed annotation data.')

        # parse binding profiles (matrix id, position, base, frequency)
        response = req.get(database_url + 'MATRIX_DATA.txt')
        response.raise_for_status()
        file = six.BytesIO(response.content)
        matrix_data_table = type_cast_matrix_positions_and_frequencies(
            type_cast_matrix_ids_to_ints(list(unicodecsv.reader(file, delimiter='\t'))))
        matrix_data = group_by_matrix_ids(matrix_data_table)
        for id, position_frequencies in matrix_data.items():
            matrix_data[id] = group_by_position(position_frequencies)
        if self.verbose:
            print('  Downloaded and parsed matrix binding data.')

        # parse transcription factors (matrix id, transcripton factor id, collection, matrix version, transcription factor name)
        response = req.get(database_url + 'MATRIX.txt')
        response.raise_for_status()
        file = six.BytesIO(response.content)
        matrices = list(unicodecsv.reader(file, delimiter='\t'))
        if self.verbose:
            print('  Downloaded and parsed transcription factor data.')

        # add transcription factors to the SQLite database
        if self.verbose:
            print('  Loading transcription factors into the SQLite database ...')
        for id, collection, tf_id, version, tf_name in matrices:
            id = int(id)

            tf = self.get_or_create_object(TranscriptionFactor, id=tf_id)
            if not tf.matrices or int(version) > max([matrix.version for matrix in tf.matrices]):
                tf.name = tf_name
            tf.collection = self.get_or_create_object(Collection, name=collection)

            if id in subunits:
                for subunit in subunits[id]:
                    assert isinstance(subunit[0], six.string_types)
                    tf.subunits.append(self.get_or_create_object(Subunit, uniprot_id=subunit[0]))

            if id in species:
                for specie in species[id]:
                    if specie[0] != '-':
                        tf.species.append(self.get_or_create_object(Species, ncbi_id=int(specie[0])))

            if id in classes:
                for cls in classes[id]:
                    assert isinstance(cls[0], six.string_types)
                    tf.classes.append(self.get_or_create_object(Class, name=cls[0]))

            if id in families:
                for family in families[id]:
                    assert isinstance(family[0], six.string_types)
                    tf.families.append(self.get_or_create_object(Family, name=family[0]))

            matrix = Matrix(id=id, transcription_factor=tf, version=int(version))

            if id in types:
                for type in types[id]:
                    assert isinstance(type[0], six.string_types)
                    matrix.type = self.get_or_create_object(Type, name=type[0])

            if id in references:
                for ref in references[id]:
                    if ref[0] not in ['-', 'unpublished']:
                        matrix.references.append(self.get_or_create_object(Resource, id=int(ref[0])))

            for position, freqs in matrix_data[id].items():
                matrix_position = MatrixPosition(matrix=matrix, position=position)

                for freq in freqs:
                    if freq[0] == 'A':
                        matrix_position.frequency_a = freq[1]
                    elif freq[0] == 'C':
                        matrix_position.frequency_c = freq[1]
                    elif freq[0] == 'G':
                        matrix_position.frequency_g = freq[1]
                    elif freq[0] == 'T':
                        matrix_position.frequency_t = freq[1]
                    else:
                        raise Exception('Invalid base {}'.format(freq[0]))

        if self.verbose:
            print('  Done.')
