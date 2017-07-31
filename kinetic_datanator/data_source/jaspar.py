"""
This code downloads the JASPAR database of transcription factor binding motifs (http://jaspar.genereg.net/)
via a seris of text files, parses them, and stores them in an SQLlite database.

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-07-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

from sqlalchemy import Column, BigInteger, Integer, String, Text, ForeignKey, Table, create_engine
from sqlalchemy.orm import relationship, backref, sessionmaker
from kinetic_datanator.core import data_source
from unicodecsv import reader
from itertools import groupby
from operator import itemgetter
from six import BytesIO
import sqlalchemy.ext.declarative
import requests
import six

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

# transcription_factor_species = Table(
#     'transcription_factor_species', Base.metadata,
#     Column('transcription_factor__id', Integer, ForeignKey('transcription_factor._id'), index=True),
#     Column('species_id', Integer, ForeignKey('species.id'), index=True),
# )
# :obj:`Table`: TranscriptionFactor:Class many-to-many association table

class TranscriptionFactor(Base):
    """
    Attributes:
        subunits (:obj:`list` of :obj:`Subunit`): subunits
        species (:obj:`Species`): species
        classes (:obj:`list` of :obj:`Class`): classes
        families (:obj:`list` of :obj:`Family`): families
        collection (:obj:`Collection`): collection from which matrix stems from
    """

    __tablename__ = 'transcription_factor'

    _id = Column(Integer, primary_key=True)
    id = Column(String(255))
    name = Column(String(255))

    subunits = relationship('Subunit', secondary=transcription_factor_subunit, backref='transcription_factors')

    species_id = Column(Integer, ForeignKey('species.id'))
    species = relationship('Species', backref='transcription_factors')

    classes = relationship('Class', secondary=transcription_factor_class, backref='transcription_factors')
    families = relationship('Family', secondary=transcription_factor_family, backref='transcription_factors')

    collection_id = Column(Integer, ForeignKey('collection.id'))
    collection = relationship('Collection', backref='transcription_factors')


class Matrix(Base):
    """ Represents the identification of the binding profile matrix and transcription factor found in Jaspar

    Attributes:
        id (:obj:`int`): ID of Matrix (INTEGER)
        transcription_factor (:obj:`TranscriptionFactor`): transcription factor
        version (:obj:`int`): Jaspar version of given matrix
        positions (:obj:`list` of :obj:`MatrixPosition`):
        type_id (:obj:`Type`): type
        references (:obj:`list` of :obj:`Resource`): references
    """
    __tablename__ = 'matrix'

    id = Column(Integer, primary_key=True)
    transcription_factor_id = Column(Integer, ForeignKey('transcription_factor.id'))
    transcription_factor = relationship('TranscriptionFactor', backref='matrices')
    version = Column(Integer)
    type_id = Column(Integer, ForeignKey('type.id'))
    type = relationship('Type', backref='matrices')
    references = relationship('Resource', secondary=matrix_resource, backref='matrices')


class MatrixPosition(Base):
    """ Represents the Matrix based transcription factor binding profile

    Attributes:
        position (:obj:`int`): integer binding position of corresponding frequency
        frequency_a (:obj:`int`): frequency of A
        frequency_c (:obj:`int`): frequency of C
        frequency_g (:obj:`int`): frequency of G
        frequency_t (:obj:`int`): frequency of T
        matrix_id (:obj:`int`): ID of Matrix
        matrix (:obj:`Matrix`)
    """
    __tablename__ = 'matrix_position'
    id = Column(Integer, primary_key=True)
    position = Column(Integer)
    frequency_a = Column(Integer)
    frequency_c = Column(Integer)
    frequency_g = Column(Integer)
    frequency_t = Column(Integer)

    matrix_id = Column(Integer, ForeignKey('matrix.id'))
    matrix = relationship('Matrix', backref='positions')


class Collection(Base):
    """ Represents a collection

    Attributes:
        id (:obj:`int`): assigned numeric id
        name (:obj:`str`): name
        transcription_factors (:obj:`list` of :obj:`TranscriptionFactor`): transcription factors
    """

    __tablename__ = 'collection'
    id = Column(Integer, primary_key=True)
    name = Column(String(255))


class Type(Base):
    """ Represents Methodology for Matrix Construction

    Attributes:
        id (:obj:`int`): id
        name (:obj:`str`): name
        matrices (:obj:`list` of :obj:`Matrix`): matrices
    """
    __tablename__ = 'type'

    id = Column(Integer, primary_key=True)
    name = Column(String(255))


class Family(Base):
    """ Represents protein family

    Attributes:
        id (:obj:`int`): id
        name (:obj:`str`): name
        transcription_factors (:obj:`list` of :obj:`TranscriptionFactor`): transcription factors
    """
    __tablename__ = 'family'

    id = Column(Integer, primary_key=True)
    name = Column(String(255))


class Class(Base):
    """ Represents protien class

    Attributes:
        id (:obj:`int`): id
        name (:obj:`str`): name
        transcription_factors (:obj:`list` of :obj:`TranscriptionFactor`): transcription factors
    """

    __tablename__ = 'class'

    id = Column(Integer, primary_key=True)
    name = Column(String(255))


class Species(Base):
    """ Represents Species of Observation

    Attributes:
        id (:obj:`int`): NCBI id
        transcription_factors (:obj:`list` of :obj:`TranscriptionFactor`): transcription factors
    """

    __tablename__ = 'species'

    id = Column(Integer, primary_key=True)


class Resource(Base):
    """ Represents external resource

    Attributes:
        id (:obj:`int`): Pubmed/Medline ID
        matrices (:obj:`list` of :obj:`Matrix`): matrices
    """

    __tablename__ = 'resource'

    id = Column(BigInteger, primary_key=True)


class Subunit(Base):
    """ Represents a Transcription Factor from Jaspar Database

    Attributes:
        id (:obj:`int`): id
        uniprot_id (:obj:`str`): uniprot id
        transcription_factors (:obj:`list` of :obj:`TranscriptionFactor`): transcription factors
    """
    __tablename__ = 'subunit'

    id = Column(Integer, primary_key=True)
    uniprot_id = Column(String(255))

""" - - - - - - - Iterative Sorting and Data Type Functions Used for Parsing - - - - - - - - - """


def make_first_col_ints(data):
    """ Converts first entry in each row of a table to ints

    Args:
        data (:obj:`list`): A list with string as first element

    Returns:
        :obj:`list`: list of with ints as first elements

    """

    data = list(data)
    for i in range(0, len(data)):
        data[i][0] = int(data[i][0])
    return data


def make_data_int(data):
    """ Converts matrix_id and position into ints and frequency into floats

    Args:
        data (:obj:`list`): Matrix binding data

    Returns:
        :obj:`list`: list of matrix_ids and position converted to int and frequency converted to float

    """
    data = list(data)
    for i in range(0, len(data)):
        data[i][0] = int(data[i][0])
        data[i][2] = int(data[i][2])
        data[i][3] = float(data[i][3])
    return data


def sort(data):
    """ Sorts data by matrix_id and/or position

    Args:
        data (:obj:`list`): A list with matrix_id as first column

    Returns:
        :obj:`list`: sorted list by matrix_id
    """
    sorted_list = sorted(data, key=itemgetter(0), reverse=False)
    return sorted_list

def group(data, round):
    """ Groups list by matrix_id

    Args:
        data (:obj:`list`): List of sorted matrix binding data

    Returns:
        :obj:`dict`: grouped by matrix_id

    """

    if round == 1:
        groups = groupby(data, key = itemgetter(0))
        return {k:[x[1:len(x)] for x in v] for k, v  in groups}
    elif round == 2:
        groups = groupby(data, key = itemgetter(1))
        return {k:[[x[0],x[2]] for x in v] for k, v  in groups}
    else:
        return TypeError('No Round Specified')


""" - - - -  Collecting and Parsing Data From Website then Adding to Session DB - -  - - """


class Jaspar(data_source.HttpDataSource):
    """ A local sqlite copy of the ECMDB database

    """
    base_model = Base
    ENDPOINT_DOMAINS = {
        'jaspar': 'http://jaspar.genereg.net/html/DOWNLOAD/database/',
    }

    def load_content(self):
        if self.verbose:
            print('Downloading and parsing Jaspar ...')
        self.parse_db()

        if self.verbose:
            print('Committing Jaspar ...')
        self.session.commit()

        if self.verbose:
            print('Done.')

    def parse_db(self):
        """ Collects and Parses all data from jaspar DB website and adds to SQLlite DB

        Args:
            session (:obj:`session object`): SQLAlchemy session object
            req (:obj:`requests object`): Requests session object
        """
        database_url = self.ENDPOINT_DOMAINS['jaspar']
        session = self.session
        req = self.requests_session

        response = req.get(database_url + 'MATRIX_PROTEIN.txt')
        response.raise_for_status()
        f = BytesIO(response.content)
        subunits = reader(f, delimiter='\t')
        subunits = make_first_col_ints(subunits)
        sorted_list = sort(subunits)
        subunits = group(sorted_list,1)
        if self.verbose:
            print('Collected Subunit Data.')

        response = req.get(database_url + 'MATRIX_SPECIES.txt')
        response.raise_for_status()
        f = BytesIO(response.content)
        species = reader(f, delimiter='\t')
        species = make_first_col_ints(species)
        sorted_list = sort(species)
        species = group(sorted_list,1)
        if self.verbose:
            print('Collected Species Data.')

        response = req.get(database_url + 'MATRIX_ANNOTATION.txt')
        response.raise_for_status()
        f = BytesIO(response.content)
        annots = reader(f, delimiter='\t')
        annots = make_first_col_ints(annots)
        sorted_list = sort(annots)
        cls = []
        family = []
        pubmed_id = []
        type_id = []
        for row in sorted_list:
            if row[1] == 'class':
                cls.append([row[0], row[2]])
            elif row[1] == 'family':
                family.append([row[0], row[2]])
            elif row[1] == 'medline':
                pubmed_id.append([row[0], row[2]])
            elif row[1] == 'type':
                type_id.append([row[0], row[2]])
        cls = group(cls,1)
        family = group(family,1)
        pubmed_id = group(pubmed_id,1)
        type_id = group(type_id,1)
        if self.verbose:
            print('Collected Annotation Data.')

        response = req.get(database_url + 'MATRIX_DATA.txt')
        response.raise_for_status()
        f = BytesIO(response.content)
        matrix_data = reader(f, delimiter='\t')
        matrix_data = make_data_int(matrix_data)
        sorted_list = sort(matrix_data)
        matrix_data = group(sorted_list,1)
        for k in matrix_data:
            matrix_data[k] = sorted(matrix_data[k], key = itemgetter(1))
            matrix_data[k] = group(matrix_data[k],2)
        if self.verbose:
            print('Collected Matrix Binding Data.')


        response = req.get(database_url + 'MATRIX.txt')
        response.raise_for_status()
        f = BytesIO(response.content)
        units = list(reader(f, delimiter='\t'))
        if self.verbose:
            print('Table Creation Started...')
        for unit in units:
            # error checking
            matrix_id = int(unit[0])
            assert isinstance(matrix_id, int)

            tf_id = unit[2]
            assert isinstance(tf_id, six.string_types)

            version = int(unit[3])
            assert isinstance(version, int)

            tf_name = unit[4]
            assert isinstance(tf_name, six.string_types)

            collection_name = unit[1]
            assert isinstance(collection_name, six.string_types)

            # append transcription factors and matrices to the database
            tf = self.get_or_create_object(TranscriptionFactor, id=tf_id, name=tf_name)
            tf.collection = self.get_or_create_object(Collection, name=collection_name)

            if matrix_id in subunits:
                for uniprot_id in subunits[matrix_id]:
                    assert isinstance(uniprot_id[0], six.string_types)
                    tf.subunits.append(self.get_or_create_object(Subunit,uniprot_id = uniprot_id[0]))
                subunits.pop(matrix_id)

            if matrix_id in species:
                for ncbi in species[matrix_id]:
                    if ncbi[0] != '-':
                        ncbi[0] = int(ncbi[0])
                        assert isinstance(ncbi[0], six.integer_types)
                        tf.species = self.get_or_create_object(Species, id = ncbi[0])
                species.pop(matrix_id)

            if matrix_id in cls:
                for name in cls[matrix_id]:
                    assert isinstance(name[0], six.string_types)
                    tf.classes.append(self.get_or_create_object(Class, name = name[0]))
                cls.pop(matrix_id)

            if matrix_id in family:
                for name in family[matrix_id]:
                    assert isinstance(name[0], six.string_types)
                    tf.families.append(self.get_or_create_object(Family, name = name[0]))
                family.pop(matrix_id)


            matrix = self.get_or_create_object(Matrix, id = matrix_id, transcription_factor = tf,
                                                version = version)

            if matrix_id in type_id:
                for name in type_id[matrix_id]:
                    assert isinstance(name[0], six.string_types)
                    matrix.type = self.get_or_create_object(Type, name = name[0])
                type_id.pop(matrix_id)

            if matrix_id in pubmed_id:
                for id_ in pubmed_id[matrix_id]:
                    if id_[0] != '-' and id_[0] != 'unpublished':
                        id_[0] =  int(id_[0])
                        assert isinstance(id_[0], six.integer_types)
                        matrix.references.append(self.get_or_create_object(Resource, id = id_[0]))
                    else:
                        matrix.references.append(self.get_or_create_object(Resource, id = 0))
                pubmed_id.pop(matrix_id)

            if matrix_id in matrix_data:
                for position in matrix_data[matrix_id]:
                    freq = matrix_data[matrix_id][position]
                    matrix_position = self.get_or_create_object(MatrixPosition, position = position,
                        frequency_a = freq[0][1], frequency_c = freq[1][1], frequency_g = freq[2][1],
                        frequency_t = freq[3][1], matrix = matrix)
                matrix_data.pop(matrix_id)

        print('Done.')
