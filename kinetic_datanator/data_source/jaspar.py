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
        type_info (:obj:`Type`): type                
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


# TODO: Figure out if there are any other attributes within the family (Ex. global IDS)
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
    """ Sorts data by matrix_id

    Args:
        data (:obj:`list`): A list with matrix_id as first column

    Returns:
        :obj:`list`: sorted list by matrix_id
    """
    sorted_list = sorted(data, key=itemgetter(0), reverse=False)
    return sorted_list


def group_by_jaspar_id(data):
    """ Groups list by matrix_id

    Args:
        data (:obj:`list`): List of sorted matrix binding data

    Returns:
        :obj:`list`: grouped by matrix_id

    """
    return {k: g in groupby(data, lambda x: x[0])}

def group_by_position(data):
    """ Groups data post-matrix_id group by position

    Args:
        data (:obj:`list`): List of sorted matrix binding data

    Returns:
        :obj:`list`: grouped by position

    """
    group = []
    for k, g in groupby(data, lambda x: x[2]):
        group.append(list(g))
    return group

""" - - - - - - - - - - - - - - - - - -  Retrieval Functions - - - - - - - - - - - - - - - - - - """


def call_attribute(key_list, data_list, jid):
    """ Calls relevant attribute from matrix_id. Returns data point and keylist

    Args:
        key_list (:obj:`list`): List of matrix_id keys for attribute
        data_list (:obj:`list`): Grouped list by matrix_id of attribute (Ex. Uniprot ID)
        jid (:obj:`int`): matrix_id used to locate corresponding attribute

    Returns:
        :obj:`str`: designated attribute for a given matrix_id
        :obj:`list`: list of keys remaining after calling for an attribute
        :obj:`list`: list of data remaining after calling for an attribute


    """
    answer = []
    ind = 0
    while True:
        if key_list[ind] == jid:
            for i in range(0, len(data_list[ind])):
                answer.append(data_list[ind][i][1])
            key_list.pop(ind)
            data_list.pop(ind)
            return answer, key_list, data_list
            break
        elif key_list[ind] < jid:
            ind += 1
            continue
        else:
            return [], key_list, data_list
            break


def call_data(key_list, data_list, jid):
    """ Calls data by matrix_id

    Args:
        key_list (:obj:`list`): List of matrix_id keys for data
        data_list (:obj:`list`): Grouped list by matrix_id of transcription binding matrix frequencies and positions
        jid (:obj:`int`): matrix_id used to locate binding matrix data

    Returns:
        :obj:`list`: frequency of A given a matrix_id
        :obj:`list`: frequency of C given a matrix_id
        :obj:`list`: frequency of G given a matrix_id
        :obj:`list`: frequency of T given a matrix_id
        :obj:`list`: list of keys remaining after calling for an attribute
        :obj:`list`: list of data remaining after calling for an attribute

    """
    A = []
    C = []
    G = []
    T = []
    ind = 0
    while True:
        if key_list[ind] == jid and len(data_list[ind][0]) == 4:
            for position in range(0, len(data_list[ind])):
                A.append(data_list[ind][position][0][3])
                C.append(data_list[ind][position][1][3])
                G.append(data_list[ind][position][2][3])
                T.append(data_list[ind][position][3][3])
            key_list.pop(ind)
            data_list.pop(ind)
            return A, C, G, T, key_list, data_list
            break
        elif len(data_list[ind][0]) != 4:
            key_list.pop(ind)
            data_list.pop(ind)
            continue
        elif key_list[ind] != jid:
            ind += 1
            continue
        else:
            return [], [], [], [], key_list, data_list
            break

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
        subunits = group_by_jaspar_id(sorted_list)
        if self.verbose:
            print('Collected subunit Data.')

        response = req.get(database_url + 'MATRIX_SPECIES.txt')
        response.raise_for_status()
        f = BytesIO(response.content)
        species = reader(f, delimiter='\t')
        species = make_first_col_ints(species)
        sorted_list = sort(species)
        species = group_by_jaspar_id(sorted_list)
        if self.verbose:
            print('Collected Species Data.')

        response = req.get(database_url + 'MATRIX_ANNOTATION.txt')
        response.raise_for_status()
        f = BytesIO(response.content)
        annots = reader(f, delimiter='\t')
        annots = make_first_col_ints(annots)
        sorted_list = sort(annots)
        class_info = []
        fam_info = []
        med_info = []
        type_info = []
        for row in sorted_list:
            if row[1] == 'class':
                class_info.append([row[0], row[2]])
            elif row[1] == 'family':
                fam_info.append([row[0], row[2]])
            elif row[1] == 'medline':
                med_info.append([row[0], row[2]])
            elif row[1] == 'type':
                type_info.append([row[0], row[2]])
        class_info = group_by_jaspar_id(class_info)
        fam_info = group_by_jaspar_id(fam_info)
        med_info = group_by_jaspar_id(med_info)
        type_info = group_by_jaspar_id(type_info)
        if self.verbose:
            print('Collected Annotation Data.')

        response = req.get(database_url + 'MATRIX_DATA.txt')
        response.raise_for_status()
        f = BytesIO(response.content)
        matrix_data = reader(f, delimiter='\t')
        matrix_data = make_data_int(matrix_data)
        sorted_list = sort(matrix_data)
        matrix_data = group_by_jaspar_id(sorted_list)
        for row in matrix_data:
            row = sorted(row, key=itemgetter(2), reverse=False)
            row = group_by_position(row)
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

            if matrix_id in subunit:
                for uniprot_id in subunits[matrix_id]:
                    assert isinstance(uniprot_id, six.string_types)
                    tf.subunits.append(Subunit(uniprot_id=uniprot_id))

            ncbi_id, key_species, species = call_attribute(key_species, species, matrix_id)        

            for n in range(0, len(ncbi_id)):
                assert isinstance(ncbi_id[n], six.string_types)
                species = Species(ncbi_id=ncbi_id[n], matrix_id=matrix_id)

            class_name, key_class, class_info = call_attribute(key_class, class_info, matrix_id)
            for n in range(0, len(class_name)):
                assert isinstance(class_name[n], six.string_types)
                protein_class = Class(class_name=class_name[n], matrix_id=matrix_id)

            fam_name, key_fam, fam_info = call_attribute(key_fam, fam_info, matrix_id)
            for n in range(0, len(fam_name)):
                assert isinstance(fam_name[n], six.string_types)
                family = Family(family_name=fam_name[n], matrix_id=matrix_id)

            medline_id, key_med, med_info = call_attribute(key_med, med_info, matrix_id)
            for n in range(0, len(medline_id)):
                assert isinstance(medline_id[n], six.string_types)
                resources = Resource(medline_id=medline_id[n], matrix_id=matrix_id)

            type_name, key_type, type_info = call_attribute(key_type, type_info, matrix_id)
            for n in range(0, len(type_name)):
                assert isinstance(type_name[n], six.string_types)
                type = Type(type_name=type_name[n], matrix_id=matrix_id)

            matrix = Matrix(tf=tf, version=version,
                            collection=collection,
                            references=[resources],
                            type_info=[type], subunits=[subunit],
                            family_info=[family], class_info=[protein_class], species=[species])

            frequency_a, frequency_c, frequency_g, frequency_t, key_data, matrix_data = call_data(key_data, matrix_data, matrix_id)
            for position in range(0, len(frequency_a)):
                a = frequency_a[position]
                assert isinstance(a, float)
                c = frequency_c[position]
                assert isinstance(c, float)
                g = frequency_g[position]
                assert isinstance(g, float)
                t = frequency_t[position]
                assert isinstance(t, float)
                position = position + 1
                assert isinstance(position, int)
                matrix_position = MatrixPosition(matrix=matrix, position=position, frequency_a=a,
                                                 frequency_c=c, frequency_g=g, frequency_t=t)
        print('Done.')
