# -*- coding: utf-8 -*-
"""
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017 July 10th
:Copyright: 2017, Karr Lab
:License: MIT
"""


##TODO: Figure out dependencies for this module -- kinetic_datanator -> wc_util -> libgit2?
from sqlalchemy import Column, Integer, String, Numeric, ForeignKey, UnicodeText, create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref

engine = create_engine('sqlite:///jaspar_prac.db')
Base = declarative_base()
# :obj:`Base`: base model for local sqlite database



class MatrixObservation(Base):
    """ Represents the identification of the binding profile matrix and transcription factor found in Jaspar
        jaspar_matrix_ID (:obj:`int`): ID of Matrix (INTEGER)
        jaspar_tf_ID (:obj:`str`): ID of Transcription Factor (STRING)
        version (:obj:`int`): Jaspar version of given matrix
        tf_name (:obj:`str`): transcription factor name (standardized Entrez gene symbols)
        jaspar_collection (:obj:`str`): collection from which matrix stems from
    """
    __tablename__ = 'matrixobservation'

    jaspar_matrix_ID = Column(Integer, primary_key = True)
    jaspar_tf_ID = Column(String(255))
    version = Column(Integer)
    tf_name = Column(String(255))
    jaspar_collection = Column(String(255))

    data = relationship('BindingMatrix', backref = 'observation')
    references = relationship('Resources', backref = 'observation')
    type_info = relationship('Type', backref = 'observation')
    family_info = relationship('Family', backref = 'observation')
    class_info = relationship('Class', backref = 'observation')
    transcriptionfactor = relationship('TranscriptionFactor', backref = 'observation')
    # resource_relation = relationship('Resources', backref = backref('matrix_relation'))
    # type_relation = relationship('Type', backref = backref('matrix_relation'))
    # tf_relation = relationship('TranscriptionFactor', backref = backref('matrix_relation'))
    # family_relation = relationship('Family', backref = backref('matrix_relation'))
    # class_relation = relationship('Class', backref = backref('matrix_relation'))

class BindingMatrix(Base):
    """ Represents the Matrix based transcription factor binding profile

    Attributes:
        position (:obj:`int`): integer binding position of corresponding frequency
        frequency_A (:obj:`float`): frequency of A
        frequency_C (:obj:`float`): frequency of C
        frequency_G (:obj:`float`): frequency of G
        frequency_T (:obj:`float`): frequency of T
        jaspar_matrix_ID (:obj:`int`): ID of Matrix
    """
    __tablename__ = 'bindingmatrix'


    position = Column(Integer, primary_key = True)
    frequency_A = Column(Numeric)
    frequency_C = Column(Numeric)
    frequency_G = Column(Numeric)
    frequency_T = Column(Numeric)

    jaspar_matrix_ID = Column(Integer, ForeignKey('matrixobservation.jaspar_matrix_ID'))
    # matrix_relation = relationship('MatrixObservation', backref= 'binding_relation', foreign_keys = [jaspar_matrix_ID])


class Type(Base):
    """ Represents Methodology for Matrix Construction

    Attributes:
        type_name (:obj:`str`): name
        jaspar_matrix_ID (:obj:`int`): ID of Matrix
    """
    __tablename__ = 'type'

    id = Column(Integer, primary_key = True)
    type_name = Column(String(255))

    jaspar_matrix_ID = Column(Integer, ForeignKey('matrixobservation.jaspar_matrix_ID'))
    # matrix_relation = relationship('MatrixObservation', backref=backref('type'), foreign_keys = [jaspar_matrix_ID])

class Family(Base):
    """ Represents protein family

    Attributes:
        familyname (:obj:`str`): name
        jaspar_matrix_ID (:obj:`int`): ID of Matrix
    """
    __tablename__ = 'family'

    id = Column(Integer, primary_key = True)
    family_name = Column(UnicodeText(64))

    jaspar_matrix_ID = Column(String, ForeignKey('matrixobservation.jaspar_matrix_ID'))
    matrix_relation = relationship('MatrixObservation', backref = backref('family'), foreign_keys= [jaspar_matrix_ID])

## TODO: Figure out if there are any other attributes within the family (Ex. global IDS)

class Class(Base):
    """ Represents protien class

    Attributes:
        class_name (:obj:`str`): name
        jaspar_matrix_ID (:obj:`int`): ID of Matrix
    """

    __tablename__ = 'class'

    id = Column(Integer, primary_key = True)
    class_name = Column(String(255))

    jaspar_matrix_ID = Column(String, ForeignKey('matrixobservation.jaspar_matrix_ID'))
    # matrix_relation = relationship('MatrixObservation', backref = backref('class'), foreign_keys= [jaspar_matrix_ID])


class Resources(Base):
    """ Represents external resource

    Attributes:
        medline_id (:obj:`int`): Pubmed/Medline Reference ID
        jaspar_matrix_ID (:obj:`int`): ID of Matrix
    """

    __tablename__ = 'resource'

    id = Column(Integer, primary_key = True)
    medline_id = Column(String(255))

    jaspar_matrix_ID = Column(Integer, ForeignKey('matrixobservation.jaspar_matrix_ID'))
    # matrix_relation = relationship('MatrixObservation', backref=backref('resource'), foreign_keys= [jaspar_matrix_ID])

class TranscriptionFactor(Base):
    """ Represents a Transcription Factor from Jaspar Database

    Attributes:
        uniprot_id (:obj:`str`): uniprot id
        NCBI_id (:obj:`int`): NCBI species ID from which TF was found

        jaspar_matrix_ID (:obj:`int`): ID of Matrix
        tf_name (:obj:`str`): transcription factor name (standardized Entrez gene symbols)
    """
    __tablename__ = 'transcriptionfactor'

    id = Column(Integer, primary_key = True)
    uniprot_id = Column(String(255))
    NCBI_id = Column(String(255))

    jaspar_matrix_ID = Column(Integer, ForeignKey('matrixobservation.jaspar_matrix_ID'))
    # matrix_relation = relationship('MatrixObservation', backref=backref('transcriptionfactor'),
    #                                 foreign_keys = [jaspar_matrix_ID])

Base.metadata.drop_all(engine)
Base.metadata.create_all(engine)
