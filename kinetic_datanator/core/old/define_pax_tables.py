"""
This codebase takes the txt files of the PaxDB protein abundance database
and inserts them into an SQL database

define_tables.py - defines the python classes corresponding to the tables in
the resulting SQL database

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2017 June 3
:Copyright: 2017, Karr Lab
:License: MIT
"""

from sqlalchemy import create_engine, ForeignKey
from sqlalchemy import Column, Integer, String, Numeric
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref

engine = create_engine('sqlite:///pax.db')
Base = declarative_base()

""" -------------------------------------------------------------------------"""
class Taxon(Base):
    """  Represents a species
    Attributes:
        ncbi_id (:obj:`int`): NCBI id
        species_name (:obj:`str`): name and possibly genetic variant
    """

    __tablename__ = 'taxon'
    ncbi_id      = Column(Integer,primary_key=True)
    species_name = Column(String(255))

""" -------------------------------------------------------------------------"""
class Dataset(Base):
    """  Represents a given dataset (typically results form a single paper)
    Attributes:
        ncbi_id     (:obj:`int`): NCBI id - linked to the 'taxon' table
        publication (:obj:`str`): URL of the corresponding publication
        file_name   (:obj:`str`): the name of text file corresponding to the dataset
        score       (:obj:`flt`): PaxDb's internal quality score
        weight      (:obj:`int`): TBA
        coverage    (:obj:`int`): what percentage of the genome is coevred by the datatset
    """

    __tablename__  = 'dataset'
    id          = Column(Integer, primary_key=True)
    publication = Column(String(255))
    file_name   = Column(String(255))
    score       = Column(Numeric)
    weight      = Column(Numeric)
    coverage    = Column(Numeric)

    taxon_ncbi_id = Column(Integer, ForeignKey('taxon.ncbi_id'))
    taxon         = relationship('Taxon', backref=backref('datasets'), foreign_keys=[taxon_ncbi_id])

""" -------------------------------------------------------------------------"""
class Protein(Base):
    """  Represents a protein
    Attributes:
        protein_id (:obj:`int`): PaxDB's internal numerical protein ID
        string_id  (:obj:`str`): Ensembl ID of protein
    """

    __tablename__ = 'protein'
    protein_id = Column(Integer, primary_key=True)
    string_id  = Column(String(255))

""" -------------------------------------------------------------------------"""
class Observation(Base):
    """  Represents a protein
    Attributes:
        protein_id (:obj:`int`): PaxDB's internal numerical protein ID
        dataset_id (:obj:`int`): ID of the database - linked to the 'dataset' table
        abundance  (:obj:`flt`): Normalized abudnance of the protein
    """

    __tablename__ = 'observation'
    id         = Column(Integer, primary_key=True)
    abundance  = Column(Numeric)

    dataset_id = Column(Integer, ForeignKey('dataset.id'))
    protein_id = Column(Integer, ForeignKey('protein.protein_id'))
    dataset    = relationship('Dataset',backref=backref('observation'),foreign_keys=[dataset_id])
    protein    = relationship('Protein',backref=backref('observation'),foreign_keys=[protein_id])

    #FORMAT: foreign_table = relationship('foreign_class',backref=backref('self_table'),foreign_keys=[self_column])

""" -------------------------------------------------------------------------"""

Base.metadata.drop_all(engine)
Base.metadata.create_all(engine)
