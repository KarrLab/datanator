"""
This codebase takes CORUM protein complexes database and formats it to an SQL database

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

engine = create_engine('sqlite:///corum.db')
Base  = declarative_base()

""" -------------------------------------------------------------------------"""
class Taxon(Base):
    """  Represents a species
    Attributes:
        ncbi_id      (:obj:`int`): NCBI id
        species_name (:obj:`str`): name and possibly genetic variant
    """

    __tablename__ = 'taxon'
    ncbi_id       = Column(Integer,primary_key=True)
    swissprot_id  = Column(String(255))

""" -------------------------------------------------------------------------"""
class Observation(Base):
    """  Represents an observation (entries in the original DB)
    Attributes:
        id            (:obj:`int`): internal ID for the observation entry
        cell line     (:obj:`str`): cell line (in whcih the measurement was done)
        pur_method    (:obj:`str`): purification method
        pubmed_id     (:obj:`str`): Pubmed ID of the associated publication
        taxon_ncbi_id (:obj:`str`): NCBI taxonomy id of the organism
    """

    __tablename__    = 'observation'
    id               = Column(Integer, primary_key=True)
    #organism         = Column(String(255))
    cell_line        = Column(String(255))
    pur_method       = Column(String(255))
    pubmed_id        = Column(Integer)

    taxon_ncbi_id    = Column(String(255), ForeignKey('taxon.ncbi_id'))
    taxon            = relationship('Taxon', backref=backref('observation'), foreign_keys=[taxon_ncbi_id])
    #FORMAT: foreign_table = relationship('foreign_class',backref=backref('self_table'),foreign_keys=[self_column])

""" -------------------------------------------------------------------------"""
class Complex(Base):
    """  Represents a protein complex
    Attributes:
        observation_id (:obj:`int`): ID of the  observation
        complex_id     (:obj:`int`): ID of the complex
        complex_name   (:obj:`str`): Complex name
        go_id          (:obj:`str`): GO funtinal annotation
        go_dsc         (:obj:`str`): Description of the annotation
        funcat_id      (:obj:`str`): FUNCAT functional annotation
        funcat_dsc     (:obj:`str`): Description of the annotation
        su_cmt         (:obj:`str`): Subunit comments
        complex_cmt    (:obj:`str`): Compex comments
        disease_cmt    (:obj:`str`): Disease comments
    """

    __tablename__  = 'complex'
    complex_id     = Column(Integer, primary_key=True)
    complex_name   = Column(String(255))
    go_id          = Column(String(255))
    go_dsc         = Column(String(255))
    funcat_id      = Column(String(255))
    funcat_dsc     = Column(String(255))
    su_cmt         = Column(String(255))
    complex_cmt    = Column(String(255))
    disease_cmt    = Column(String(255))

    observation_id = Column(Integer, ForeignKey('observation.id'))
    observation    = relationship('Observation', backref=backref('complex'), foreign_keys=[observation_id])

""" -------------------------------------------------------------------------"""
class Subunit(Base):
    """  Represents subunits of complexes
    Attributes:
        id           (:obj:`int`): Internal subunit ID
        complex_id   (:obj:`int`): ID of the complex to which the subunit belongs
        su_uniprot   (:obj:`int`): UNIPROT ID
        su_entrezs   (:obj:`int`): ENTREZS ID
        protein_name (:obj:`str`): Name of the protein
        gene_name    (:obj:`str`): Gene name
        gene_syn     (:obj:`str`): Synonyms of the gene name
    """

    __tablename__ = 'subunits'
    id            = Column(Integer, primary_key=True)
    su_uniprot    = Column(String(255))
    su_entrezs    = Column(Integer)
    protein_name  = Column(String(255))
    gene_name     = Column(String(255))
    gene_syn      = Column(String(255))

    complex_id    = Column(Integer, ForeignKey('complex.complex_id'))
    complex       = relationship('Complex', backref=backref('subunits'), foreign_keys=[complex_id])

Base.metadata.drop_all(engine)
Base.metadata.create_all(engine)
