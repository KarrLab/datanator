"""
This module downloads the JASPAR database of transcription factor binding motifs (http://jaspar.genereg.net/)
via a seris of text files, parses them, and stores them in an SQLlite database.

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-08-01
:Copyright: 2017, Karr Lab
:License: MIT
"""

from datanator.core import data_source
from sqlalchemy import Column, Integer, String,Numeric, ForeignKey
from sqlalchemy.orm import relationship, backref
import requests
import six
import sqlalchemy.ext.declarative
import tarfile
import os

Base = sqlalchemy.ext.declarative.declarative_base()

class Matrix(Base):

    __tablename__ = 'MATRIX'

    ID = Column(Integer, primary_key = True)
    COLLECTION = Column(String(255))
    BASE_ID = Column(String(255))
    VERSION = Column(Integer)
    NAME = Column(String(255))

class Annotation(Base):
    __tablename__ = 'MATRIX_ANNOTATION'

    ID = Column(Integer, primary_key = True)
    TAG = Column(String(255))
    VAL =  Column(String(255), primary_key = True)

class Data(Base):
    __tablename__ = 'MATRIX_DATA'

    ID = Column(Integer, primary_key = True)
    row = Column(String(255), primary_key = True)
    col = Column(Integer, primary_key = True)
    val = Column(Integer)

class Protein(Base):
    __tablename__ = 'MATRIX_PROTEIN'

    ID = Column(Integer, primary_key = True)
    ACC = Column(String(255), primary_key = True)

class Species(Base):
    __tablename__ = 'MATRIX_SPECIES'

    ID = Column(Integer, primary_key = True)
    TAX_ID = Column(Integer, primary_key = True)

class Taxon(Base):
    __tablename__ = 'TAX'

    TAX_ID = Column(Integer, primary_key = True)
    SPECIES= Column(String(255), primary_key = True)

class TaxonExtension(Base):
    __tablename__ = 'TAX_EXT'

    TAX_ID = Column(Integer, primary_key = True)
    NAME = Column(String(255), primary_key = True)

class Tffm(Base):
    __tablename__ = 'TFFM'

    ID =  Column(Integer, primary_key = True)
    BASE_ID = Column(String(255))
    VERSION = Column(Integer)
    MATRIX_BASE_ID = Column(String(255))
    MATRIX_VERSION = Column(Integer)
    NAME = Column(String(255))
    LOG_P_1ST_ORDER = Column(Numeric)
    LOG_P_DETAILED = Column(Numeric)
    EXPERIMENT_NAME = Column(String(255))

class Jaspar(data_source.HttpDataSource):
    """ A local SQLite copy of the `JASPAR <http://jaspar.genereg.net>`_ database of transcription factor binding profiles """
    base_model = Base
    ENDPOINT_DOMAINS = {
        'jaspar': 'http://jaspar.genereg.net/download/database/JASPAR2018.sqlite.tar.gz',
    }

    def load_content(self):
        database_url = self.ENDPOINT_DOMAINS['jaspar']
        req = self.requests_session
        response = req.get(database_url)
        response.raise_for_status()
        tarred = tarfile.TarFile(fileobj = six.BytesIO(response.content))
        tarred.extractall(self.cache_dirname)
        os.rename(self.cache_dirname+'/JASPAR2018.sqlite', self.cache_dirname+'/Jaspar.sqlite')

        self.session.commit()
