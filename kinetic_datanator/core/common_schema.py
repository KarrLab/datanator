# -*- coding: utf-8 -*-

# """
# This code is a common schema for all the kinetic_datanator modules
#
# :Author: Saahith Pochiraju <saahith116@gmail.com>
# :Date: 2017-07-31
# :Copyright: 2017, Karr Lab
# :License: MIT
# """
#
# ##TODO: Create the definition of tables for the format and then include with the query function in sql
#
# from sqlalchemy import Column, BigInteger, Integer, String, Text, ForeignKey, Table, create_engine
# from sqlalchemy.orm import relationship, backref, sessionmaker
# from kinetic_datanator.core import data_source
# import sqlalchemy.ext.declarative
# from six import BytesIO
# import six
#
# Base = sqlalchemy.ext.declarative.declarative_base()
#
#
# class Taxon(Base):
#     """
#
#
#
#     """
