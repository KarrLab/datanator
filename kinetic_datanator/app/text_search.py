#-*- coding: utf-8 -*-

""" Test of common schema

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-12-18
:Copyright: 2017, Karr Lab
:License: MIT
"""

import abc
import six
import os
from kinetic_datanator.data_query import dna_protein_interactions, metabolite_concentrations, protein_concentrations, protein_protein_interactions, reaction_kinetics
from kinetic_datanator.app import flask_common_schema, models
import flask_whooshalchemy


class TextSearchSession(six.with_metaclass(abc.ABCMeta, object)):
    """
    Represents a text search through the database grabbing relevant information

    Attributes:
        db_cache_dirname (:obj:`str`): path location for DB
    """

    def __init__(self, db_cache_dirname = os.getcwd()):

        flaskdb = flask_common_schema.FlaskCommonSchema(cache_dirname = db_cache_dirname)

        for item in flaskdb.text_indicies:
            flask_whooshalchemy.whoosh_index(models.app, item)

    def return_search(self, string):
        """
        Collects and creates dictionary of all Database objects coming up in a full
        text search for related string

        Args:
            string (:obj:`str`): item to be searched for

        Returns:
            db_models (:obj:`dict`): Dictionary of collected items from text search

        """

        db_models = {}
        db_models['Compound'] = models.Compound.query.whoosh_search(string).all()
        db_models['ProteinComplex'] = models.ProteinComplex.query.whoosh_search(string).all()
        db_models['ProteinInteractions'] = models.ProteinInteractions.query.whoosh_search(string).all()
        db_models['Taxon'] = models.Taxon.query.whoosh_search(string).all()
        db_models['Synonym'] = models.Synonym.query.whoosh_search(string).all()
        db_models['CellLine'] = models.CellLine.query.whoosh_search(string).all()
        db_models['CellCompartment'] = models.CellCompartment.query.whoosh_search(string).all()
        db_models['ProteinSubunit'] = models.ProteinSubunit.query.whoosh_search(string).all()

        return db_models
