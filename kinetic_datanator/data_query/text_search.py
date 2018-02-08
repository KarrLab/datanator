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
from kinetic_datanator.core import flask_common_schema, models
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

        compound = models.Compound.query.whoosh_search(string).all()
        complex = models.ProteinComplex.query.whoosh_search(string).all()
        subunit = models.ProteinSubunit.query.whoosh_search(string).all()

        db_models = self.rank([compound, complex, subunit])

        return db_models


    def rank(self, lists):
        """
        Ranks based on number of results in the list
        """
        ans = []
        len_lists = [len(l) for l in lists]
        sorted_index = sorted(range(len(len_lists)), key=lambda k: len_lists[k], reverse=True)

        ## Put top results
        for i in sorted_index:
            if lists[i]:
                ans += lists[i][:3]
            else: continue

        ## Fill in bottom results
        for i in sorted_index:
            if lists[i]:
                ans += lists[i][3:]
            else: continue

        return ans
