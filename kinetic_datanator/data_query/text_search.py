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
from kinetic_datanator.data_query import dna_protein_interactions, metabolite_concentrations, protein_abundance, protein_protein_interactions, reaction_kinetics
from kinetic_datanator.core import common_schema, models


class TextSearchSession(object):
    """
    Represents a text search through the database grabbing relevant information

    Attributes:
        db_cache_dirname (:obj:`str`): path location for DB
    """

    def __init__(self, db_cache_dirname = os.getcwd()):

        flaskdb = common_schema.CommonSchema(cache_dirname = db_cache_dirname)
        self.q = reaction_kinetics.ReactionKineticsQuery(cache_dirname=db_cache_dirname, include_variants=True)


    def return_search(self, string):
        """
        Collects and creates dictionary of all Database objects coming up in a full
        text search for related string

        Args:
            string (:obj :`str`): item to be searched for

        Returns:
            list_db_models (:obj:`list`): List of ranked collected items from text search
            dict_db_models (:obj:`dict`): Dictionary of collected items from text search
        """

        compound = models.Compound.query.search(string).all()
        complex = models.ProteinComplex.query.search(string).all()
        subunit = models.ProteinSubunit.query.search(string).all()

        rxns = []
        for item in compound:
            rxn_list = self.q.get_reaction_by_compound(item)
            if rxn_list:
                rxns.append(rxn_list)
        reactions = [y for x in rxns for y in x]

        dict_db_models = {'Compound':compound, 'Reaction':reactions, 'ProteinSubunit':subunit, 'ProteinComplex':complex}
        list_db_models = [item for sublist in [compound, reactions, subunit, complex] for item in sublist]

        return list_db_models, dict_db_models
