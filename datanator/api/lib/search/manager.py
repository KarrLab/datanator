""" Search manager for full text searching

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2018-08-08
:Copyright: 2018, Karr Lab
:License: MIT
"""

import abc
import six
import os
# from datanator.api.query import reaction_kinetics
from datanator.core import models, common_schema
from datanator.util.constants import DATA_CACHE_DIR, METABOLITE_REACTION_LIMIT
from datanator.api.lib.data_manager import BaseManager
from datanator.api.lib.metabolite.manager import metabolite_manager
from datanator.api.lib.subunit.manager import subunit_manager
from datanator.api.lib.complex.manager import complex_manager
from datanator.api.lib.reaction.manager import reaction_manager

class SearchManager(BaseManager):
    """
    Represents a text search through the database grabbing relevant information

    Attributes:
        db_cache_dirname (:obj:`str`): path location for DB
    """

    def __init__(self, cache_dirname = DATA_CACHE_DIR):
        self.data_source = common_schema.CommonSchema(cache_dirname=cache_dirname)
        # self.data_source = common_schema.CommonSchema(cache_dirname= cache_dirname)
        # self.q = reaction_kinetics.ReactionKineticsQuery(cache_dirname=cache_dirname, include_variants=True)

    def search(self, string):
        """
        Collects and creates dictionary of all Database objects coming up in a full
        text search for related string

        Args:
            string (:obj :`str`): item to be searched for

        Returns:
            list_db_models (:obj:`list`): List of ranked collected items from text search
            dict_db_models (:obj:`dict`): Dictionary of collected items from text search
        """

        metabolite = metabolite_manager._search_simple(string)
        complex = complex_manager._search(string)
        subunit = subunit_manager._search_simple(string)

        rxns = []
        print(metabolite[:METABOLITE_REACTION_LIMIT])
        for _metabolite in metabolite[:METABOLITE_REACTION_LIMIT]:
            rxn_list = reaction_manager.get_reaction_by_metabolite(_metabolite)
            if rxn_list:
                rxns.append(rxn_list)

        reactions = [y for x in rxns for y in x]

        found_dict = {'Metabolite':metabolite, 'ProteinSubunit':subunit, 'ProteinComplex':complex, 'Reaction': reactions}

        return found_dict

    def get_object_by_id(self, id):
        return self.data_source.session.query(models.Observation).get(id)

search_manager = SearchManager()
