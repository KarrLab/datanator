""" Search manager for full text searching

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2018-08-08
:Copyright: 2018, Karr Lab
:License: MIT
"""

import abc
import six
import os
# from kinetic_datanator.api.query import reaction_kinetics
from kinetic_datanator.core import models, common_schema
from kinetic_datanator.util.constants import DATA_CACHE_DIR
from kinetic_datanator.api.lib.data_manager import BaseManager
from kinetic_datanator.api.lib.metabolite.manager import metabolite_manager
from kinetic_datanator.api.lib.subunit.manager import subunit_manager
from kinetic_datanator.api.lib.complex.manager import complex_manager

class SearchManager(BaseManager):
    """
    Represents a text search through the database grabbing relevant information

    Attributes:
        db_cache_dirname (:obj:`str`): path location for DB
    """

    def __init__(self, cache_dirname = DATA_CACHE_DIR):
        pass
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

        compound = metabolite_manager._search(string)
        complex = complex_manager._search(string)
        subunit = subunit_manager._search(string)

        # rxns = []
        # for item in compound:
        #     rxn_list = self.q.get_reaction_by_compound(item)
        #     if rxn_list:
        #         rxns.append(rxn_list)
        # reactions = [y for x in rxns for y in x]

        found_dict = {'Compound':compound, 'ProteinSubunit':subunit, 'ProteinComplex':complex}

        return found_dict

search_manager = SearchManager()
