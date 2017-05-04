# -*- coding: utf-8 -*-

"""
:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-06
:Copyright: 2017, Karr Lab
:License: MIT
"""

from .util import molecule_util
from .util import reaction_util
from .util import taxonomy_util
import numpy


class Query:
    """ Represents a SABIO-RK query

    Attributes:
        max_taxon_dist (:obj:`float`): maximum taxonomic distance from the target taxon to its latest
            common ancestor with the observed taxon
        include_mutants (:obj:`bool`): if :obj:`True`, include observations from mutants
        min_temp (:obj:`float`): minimum observed temperature
        max_temp (:obj:`float`): maximum observed temperature
        min_ph (:obj:`float`): minimum observed pH
        max_ph (:obj:`float`): maximum observed pH
    """

    def __init__(self, max_taxon_dist=None, include_mutants=False,
                 min_temp=None, max_temp=None, min_ph=None, max_ph=None):
        """
        Args:
            max_taxon_dist (:obj:`float`, optional): maximum taxonomic distance from the target taxon to its
                latest common ancestor with the observed taxon
            include_mutants (:obj:`bool`, optional): if :obj:`True`, include observations from mutants
            min_temp (:obj:`float`, optional): minimum observed temperature
            max_temp (:obj:`float`, optional): maximum observed temperature
            min_ph (:obj:`float`, optional): minimum observed pH
            max_ph (:obj:`float`, optional): maximum observed pH
        """
        self.max_taxon_dist = max_taxon_dist
        self.include_mutants = include_mutants
        self.min_temp = min_temp
        self.max_temp = max_temp
        self.min_ph = min_ph
        self.max_ph = max_ph

    def run(self, reaction):
        """
        Args:
            reaction (:obj:`reaction_util.Reaction`): reaction to find kinetic data for
        """

        """ search for relevant kinetic laws """
        results = []

        # search for reactions by their participants

        # search for reactions by their EC numbers

        """ filter out the most relevant kinetic laws """
        self._filter_results(results)

        """ calculate a statistical representation of the most relevant rate laws """
        consensus = self._calc_consensus(results)

        return (consensus, results)

    def _make_query(self, reaction, query_participants=True):
        """ Create a query for the SABIO-RK web service

        Args:
            reaction (:obj:`reaction_util.Reaction`): reaction to find kinetic data for
            query_participants (:obj:`bool`, optional): if :obj:`True`, query reactions by the names of their
                participants. Otherwise, query reactions by their EC number(s)
        """
        query_fragments = [
            self._make_query_mutant(),
            self._make_query_temperature(),
            self._make_query_ph(),
        ]
        if query_participants:
            query_fragments.append(self._make_query_participant(reaction))
        else:
            query_fragments.append(self._make_query_ec_number(reaction))

        return ' AND '.join(filter(lambda x: x, query_fragments))

    def _parse_results(self, text):
        compartments, compounds, reactions, kinetic_laws = SabioRkIo().parse_kinetic_laws(text)
        return kinetic_laws

    def _filter_results(self, results):
        # same reaction
        # todo: same participants, excluding H+, OH-

        # taxonomic distance
        base_taxon = taxonomy_util.Taxon(base_species)
        for result in reversed(results):
            entry_taxon = taxonomy_util.Taxon(result.taxon)
            entry.taxon_dist = base_taxon.get_distance_to_common_ancestor(entry_taxon)

            if entry.taxon_dist > max_taxon_dist:
                results.remove()

        return results

    def _calc_consensus(self, results):
        values = {}
        for result in results:
            for param in result.parameters:
                if (param.type, param.subtype) not in values:
                    values[(param.type, param.subtype)] = []
                values[(param.type, param.subtype)].append((param.value, param.units))

        consensus = []
        for (type, subtype), vals_units in values.items():
            vals = [vu[0] for vu in vals_units]
            units = [vu[1] for vu in vals_units]

            consensus.append({
                'type': type,
                'subtype': subtype,
                'min': min(vals),
                'max': max(vals),
                'mean': numpy.mean(vals),
                'std': numpy.std(vals),
            })

        return consensus
