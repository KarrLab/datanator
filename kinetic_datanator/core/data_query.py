"""
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-05-16
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import filter
from kinetic_datanator.core import data_model
import abc
import kinetic_datanator.core.data_source
import math
import numpy
import pint
import six

unit_registry = pint.UnitRegistry()


class DataQueryEngine(six.with_metaclass(abc.ABCMeta, object)):
    """ Represents a query of a data source

    1. Find observations for the exact or similar model components

    2. Filter out observations from disimilar genetic and environmental conditions and 
      rank the remaing observations by their similarity to the desired genetic and environmental 
      conditions

      * Taxonomy
      * Genetic variation (wildtype/mutant)
      * Temperature
      * pH

    3. Calculate a statistical representation of the relevant observations

    Attributes:
        filters (:obj:`list` of :obj:`filter.Filter`): list of filters        
    """

    def __init__(self,
                 taxon=None, max_taxon_dist=None, taxon_dist_scale=None, include_variants=False,
                 temperature=37., temperature_std=1.,
                 ph=7.5, ph_std=0.3):
        """
        Args:
            taxon (:obj:`str`, optional): target taxon
            max_taxon_dist (:obj:`int`, optional): maximum taxonomic distance to include
            taxon_dist_scale (:obj:`float`, optional): The scale of the taxonomic distance scoring distribution. 
                This determines how quickly the score falls to zero away from zero.
            include_variants (:obj:`bool`, optional): if :obj:`True`, also include observations from mutant taxa
            temperature (:obj:`float`, optional): desired temperature to search for
            temperature_std (:obj:`float`, optional): how much to penalize observations from other temperatures
            ph (:obj:`float`, optional): desired pH to search for
            ph_std (:obj:`float`, optional): how much to penalize observations from other pHs
        """
        # todo: filter media

        self.filters = filters = []

        if taxon:
            filters.append(filter.TaxonomicDistanceFilter(taxon, max=max_taxon_dist, scale=taxon_dist_scale))

        if not include_variants:
            filters.append(filter.WildtypeFilter())

        if not math.isnan(temperature) and not math.isnan(temperature_std):
            filters.append(filter.TemperatureNormalFilter(temperature, temperature_std))

        if not math.isnan(ph) and not math.isnan(ph_std):
            filters.append(filter.PhNormalFilter(ph, ph_std))

    def run(self, component):
        """ 

        1. Find observations for the exact or similar model components and genetic and environmental conditions
        2. Rank the results by their similarity to the model component and the genetic and environmental conditions
        3. Calculate a consensus statistical representation of the relevant observations

        Args:
            component (:obj:`str`): model component to find data for

        Returns:
            :obj:`list` of :obj:`data_model.Consensus`: statistical consensus of the relevant observations of 
                :obj:`component` and the observations it was based on
        """
        obs = self.get_observed_values(component)
        filter_result = self.filter_observed_values(component, obs)
        return self.get_consensus(component, filter_result)

    @abc.abstractmethod
    def get_observed_values(self, component):
        """ Find the observations relevant to :obj:`component`

        Args:
            component (:obj:`str`): model component to find data for

        Returns:
            :obj:`list` of :obj:`data_model.ObservedValue`: list of relevant observations
        """
        pass

    def filter_observed_values(self, component, observations):
        """ Filter out observations from dissimilar genetic and environmental conditions and
        order the remaining observations by their similarity to specified genetic and 
        environmental conditions.

        Args:
            component (:obj:`str`): model component to find data for
            observations (:obj:`list` of :obj:`data_model.ObservedValue`): list of observations

        Returns:
            :obj:`filter.FilterResult`: filter result
        """
        return filter.FilterRunner(self.filters) \
            .run(observations, return_info=True)

    def get_consensus(self, component, filter_result):
        """ Calculate a consensus statistical representation of the one or more observations

        Args:
            component (:obj:`str`): model component to find data for
            filter_result (:obj:`filter.FilterResult`): filter result

        Returns:
            :obj:`list` of :obj:`data_model.Consensus`: statistical consensus of the relevant observations of 
                :obj:`component` and the observations it was based on
        """
        # group observations by their subcomponents, attributes
        grouped_obs = {}
        for obs, score in zip(filter_result.observations, filter_result.scores):
            if obs.attribute not in grouped_obs:
                grouped_obs[obs.attribute] = {}
            if obs.subcomponent:
                if obs.subcomponent not in grouped_obs[obs.attribute]:
                    grouped_obs[obs.attribute][obs.subcomponent] = []
                grouped_obs[obs.attribute][obs.subcomponent].append((obs, score))
            else:
                if '__default__' not in grouped_obs[obs.attribute]:
                    grouped_obs[obs.attribute]['__default__'] = []
                grouped_obs[obs.attribute]['__default__'].append((obs, score))

        # calculate weighted mean
        consensus = []
        for attribute, attribute_obs in grouped_obs.items():
            for subcomponent, subcomponent_obs in attribute_obs.items():
                obs = []
                scores = []
                values = []
                errors = []
                for ob, score in subcomponent_obs:
                    obs.append(data_model.Evidence(observation=ob, relevance=score))
                    scores.append(score)

                    value = ob.value * unit_registry(ob.units)
                    error = ob.error * unit_registry(ob.units)
                    values.append(value.to_base_units().magnitude)
                    errors.append(error.to_base_units().magnitude)

                value = numpy.average(values, weights=scores)
                error = numpy.sqrt(numpy.average(numpy.power(errors, 2), weights=scores))

                consensus.append(data_model.Consensus(
                    component=component,
                    attribute=attribute,
                    value=value,
                    error=error,
                    units=units,
                    method=data_model.ConsensusMethod.weighted_mean,
                    evidence=obs))

        return consensus


class CachedDataSourceQueryEngine(DataQueryEngine):
    """ Represents a query of a cached data source

    Attributes:
        data_source (:obj:`kinetic_datanator.core.data_source.CachedDataSource`): cached data source
    """

    def __init__(self,
                 taxon=None, max_taxon_dist=None, taxon_dist_scale=None, include_variants=False,
                 temperature=37., temperature_std=1.,
                 ph=7.5, ph_std=0.3,
                 data_source=None):
        """
        Args:            
            taxon (:obj:`str`, optional): target taxon
            max_taxon_dist (:obj:`int`, optional): maximum taxonomic distance to include
            taxon_dist_scale (:obj:`float`, optional): The scale of the taxonomic distance scoring distribution. 
                This determines how quickly the score falls to zero away from zero.
            include_variants (:obj:`bool`, optional): if :obj:`True`, also include observations from mutant taxa
            temperature (:obj:`float`, optional): desired temperature to search for
            temperature_std (:obj:`float`, optional): how much to penalize observations from other temperatures
            ph (:obj:`float`, optional): desired pH to search for
            ph_std (:obj:`float`, optional): how much to penalize observations from other pHs
            data_source (:obj:`kinetic_datanator.core.data_source.CachedDataSource`, optional): cached data source
        """
        super(CachedDataSourceQueryEngine, self).__init__(
            taxon=taxon, max_taxon_dist=max_taxon_dist, taxon_dist_scale=taxon_dist_scale, include_variants=include_variants,
            temperature=temperature, temperature_std=temperature_std,
            ph=ph, ph_std=ph_std)

        self.data_source = data_source
