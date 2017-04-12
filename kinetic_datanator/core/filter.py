""" Methods to filter and prioritize observations according to multiple criteria

* Taxonomic similarity
* Molecular similarity (similar chemical structure)
* Reaction similarity (similar chemical transformation)
* Genetic perturbation similarity (wildtype / mutant)
* Environmental similarity
    
    * Temperature
    * pH
    * Media

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-10
:Copyright: 2017, Karr Lab
:License: MIT
"""

from . import observation
from kinetic_datanator.util import taxonomy_util
from abc import ABCMeta, abstractmethod
from scipy.stats import norm
from six import with_metaclass
import numpy as np


class FilterRunner(object):
    """ Filter and order a list of observations according to a list of filters.

    Attributes:
        filters (:obj:`list` of :obj:`Filter`): list of filters
    """

    def __init__(self, filters):
        """
        Args:
            filters (:obj:`Filter` or :obj:`list` of :obj:`Filter`): filter or list of filters
        """
        if not isinstance(filters, list):
            filters = [filters]
        self.filters = filters

    def run(self, observations, return_info=False):
        """ Filter and order a list of observations according to a list of filters. Optionally, return additional
        information about the filtering including the scores of the observations and the indices of the prioritized
        observations in the input list of observations.

        1. Calculate the score of each observation for each filter

            * Scores equal to -1, indicate that the observation should be discarded
            * Scores between 0 and 1, indicate how much the observation should be prioritized

        2. Discard any observation which has at least one score equal to -1
        3. Order the observations by their mean score

        Args:
            observations (:obj:`list` of :obj:`observation.Observation`): list of experimental and/or computational observations            
            return_info (:obj:`bool`, optional): if `True`, also return the scores and indices of the ordered observations in the input list

        Returns:
            :obj:`list` of :obj:`observation.Observation` or :obj:`FilterResult`:

                * If `return_info` is `False`: return a list of the observations which matches the filters, ordered by their mean score
                * If `return_info` is `True`: return a list of the observations which matches the filters, ordered by their mean score plus additional diagnostic information
        """

        # score observations against the filters
        all_scores = self.score(observations)

        # filter out observations that must be discarded (observations with score = -1)
        obs, scores, i_obs = self.filter(observations, all_scores)

        # order remaining observations by their mean score
        obs, scores, i_obs = self.order(obs, scores, i_obs)

        # return
        if return_info:
            # return ordered list of observations and additional information
            return FilterResult(observations, self.filters, all_scores, obs, scores, i_obs)
        else:
            # return ordered list of observations
            return obs

    def score(self, observations):
        """ Score observations against the filters

        Args:
            observations (:obj:`list` of :obj:`observation.Observation`): list of experimental and/or computational observations

        Returns:
            :obj:`list` of :obj:`float`: list of scores
        """

        n_obs = len(observations)
        n_filt = len(self.filters)
        scores = np.full((n_obs, n_filt, ), np.nan)
        for i_filter, filter in enumerate(self.filters):
            for i_obs, obs in enumerate(observations):
                scores[i_obs, i_filter] = filter.score(obs)

        return scores

    def filter(self, observations, scores):
        """ Filter out observations that must be discarded (observations with score = -1)

        Args:
            observations (:obj:`list` of :obj:`observation.Observation`): list of experimental and/or computational observations
            scores (:obj:`list` of :obj:`float`): list of scores

        Returns:
            :obj:`tuple`: 

                * :obj:`list` of :obj:`observation.Observation`: list of acceptable observations (observations without scores = -1)
                * :obj:`list` of :obj:`float`: list of scores of the acceptable observations
                * :obj:`list` of :obj:`int`: list of indices of the ordered observations within the original list of observations

        """
        ok_observations = np.extract(np.all(scores >= 0, 1).transpose(), observations).tolist()
        i_ok_observations = np.flatnonzero(np.all(scores >= 0, 1)).tolist()
        ok_scores = scores[i_ok_observations, :]

        return (ok_observations, ok_scores, i_ok_observations, )

    def order(self, observations, scores, i_observations=None):
        """ Order observations by their mean score

        Args:
            observations (:obj:`list` of :obj:`observation.Observation`): list of observations
            scores (:obj:`list` of :obj:`float`): list of scores
            i_observations (:obj:`list` of :obj:`int`, optional): list of indices within the original list of observations

        Returns:
            :obj:`tuple`: 

                * :obj:`list` of :obj:`observation.Observation`: ordered list of observations
                * :obj:`list` of :obj:`float`: list of scores of the ordered observations
                * :obj:`list` of :obj:`int`: list of indices of the ordered observations within the original list of observations
        """

        if not i_observations:
            i_observations = range(len(observations))

        order = np.argsort(np.mean(scores, 1))[::-1]

        ordered_observations = [observations[i] for i in order]
        ordered_scores = scores[order, :]
        i_ordered_observations = [i_observations[i] for i in order]

        return (ordered_observations, ordered_scores, i_ordered_observations, )


class FilterResult(object):
    """ Represents the results of applying a list of filters to a dataset

    Attributes:
        observations (:obj:`list` of `observation.Observation`): input list of observations
        filters (:obj:`list` of `Filter`): list of filters applied to observations
        scores (:obj:`np.ndarray`): matrix of scores (rows: observations in same order as in `observations`; columns: filters, in same orders as in `filters`)
        ordered_observations (:obj:`list` of `observation.Observation`): prioritized list of observations
        ordered_scores (:obj:`np.ndarray`): matrix of scores (rows: observations in same order as in `ordered_observations`; columns: filters, in same orders as in `filters`)
        ordered_observation_indices (:obj:`list` of :obj:`int`): indices of the ordered observations in the input list of observations
    """

    def __init__(self, observations, filters, scores, ordered_observations, ordered_scores, ordered_observation_indices):
        """
        Args:
            observations (:obj:`list` of `observation.Observation`): input list of observations
            filters (:obj:`list` of `Filter`): list of filters applied to observations
            scores (:obj:`np.ndarray`): matrix of scores (rows: observations in same order as in `observations`; columns: filters, in same orders as in `filters`)
            ordered_observations (:obj:`list` of `observation.Observation`): prioritized list of observations
            ordered_scores (:obj:`np.ndarray`): matrix of scores (rows: observations in same order as in `ordered_observations`; columns: filters, in same orders as in `filters`)
            ordered_observation_indices (:obj:`list` of :obj:`int`): indices of the ordered observations in the input list of observations
        """
        self.observations = observations
        self.filters = filters
        self.scores = scores
        self.ordered_observations = ordered_observations
        self.ordered_scores = ordered_scores
        self.ordered_observation_indices = ordered_observation_indices


class Filter(with_metaclass(ABCMeta, object)):
    """ Calculate a numeric score which indicates how well an observation matches one or more criteria.
    Please see :obj:`FilterRunner` to see how these scores are used to filter and order observations.

    Attributes:
        attribute (:obj:`tuple`): list of nested attribute names to score on
    """

    def __init__(self, attribute):
        """
        Args:
            attribute (:obj:`str`): name of attribute to score

        Raises:
            :obj:`ValueError`: if attribute is not defined
        """
        # check that attribute exists
        cls = observation.Observation
        for attr in attribute[0:-1]:
            if attr in cls.Meta.attributes:
                cls = cls.Meta.attributes[attr].related_class
            elif attr in cls.Meta.related_attributes:
                cls = cls.Meta.related_attributes[attr].primary_class
            else:
                raise ValueError('Cannot filter on attribute "{}": Attribute is not defined'.format(attribute))
        if attribute[-1] not in cls.Meta.attributes:
            raise ValueError('Cannot filter on attribute "{}": Attribute is not defined'.format(attribute))

        # store attribute
        self.attribute = attribute

    def get_attribute_value(self, obs):
        """ Get the value of the attribute of observation :obj:`obs`

        Args:
            obs (:obj:`observation.Observation`): observation

        Returns:
            :obj:`object`: value of the attribute of the observation
        """
        val = obs
        for attr in self.attribute:
            val = getattr(val, attr)
        return val

    @abstractmethod
    def score(self, observation):
        """ Calculate a numeric score which indicates how well the observation matches one or more criteria
        Please see :obj:`FilterRunner` to see how these scores are used to filter and order observations.

        Args:
            observation (:obj:`observation.Observation`): experimental and/or computational observation

        Returns:
            :obj:`float`: score which indicates how well the observation matches the criteria
        """
        pass


class OptionsFilter(Filter):
    """ Filters out observations whose attributes have values that are not in a list of acceptable options.

    Attributes:
        attribute (:obj:`str`): name of attribute to score
        options (:obj:`list` of :obj:`object`): list of acceptable values
    """

    def __init__(self, attribute, options):
        """
        Args:
            attribute (:obj:`str`): name of attribute to score
            options (:obj:`list` of :obj:`object`): list of acceptable values
        """
        super(OptionsFilter, self).__init__(attribute)
        self.options = options

    def score(self, obs):
        """ Calculate a numeric score which indicates if the attribute of the observation is one of the acceptable values.

        Args:
            obs (:obj:`observation.Observation`): experimental and/or computational observation

        Returns:
            :obj:`float`: score which indicates if the value of the attribute is one of the acceptable values.
        """
        val = self.get_attribute_value(obs)

        if val in self.options:
            return 1
        return -1


class RangeFilter(Filter):
    """ Filters out observations whose attributes have values that fall outside a specified range.

    Attributes:
        attribute (:obj:`str`): name of attribute to score
        min (:obj:`float`): minimum value
        max (:obj:`float`): maximum value
    """

    def __init__(self, attribute, min=float('nan'), max=float('nan')):
        """
        Args:
            attribute (:obj:`str`): name of attribute to score
            min (:obj:`float`, optional): minimum value
            max (:obj:`float`, optional): maximum value
        """
        super(RangeFilter, self).__init__(attribute)
        self.min = float(min)
        self.max = float(max)

    def score(self, obs):
        """ Calculate a numeric score which indicates if the attribute of the observation falls within the specified range.

        Args:
            obs (:obj:`observation.Observation`): experimental and/or computational observation

        Returns:
            :obj:`float`: score which indicates if the value of the attribute is inside or outside the specified range
        """

        val = self.get_attribute_value(obs)

        if not np.isnan(self.min) and (np.isnan(val) or val < self.min):
            return -1
        if not np.isnan(self.max) and (np.isnan(val) or val > self.max):
            return -1
        return 1


class NormalFilter(Filter):
    """ Prioritizes observations whose attributes have values that are closed to `mean`.

    Attributes:
        attribute (:obj:`str`): Name of an attribute to score
        mean (:obj:`float`): The mean of the distribution. This indicates the value at which the score will be 1.
        std (:obj:`float`): The standard deviation of the distribution. This determines how quickly the score falls to zero away from the mean.
    """

    def __init__(self, attribute, mean, std):
        """
        Args:
            attribute (:obj:`str`): Name of an attribute to score
            mean (:obj:`float`): The mean of the distribution. This indicates the value at which the score will be 1.
            std (:obj:`float`): The standard deviation of the distribution. This determines how quickly the score falls to zero away from the mean.
        """
        super(NormalFilter, self).__init__(attribute)
        self.mean = mean
        self.std = std

    def score(self, obs):
        """ Calculate a numeric score which indicates how well the attribute of the observation matches the specified
        normal distribution (mean, std).

        Args:
            obs (:obj:`observation.Observation`): experimental and/or computational observation

        Returns:
            :obj:`float`: score which indicates how well the observation matches the normal distribution (mean, std)
        """

        val = self.get_attribute_value(obs)

        return 1 - 2 * abs(norm.cdf(val, loc=self.mean, scale=self.std) - 0.5)


class TaxonomicDistanceFilter(Filter):
    """ Prioritizes observations that are from taxonomically close taxa

    Attributes:
        taxon (:obj:`str`): name of the taxon to find data for
        max_dist (:obj:`int`): maximum acceptable number of taxonomic ranks to the latest common ancestor between the 
            target and observed taxa
    """

    def __init__(self, taxon, max_dist=None):
        """
        Args:
            taxon (:obj:`str`): name of the taxon to find data for
            max_dist (:obj:`int`, optional): maximum acceptable number of taxonomic ranks to the latest common ancestor between the target and observed taxa
        """
        super(TaxonomicDistanceFilter, self).__init__(('taxon', 'name', ), )
        self.taxon = taxon

        if max_dist is None:
            taxon = taxonomy_util.Taxon(self.taxon)
            max_dist = taxon.get_max_distance_to_common_ancestor()

        self.max_dist = max_dist

    def score(self, observation):
        """ Score the taxonomic distance from the target taxon to its least common ancestor with the observed taxon.

        Returns:
            :obj:`float`:

                * If the distance to the least common ancestor is greater than `max_dist`, return -1
                * Else, return 1 - {the distance to the least common ancestor} / `max_dist`
        """
        self_taxon = taxonomy_util.Taxon(self.taxon)
        other_taxon = taxonomy_util.Taxon(self.get_attribute_value(observation))

        dist = self_taxon.get_distance_to_common_ancestor(other_taxon)        

        if dist > self.max_dist:
            return -1

        return 1 - dist / self.max_dist


class WildtypeFilter(OptionsFilter):
    """ Filter out observations which were observed for taxa with genetic perturbations """

    def __init__(self):
        super(WildtypeFilter, self).__init__(('taxon', 'perturbations', ), [''])


class ChemicalSimilarityFilter(Filter):
    pass  # todo


class ReactionSimilarityFilter(Filter):
    pass  # todo


class TemperatureRangeFilter(RangeFilter):
    """ Filters out observations with temperatures that fall outside a specified range. """

    def __init__(self, min=float('nan'), max=float('nan')):
        """
        Args:
            min (:obj:`float`, optional): minimum value
            max (:obj:`float`, optional): maximum value
        """
        super(TemperatureRangeFilter, self).__init__(('environment', 'temperature', ), min=min, max=max)


class TemperatureNormalFilter(NormalFilter):
    """ Prioritizes observations with temperatures that are close to `mean`. """

    def __init__(self, mean, std):
        """
        Args:
            mean (:obj:`float`): The mean of the distribution. This indicates the value at which the score will be 1.
            std (:obj:`float`): The standard deviation of the distribution. This determines how quickly the score falls to zero away from the mean.
        """
        super(TemperatureNormalFilter, self).__init__(('environment', 'temperature', ), mean, std)


class PhRangeFilter(RangeFilter):
    """ Filters out observations with pHs that fall outside a specified range. """

    def __init__(self, min=float('nan'), max=float('nan')):
        """
        Args:
            min (:obj:`float`, optional): minimum value
            max (:obj:`float`, optional): maximum value
        """
        super(PhRangeFilter, self).__init__(('environment', 'ph', ), min=min, max=max)


class PhNormalFilter(NormalFilter):
    """ Prioritizes observations with pHs that are close to `mean`. """

    def __init__(self, mean, std):
        """
        Args:
            mean (:obj:`float`): The mean of the distribution. This indicates the value at which the score will be 1.
            std (:obj:`float`): The standard deviation of the distribution. This determines how quickly the score falls to zero away from the mean.
        """
        super(PhNormalFilter, self).__init__(('environment', 'ph', ), mean, std)
