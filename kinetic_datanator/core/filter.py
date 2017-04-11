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

from abc import ABCMeta, abstractmethod
from scipy.stats import norm
from six import with_metaclass
import numpy as np


class FilterRunner(object):
    """ Filter and order a list of observations according to a list of filters. """

    def run(self, observations, filters, return_info=False):
        """ Filter and order a list of observations according to a list of filters. Optionally, return additional
        information about the filtering including the scores of the observations and the indices of the prioritized
        observations in the input list of observations.

        1. Calculate the score of each observation for each filter

            * Scores equal to -1, indicate that the observation should be discarded
            * Scores between 0 and 1, indicate how much the observation should be prioritized

        2. Discard any observation which has at least one score equal to -1
        3. Order the observations by their mean score

        Args:
            observations (:obj:`list` of :obj:`Observation`): list of experimental and/or computational observations
            filters (:obj:`list` of :obj:`Filter`): list of filters
            return_info (:obj:`bool`, optional): if `True`, also return the scores and indices of the ordered observations in the input list

        Returns:
            :obj:`list` of :obj:`Observation` or :obj:`FilterResult`:

                * If `return_info` is `False`: return a list of the observations which matches the filters, ordered by their mean score
                * If `return_info` is `True`: return a list of the observations which matches the filters, ordered by their mean score plus additional diagnostic information
        """
        n_obs = len(observations)
        n_filt = len(filters)

        # score observations against the filters
        all_scores = np.fill((n_obs, n_filt, ), np.nan)
        for i_filter, filter in enumerate(filters):
            for i_obs, obs in enumerate(observations):
                all_scores[i_obs, i_filter] = filter.score(obs)

        # filter out observations that must be discarded (observations with score = -1)
        i_obs = np.nonzero(np.all(all_scores >= 0, 1))
        obs = np.extract(np.all(all_scores >= 0, 1), observations)
        scores = np.extract(np.all(all_scores >= 0, 1), all_scores)

        # order remaining observations by their mean score
        order = np.argsort(np.mean(scores, 1))

        i_obs = [i_obs[i] for i in order]
        obs = [obs[i] for i in order]
        scores = [scores[i, :] in order]

        # return
        if return_info:
            # return ordered list of observations and additional information
            return FilterResult(observations, filters, all_scores, obs, scores, i_obs)
        else:
            # return ordered list of observations
            return obs


class FilterResult(object):
    """ Represents the results of applying a list of filters to a dataset

    Attributes:
        observations (:obj:`list` of `Observation`): input list of observations
        filters (:obj:`list` of `Filter`): list of filters applied to observations
        scores (:obj:`np.ndarray`): matrix of scores (rows: observations in same order as in `observations`; columns: filters, in same orders as in `filters`)
        ordered_observations (:obj:`list` of `Observation`): prioritized list of observations
        ordered_scores (:obj:`np.ndarray`): matrix of scores (rows: observations in same order as in `ordered_observations`; columns: filters, in same orders as in `filters`)
        ordered_observation_indices (:obj:`list` of :obj:`int`): indices of the ordered observations in the input list of observations
    """

    def __init__(self, observations, filters, scores, ordered_observations, ordered_scores, ordered_observation_indices):
        """
        Args:
            observations (:obj:`list` of `Observation`): input list of observations
            filters (:obj:`list` of `Filter`): list of filters applied to observations
            scores (:obj:`np.ndarray`): matrix of scores (rows: observations in same order as in `observations`; columns: filters, in same orders as in `filters`)
            ordered_observations (:obj:`list` of `Observation`): prioritized list of observations
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
    """

    @abstractmethod
    def score(self, observation):
        """ Calculate a numeric score which indicates how well the observation matches one or more criteria
        Please see :obj:`FilterRunner` to see how these scores are used to filter and order observations.

        Args:
            observation (:obj:`Observation`): experimental and/or computational observation

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
        Attributes:
            attribute (:obj:`str`): name of attribute to score
            options (:obj:`list` of :obj:`object`): list of acceptable values
        """
        self.attribute = attribute
        self.options = options

    def score(self, obs):
        """ Calculate a numeric score which indicates if the attribute of the observation is one of the acceptable values.

        Args:
            obs (:obj:`Observation`): experimental and/or computational observation

        Returns:
            :obj:`float`: score which indicates if the value of the attribute is one of the acceptable values.
        """
        val = obs
        for attr in self.attribute.split('.'):
            val = getattr(val, attr)

        if val in self.options:
            return 1
        return -1


class RangeFilter(Filter):
    """ Filters out observations whose attributes have values that fall outside a specified range.

    Attributes:
        attribute (:obj:`str`): name of attribute to score
        min_val (:obj:`float`): minimum value
        max_val (:obj:`float`): maximum value
    """

    def __init__(self, attribute, min_val=float('nan'), max_val=float('nan')):
        """
        Attributes:
            attribute (:obj:`str`): name of attribute to score
            min_val (:obj:`float`, optional): minimum value
            max_val (:obj:`float`, optional): maximum value
        """
        self.attribute = attribute
        self.min_val = float(min_val)
        self.max_val = float(max_val)

    def score(self, obs):
        """ Calculate a numeric score which indicates if the attribute of the observation falls within the specified range.

        Args:
            obs (:obj:`Observation`): experimental and/or computational observation

        Returns:
            :obj:`float`: score which indicates if the value of the attribute is inside or outside the specified range
        """

        val = obs
        for attr in self.attribute.split('.'):
            val = getattr(val, attr)

        if not np.isnan(self.min_val) and (np.isnan(val) or val < self.min_val):
            return -1
        if not np.isnan(self.max_val) and (np.isnan(val) or val > self.max_val):
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
        self.attribute = attribute
        self.mean = mean
        self.std = std

    def score(self, obs):
        """ Calculate a numeric score which indicates how well the attribute of the observation matches the specified normal distribution (mean, std).

        Args:
            obs (:obj:`Observation`): experimental and/or computational observation

        Returns:
            :obj:`float`: score which indicates how well the observation matches the normal distribution (mean, std)
        """

        val = obs
        for attr in self.attribute.split('.'):
            val = getattr(val, attr)

        return 1 - 2 * abs(norm.cdf(val, loc=self.mean, scale=self.std) - 0.5)


class TaxonomicDistanceFilter(Filter):
    """ Prioritizes observations that are from taxonomically close strains

    Attributes:
        species (:obj:`str`): name of the species to find data for
        max_dist (:obj:`int`): maximum acceptable number of taxonomic ranks to the latest common ancestor between the target and observed species
    """

    def __init__(self, species, max_dist=8):
        """
        Args:
            species (:obj:`str`): name of the species to find data for
            max_dist (:obj:`int`, optional): maximum acceptable number of taxonomic ranks to the latest common ancestor between the target and observed species
        """
        self.species = species
        self.max_dist = max_dist

    def score(self, observation):
        """ Score the taxonomic distance from the target species to its least common ancestor with the observed species.

        Returns:
            :obj:`float`:

                * If the distance to the least common ancestor is greater than `max_dist`, return -1
                * Else, return 1 - {the distance to the least common ancestor} / 8
        """
        tax_ranks_to_lca = 0.  # todo: calculate

        if tax_ranks_to_lca > max_dist:
            return -1

        return 1 - tax_ranks_to_lca / 8


class WildtypeFilter(OptionsFilter):
    """ Filter out observations which were observed for strains with genetic perturbations """

    def __init__(self):
        super(WildtypeFilter, self).__init__('strain.perturbations', [''])


class ChemicalSimilarityFilter(Filter):
    pass  # todo


class ReactionSimilarityFilter(Filter):
    pass  # todo


class TemperatureRangeFilter(RangeFilter):
    """ Filters out observations with temperatures that fall outside a specified range. """

    def __init__(self, min_val=float('nan'), max_val=float('nan')):
        """
        Attributes:
            min_val (:obj:`float`, optional): minimum value
            max_val (:obj:`float`, optional): maximum value
        """
        super(TemperatureRangeFilter, self).__init__('environment.temperature', min_val=min_val, max_val=max_val)


class TemperatureNormalFilter(NormalFilter):
    """ Prioritizes observations with temperatures that are close to `mean`. """

    def __init__(self, mean, std):
        """
        Args:
            mean (:obj:`float`): The mean of the distribution. This indicates the value at which the score will be 1.
            std (:obj:`float`): The standard deviation of the distribution. This determines how quickly the score falls to zero away from the mean.
        """
        super(TemperatureNormalFilter, self).__init__('environment.temperature', mean, std)


class PhRangeFilter(RangeFilter):
    """ Filters out observations with pHs that fall outside a specified range. """

    def __init__(self, min_val=float('nan'), max_val=float('nan')):
        """
        Attributes:
            min_val (:obj:`float`, optional): minimum value
            max_val (:obj:`float`, optional): maximum value
        """
        super(PhRangeFilter, self).__init__('environment.ph', min_val=min_val, max_val=max_val)


class PhNormalFilter(NormalFilter):
    """ Prioritizes observations with pHs that are close to `mean`. """

    def __init__(self, mean, std):
        """
        Args:
            mean (:obj:`float`): The mean of the distribution. This indicates the value at which the score will be 1.
            std (:obj:`float`): The standard deviation of the distribution. This determines how quickly the score falls to zero away from the mean.
        """
        super(PhNormalFilter, self).__init__('environment.ph', mean, std)
