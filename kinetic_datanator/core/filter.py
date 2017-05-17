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
from kinetic_datanator.util import molecule_util
from kinetic_datanator.util import taxonomy_util
from six import with_metaclass
import math
import numpy
import scipy.stats


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
        scores = numpy.full((n_obs, n_filt, ), numpy.nan)
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
        ok_observations = numpy.extract(numpy.all(scores >= 0, 1).transpose(), observations).tolist()
        i_ok_observations = numpy.flatnonzero(numpy.all(scores >= 0, 1)).tolist()
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

        order = numpy.argsort(numpy.mean(scores, 1))[::-1]

        ordered_observations = [observations[i] for i in order]
        ordered_scores = scores[order, :]
        i_ordered_observations = [i_observations[i] for i in order]

        return (ordered_observations, ordered_scores, i_ordered_observations, )


class FilterResult(object):
    """ Represents the results of applying a list of filters to a dataset

    Attributes:
        observations (:obj:`list` of `observation.Observation`): input list of observations
        filters (:obj:`list` of `Filter`): list of filters applied to observations
        scores (:obj:`numpy.ndarray`): matrix of scores (rows: observations in same order as in `observations`; columns: filters, in same orders as in `filters`)
        ordered_observations (:obj:`list` of `observation.Observation`): prioritized list of observations
        ordered_scores (:obj:`numpy.ndarray`): matrix of scores (rows: observations in same order as in `ordered_observations`; columns: filters, in same orders as in `filters`)
        ordered_observation_indices (:obj:`list` of :obj:`int`): indices of the ordered observations in the input list of observations
    """

    def __init__(self, observations, filters, scores, ordered_observations, ordered_scores, ordered_observation_indices):
        """
        Args:
            observations (:obj:`list` of `observation.Observation`): input list of observations
            filters (:obj:`list` of `Filter`): list of filters applied to observations
            scores (:obj:`numpy.ndarray`): matrix of scores (rows: observations in same order as in `observations`; columns: filters, in same orders as in `filters`)
            ordered_observations (:obj:`list` of `observation.Observation`): prioritized list of observations
            ordered_scores (:obj:`numpy.ndarray`): matrix of scores (rows: observations in same order as in `ordered_observations`; columns: filters, in same orders as in `filters`)
            ordered_observation_indices (:obj:`list` of :obj:`int`): indices of the ordered observations in the input list of observations
        """
        self.observations = observations
        self.filters = filters
        self.scores = scores
        self.ordered_observations = ordered_observations
        self.ordered_scores = ordered_scores
        self.ordered_observation_indices = ordered_observation_indices


class Filter(object):
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
        self.attribute = attribute

    def get_attribute_value(self, observation):
        """ Get the value of the attribute of observation :obj:`observation`

        Args:
            observation (:obj:`observation.Observation`): observation

        Returns:
            :obj:`object`: value of the attribute of the observation
        """
        val = observation
        for attr in self.attribute:
            val = getattr(val, attr)
        return val

    def transform_attribute_value(self, value):
        """ Transform an attribute value

        Args:
            value (:obj:`object`): raw value

        Returns:
            :obj:`object`: transformed value
        """
        return value

    def score(self, observation):
        """ Calculate a scaled numeric score betwen 0 and 1 which indicates how well the observation matches 
        one or more criteria. Please see :obj:`FilterRunner` to see how these scores are used to filter and 
        order observations.

        Args:
            observation (:obj:`observation.Observation`): experimental and/or computational observation

        Returns:
            :obj:`float`: score which indicates how well the observation matches the criteria
        """
        return self.transform_attribute_value(self.get_attribute_value(observation))


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

    def score(self, observation):
        """ Calculate a numeric score which indicates how well the observation matches 
        one or more criteria. Please see :obj:`FilterRunner` to see how these scores are used to filter and 
        order observations.

        Args:
            observation (:obj:`observation.Observation`): experimental and/or computational observation

        Returns:
            :obj:`float`: score which indicates how well the observation matches the criteria
        """
        val = self.transform_attribute_value(self.get_attribute_value(observation))
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

    def get_raw_score(self, observation):
        """ Calculate a numeric score which indicates how well the observation matches one or more criteria
        Please see :obj:`FilterRunner` to see how these scores are used to filter and order observations.

        Args:
            observation (:obj:`observation.Observation`): experimental and/or computational observation

        Returns:
            :obj:`float`: score which indicates how well the observation matches the criteria
        """
        return self.get_attribute_value(observation)

    def score(self, observation):
        """ Calculate a scaled numeric score betwen 0 and 1 which indicates how well the observation matches 
        one or more criteria. Please see :obj:`FilterRunner` to see how these scores are used to filter and 
        order observations.

        Args:
            observation (:obj:`observation.Observation`): experimental and/or computational observation

        Returns:
            :obj:`float`: score which indicates how well the observation matches the criteria
        """
        val = self.transform_attribute_value(self.get_attribute_value(observation))

        if not numpy.isnan(self.min) and (numpy.isnan(val) or val < self.min):
            return -1
        if not numpy.isnan(self.max) and (numpy.isnan(val) or val > self.max):
            return -1
        return 1


class NormalFilter(Filter):
    """ Prioritizes observations whose attributes have values that are closed to `mean`.

    Attributes:
        attribute (:obj:`str`): Name of an attribute to score
        mean (:obj:`float`): The mean of the distribution. This indicates the value at which the score will be 1.
        std (:obj:`float`): The standard deviation of the distribution. This determines how quickly the score falls to zero away from the mean.
    """

    def __init__(self, attribute, mean=0., std=1.):
        """
        Args:
            attribute (:obj:`str`): Name of an attribute to score
            mean (:obj:`float`, optional): The mean of the distribution. This indicates the value at which the score will be 1.
            std (:obj:`float`, optional): The standard deviation of the distribution. This determines how quickly the score falls to zero away from the mean.
        """
        super(NormalFilter, self).__init__(attribute)
        self.mean = mean
        self.std = std

    def score(self, observation):
        """ Calculate a numeric score which indicates how well the attribute of the observation matches the specified
        normal distribution (mean, std).

        Args:
            observation (:obj:`observation.Observation`): experimental and/or computational observation

        Returns:
            :obj:`float`: score which indicates how well the observation matches the normal distribution (mean, std)
        """
        val = self.transform_attribute_value(self.get_attribute_value(observation))

        return 1 - 2 * abs(scipy.stats.norm.cdf(val, loc=self.mean, scale=self.std) - 0.5)


class ExponentialFilter(Filter):
    """ Prioritizes observations based on an exponential scale

    Attributes:
        attribute (:obj:`str`): Name of an attribute to score
        center (:obj:`float`): The center of the distribution. This indicates the value at which the score will be 1.
        scale (:obj:`float`): The scale of the distribution. This determines how quickly the score falls to zero away from the center.
    """

    def __init__(self, attribute, center=0., scale=1.):
        """
        Args:
            attribute (:obj:`str`): Name of an attribute to score
            center (:obj:`float`, optional): The center of the distribution. This indicates the value at which the score will be 1.
            scale (:obj:`float`, optional): The scale of the distribution. This determines how quickly the score falls to zero away from the center.
        """
        super(ExponentialFilter, self).__init__(attribute)
        self.center = center
        self.scale = scale

    def score(self, observation):
        """ Calculate a numeric score which indicates how well the attribute of the observation matches the specified
        distribution (center, scale).

        Args:
            observation (:obj:`observation.Observation`): experimental and/or computational observation

        Returns:
            :obj:`float`: score which indicates how well the observation matches the distribution (center, scale)
        """
        val = self.transform_attribute_value(self.get_attribute_value(observation))

        return math.exp(-(val - self.center) / self.scale)


class TaxonomicDistanceFilter(ExponentialFilter):
    """ Prioritizes observations that are from taxonomically close taxa

    Attributes:
        taxon (:obj:`str`): name of the taxon to find data for
        scale (:obj:`float`): The scale of the distribution. This determines how quickly the score falls to zero away from the center.
        max_dist (:obj:`float`): maximum distance to the latest common ancestor with the observed taxon
    """

    def __init__(self, taxon, scale=float('nan'), max_dist=float('nan')):
        """
        Args:
            taxon (:obj:`str`): name of the taxon to find data for
            scale (:obj:`float`, optional): The scale of the distribution. This determines how quickly the score falls to zero away from the center.
            max_dist (:obj:`float`, optional): maximum distance to the latest common ancestor with the observed taxon
        """

        if numpy.isnan(scale):
            taxon_obj = taxonomy_util.Taxon(name=taxon)
            scale = (taxon_obj.get_max_distance_to_common_ancestor() - 2) / 5

        if numpy.isnan(max_dist):
            taxon_obj = taxonomy_util.Taxon(name=taxon)
            max_dist = taxon_obj.get_max_distance_to_common_ancestor() - 2

        super(TaxonomicDistanceFilter, self).__init__(('taxon', 'name', ))

        self.taxon = taxon
        self.scale = scale
        self.max_dist = max_dist

    def transform_attribute_value(self, taxon):
        """ Transform an attribute value

        Args:
            taxon (:obj:`str`): observed taxon

        Returns:
            :obj:`int`: distance to latest common ancestor with the observed taxon
        """
        self_taxon = taxonomy_util.Taxon(name=self.taxon)
        obs_taxon = taxonomy_util.Taxon(name=taxon)

        return self_taxon.get_distance_to_common_ancestor(obs_taxon)

    def score(self, observation):
        """ Calculate a scaled numeric score betwen 0 and 1 which indicates how well the observation matches 
        one or more criteria. Please see :obj:`FilterRunner` to see how these scores are used to filter and 
        order observations.

        Args:
            observation (:obj:`observation.Observation`): experimental and/or computational observation

        Returns:
            :obj:`float`: score which indicates how well the observation matches the criteria
        """
        val = self.transform_attribute_value(self.get_attribute_value(observation))

        if not numpy.isnan(self.max_dist) and (numpy.isnan(val) or val > self.max_dist):
            return -1

        return math.exp(-val / self.scale)


class WildtypeFilter(OptionsFilter):
    """ Filter out observations which were observed for taxa with genetic perturbations """

    def __init__(self):
        super(WildtypeFilter, self).__init__(('taxon', 'perturbations', ), [''])


class ComponentMolecularSimilarityFilter(Filter):
    """ Score similarity with observed molecules

    Attributes:
        structure (:obj:`str`): 
    """

    def __init__(self, structure):
        """
        Args:
            structure (:obj:`str`): structure
        """
        super(ComponentMolecularSimilarityFilter, self).__init__(('component', 'structure', ), )
        self.structure = structure

    def transform_attribute_value(self, structure):
        """ Transform an attribute value

        Args:
            structure (:obj:`str`): observed structure

        Returns:
            :obj:`float`: similarity between the observed and query structure
        """
        self_mol = molecule_util.Molecule(structure=self.structure)
        obs_mol = molecule_util.Molecule(structure=structure)
        return self_mol.get_similarity(obs_mol)


class MolecularSimilarityFilter(Filter):
    pass  # todo


class ComponentReactionSimilarityFilter(Filter):
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
