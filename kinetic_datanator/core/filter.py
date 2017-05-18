""" Methods to filter and prioritize observed values according to multiple criteria

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

from kinetic_datanator.core import data_model
from kinetic_datanator.util import molecule_util
from kinetic_datanator.util import taxonomy_util
from six import with_metaclass
from wc_utils.util import string
import math
import numpy
import scipy.stats


class FilterRunner(object):
    """ Filter and order a list of observed values according to a list of filters.

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

    def run(self, observed_values, return_info=False):
        """ Filter and order a list of observed values according to a list of filters. Optionally, return additional
        information about the filtering including the scores of the observed values and the indices of the prioritized
        observed values in the input list of observed values.

        1. Calculate the score of each observed value for each filter

            * Scores equal to -1, indicate that the observed value should be discarded
            * Scores between 0 and 1, indicate how much the observed value should be prioritized

        2. Discard any observed value which has at least one score equal to -1
        3. Order the observed values by their mean score

        Args:
            observed_values (:obj:`list` of :obj:`data_model.ObservedValue`): list of experimental and/or computational observed values            
            return_info (:obj:`bool`, optional): if `True`, also return the scores and indices of the ordered observed values in the input list

        Returns:
            :obj:`list` of :obj:`data_model.ObservedValue` or :obj:`FilterResult`:

                * If `return_info` is `False`: return a list of the observed values which matches the filters, ordered by their mean score
                * If `return_info` is `True`: return a list of the observed values which matches the filters, ordered by their mean score plus additional diagnostic information
        """

        # score observed values against the filters
        all_scores = self.score(observed_values)

        # filter out observed values that must be discarded (observed values with score = -1)
        obs_vals, scores, i_obs_values = self.filter(observed_values, all_scores)

        # order remaining observed values by their mean score
        obs_vals, scores, i_obs_values = self.order(obs_vals, scores, i_obs_values)

        # return
        if return_info:
            # return ordered list of observed values and additional information
            return FilterResult(observed_values, self.filters, all_scores, obs_vals, scores, i_obs_values)
        else:
            # return ordered list of observed values
            return obs_vals

    def score(self, observed_values):
        """ Score observed values against the filters

        Args:
            observed_values (:obj:`list` of :obj:`data_model.ObservedValue`): list of experimental and/or computational observed values

        Returns:
            :obj:`list` of :obj:`float`: list of scores
        """
        n_obs = len(observed_values)
        n_filt = len(self.filters)
        scores = numpy.full((n_obs, n_filt, ), numpy.nan)
        for i_filter, filter in enumerate(self.filters):
            for i_obs, obs in enumerate(observed_values):
                scores[i_obs, i_filter] = filter.score(obs)

        return scores

    def filter(self, observed_values, scores):
        """ Filter out observed values that must be discarded (observed values with score = -1)

        Args:
            observed_values (:obj:`list` of :obj:`data_model.ObservedValue`): list of experimental and/or computational observed values
            scores (:obj:`list` of :obj:`float`): list of scores

        Returns:
            :obj:`tuple`: 

                * :obj:`list` of :obj:`data_model.ObservedValue`: list of acceptable observed values (observed values without scores = -1)
                * :obj:`list` of :obj:`float`: list of scores of the acceptable observed values
                * :obj:`list` of :obj:`int`: list of indices of the ordered observed values within the original list of observed values

        """
        ok_observations = numpy.extract(numpy.all(scores >= 0, 1).transpose(), observed_values).tolist()
        i_ok_observations = numpy.flatnonzero(numpy.all(scores >= 0, 1)).tolist()
        ok_scores = scores[i_ok_observations, :]

        return (ok_observations, ok_scores, i_ok_observations, )

    def order(self, observed_values, scores, i_observations=None):
        """ Order observed values by their mean score

        Args:
            observed_values (:obj:`list` of :obj:`data_model.ObservedValue`): list of observed values
            scores (:obj:`list` of :obj:`float`): list of scores
            i_observations (:obj:`list` of :obj:`int`, optional): list of indices within the original list of observed values

        Returns:
            :obj:`tuple`: 

                * :obj:`list` of :obj:`data_model.ObservedValue`: ordered list of observed values
                * :obj:`list` of :obj:`float`: list of scores of the ordered observed values
                * :obj:`list` of :obj:`int`: list of indices of the ordered observed values within the original list of observed values
        """

        if not i_observations:
            i_observations = range(len(observed_values))

        order = numpy.argsort(numpy.mean(scores, 1))[::-1]

        ordered_observed_values = [observed_values[i] for i in order]
        ordered_scores = scores[order, :]
        i_ordered_observed_values = [i_observations[i] for i in order]

        return (ordered_observed_values, ordered_scores, i_ordered_observed_values, )


class FilterResult(object):
    """ Represents the results of applying a list of filters to a dataset

    Attributes:
        observed_values (:obj:`list` of `data_model.ObservedValue`): input list of observed values
        filters (:obj:`list` of `Filter`): list of filters applied to observed values
        scores (:obj:`numpy.ndarray`): matrix of scores (rows: observed values in same order as in `observed_values`; columns: filters, in same orders as in `filters`)
        ordered_observed_values (:obj:`list` of `data_model.ObservedValue`): prioritized list of observed values
        ordered_scores (:obj:`numpy.ndarray`): matrix of scores (rows: observed values in same order as in `ordered_observed_values`; columns: filters, in same orders as in `filters`)
        ordered_observed_value_indices (:obj:`list` of :obj:`int`): indices of the ordered observed values in the input list of observed values
    """

    def __init__(self, observed_values, filters, scores, ordered_observed_values, ordered_scores, ordered_observed_value_indices):
        """
        Args:
            observed_values (:obj:`list` of `data_model.ObservedValue`): input list of observed values
            filters (:obj:`list` of `Filter`): list of filters applied to observed values
            scores (:obj:`numpy.ndarray`): matrix of scores (rows: observed values in same order as in `observed_values`; columns: filters, in same orders as in `filters`)
            ordered_observed_values (:obj:`list` of `data_model.ObservedValue`): prioritized list of observed values
            ordered_scores (:obj:`numpy.ndarray`): matrix of scores (rows: observed values in same order as in `ordered_observed_values`; columns: filters, in same orders as in `filters`)
            ordered_observed_value_indices (:obj:`list` of :obj:`int`): indices of the ordered observed values in the input list of observed values
        """
        self.observed_values = observed_values
        self.filters = filters
        self.scores = scores
        self.ordered_observed_values = ordered_observed_values
        self.ordered_scores = ordered_scores
        self.ordered_observed_value_indices = ordered_observed_value_indices


class Filter(object):
    """ Calculate a numeric score which indicates how well an observed value matches one or more criteria.
    Please see :obj:`FilterRunner` to see how these scores are used to filter and order observed values.

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

    def get_attribute_value(self, observed_value):
        """ Get the value of the attribute of observed value :obj:`observed_value`

        Args:
            observed_value (:obj:`data_model.ObservedValue`): observed value

        Returns:
            :obj:`object`: value of the attribute of the observed value
        """
        val = observed_value
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

    def score(self, observed_value):
        """ Calculate a scaled numeric score betwen 0 and 1 which indicates how well the observed value matches 
        one or more criteria. Please see :obj:`FilterRunner` to see how these scores are used to filter and 
        order observed values.

        Args:
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`float`: score which indicates how well the observed value matches the criteria
        """
        return self.transform_attribute_value(self.get_attribute_value(observed_value))


class OptionsFilter(Filter):
    """ Filters out observed values whose attributes have values that are not in a list of acceptable options.

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

    def score(self, observed_value):
        """ Calculate a numeric score which indicates how well the observed value matches 
        one or more criteria. Please see :obj:`FilterRunner` to see how these scores are used to filter and 
        order observed values.

        Args:
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`float`: score which indicates how well the observed value matches the criteria
        """
        val = self.transform_attribute_value(self.get_attribute_value(observed_value))
        if val in self.options:
            return 1
        return -1


class RangeFilter(Filter):
    """ Filters out observed values whose attributes have values that fall outside a specified range.

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

    def get_raw_score(self, observed_value):
        """ Calculate a numeric score which indicates how well the observed value matches one or more criteria
        Please see :obj:`FilterRunner` to see how these scores are used to filter and order observed values.

        Args:
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`float`: score which indicates how well the observed value matches the criteria
        """
        return self.get_attribute_value(observed_value)

    def score(self, observed_value):
        """ Calculate a scaled numeric score betwen 0 and 1 which indicates how well the observed value matches 
        one or more criteria. Please see :obj:`FilterRunner` to see how these scores are used to filter and 
        order observed values.

        Args:
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`float`: score which indicates how well the observed value matches the criteria
        """
        val = self.transform_attribute_value(self.get_attribute_value(observed_value))

        if not numpy.isnan(self.min) and (numpy.isnan(val) or val < self.min):
            return -1
        if not numpy.isnan(self.max) and (numpy.isnan(val) or val > self.max):
            return -1
        return 1


class NormalFilter(Filter):
    """ Prioritizes observed values whose attributes have values that are closed to `mean`.

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

    def score(self, observed_value):
        """ Calculate a numeric score which indicates how well the attribute of the observed value matches the specified
        normal distribution (mean, std).

        Args:
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`float`: score which indicates how well the observed value matches the normal distribution (mean, std)
        """
        val = self.transform_attribute_value(self.get_attribute_value(observed_value))

        return 1 - 2 * abs(scipy.stats.norm.cdf(val, loc=self.mean, scale=self.std) - 0.5)


class ExponentialFilter(Filter):
    """ Prioritizes observed values based on an exponential scale

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

    def score(self, observed_value):
        """ Calculate a numeric score which indicates how well the attribute of the observed value matches the specified
        distribution (center, scale).

        Args:
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`float`: score which indicates how well the observed value matches the distribution (center, scale)
        """
        val = self.transform_attribute_value(self.get_attribute_value(observed_value))

        return math.exp(-(val - self.center) / self.scale)


class TaxonomicDistanceFilter(Filter):
    """ Prioritizes observed values that are from taxonomically close taxa

    Attributes:
        taxon (:obj:`str`): name of the taxon to find data for
    """

    def __init__(self, taxon, max=None, scale=None):
        """
        Args:
            taxon (:obj:`str`): name of the taxon to find data for            
            max (:obj:`float`, optional): maximum distance to the latest common ancestor with the observed taxon
            scale (:obj:`float`, optional): The scale of the distribution. This determines how quickly the score falls to zero away from the center.
        """

        if max is None or numpy.isnan(max):
            taxon_obj = taxonomy_util.Taxon(name=taxon)
            max = taxon_obj.get_max_distance_to_common_ancestor() - 2

        if scale is None or numpy.isnan(scale):
            taxon_obj = taxonomy_util.Taxon(name=taxon)
            scale = (taxon_obj.get_max_distance_to_common_ancestor() - 2) / 5.

        super(TaxonomicDistanceFilter, self).__init__(('observation', 'genetics', 'taxon', ))

        self.taxon = taxon
        self.max = max
        self.scale = scale

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

    def score(self, observed_value):
        """ Calculate a scaled numeric score betwen 0 and 1 which indicates how well the observed value matches 
        one or more criteria. Please see :obj:`FilterRunner` to see how these scores are used to filter and 
        order observed values.

        Args:
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`float`: score which indicates how well the observed value matches the criteria
        """
        val = self.transform_attribute_value(self.get_attribute_value(observed_value))

        if not numpy.isnan(self.max) and (numpy.isnan(val) or val > self.max):
            return -1

        return math.exp(-val / self.scale)


class WildtypeFilter(OptionsFilter):
    """ Filter out observed values which were observed for taxa with genetic perturbations """

    def __init__(self):
        super(WildtypeFilter, self).__init__(('observation', 'genetics', 'variation', ), [''])


class SpecieMolecularSimilarityFilter(Filter):
    """ Score similarity with observed molecules

    Attributes:
        structure (:obj:`str`): 
    """

    def __init__(self, structure):
        """
        Args:
            structure (:obj:`str`): structure
        """
        super(SpecieMolecularSimilarityFilter, self).__init__(('observable', 'structure', ), )
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
    pass  # todo: implement for small molecules, sequences


class ReactionSimilarityFilter(Filter):
    """ Prioritize reactions by their chemical similarity, as judged by (a) having the same participants and (b)
    belonging to the same EC class (or superclass).

    * 1: Reactions have the same participants
    * <0, 1>: Reactions belong to the same EC class (or superclass), but don't have different participants
    * Score=-1: Reactions belong to different EC classes

    Attributes:
        reaction (:obj:`data_model.Reaction`): reaction
        min_ec_level (:obj:`int`): minimum EC level that must be common to the observed and target reaction 
        scale (:obj:`float`): How to exponentially scale of the scores. This determines how quickly the score 
            falls to zero.
    """

    def __init__(self, reaction, min_ec_level=3, scale=2./5.):
        """
        Args:
            reaction (:obj:`data_model.Reaction`): reaction
            min_ec_level (:obj:`int`, optional): minimum EC level that must be common to the observed and target reaction 
            scale (:obj:`float`, optional): How to exponentially scale of the scores. This determines how quickly the score 
                falls to zero.
        """
        super(ReactionSimilarityFilter, self).__init__(('observable',))
        self.reaction = reaction
        self.min_ec_level = min_ec_level
        self.scale = scale

    def transform_attribute_value(self, reaction):
        """ Transform an attribute value

        Args:
            reaction (:obj:`data_model.Reaction`): observed reaction

        Returns:
            :obj:`int`: transformed value
        """
        # todo: use weights associated with predicted EC numbers

        # check participants are the same
        def get_formula_connectivities(participants):
            vals = set()
            for part in participants:
                if not part.specie.structure:
                    return None
                inchi = molecule_util.Molecule(structure=part.specie.structure).to_inchi()
                val = molecule_util.InchiMolecule(inchi).get_formula_and_connectivity()
                if val not in [None, '', 'H2O']:
                    vals.add(val)
            return vals

        reactants = get_formula_connectivities(reaction.get_reactants())
        products = get_formula_connectivities(reaction.get_products())
        self_reactants = get_formula_connectivities(self.reaction.get_reactants())
        self_products = get_formula_connectivities(self.reaction.get_products())
        if reactants is not None and self_reactants is not None and reactants == self_reactants and \
                products is not None and self_products is not None and products == self_products:
            return 5

        # check membership to same EC class
        ecs = reaction.get_ec_numbers()
        self_ecs = self.reaction.get_ec_numbers()

        for level in range(4, 0, -1):
            part_ecs = []
            for ec in ecs:
                if ec.id.count('.') >= level - 1:
                    part_ec, _, _ = string.partition_nth(ec.id, '.', level)
                    if part_ec[-1] != '.':
                        part_ecs.append(part_ec)

            self_part_ecs = []
            for self_ec in self_ecs:
                if self_ec.id.count('.') >= level - 1:
                    self_part_ec, _, _ = string.partition_nth(self_ec.id, '.', level)
                    if self_part_ec[-1] != '.':
                        self_part_ecs.append(self_part_ec)

            if part_ecs and set(part_ecs).intersection(self_part_ecs):
                return level

        return 0

    def score(self, observed_value):
        """ Calculate a scaled numeric score betwen 0 and 1 which indicates how well the observed value matches 
        one or more criteria. Please see :obj:`FilterRunner` to see how these scores are used to filter and 
        order observed values.

        Args:
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`float`: score which indicates how well the observed_value matches the criteria
        """
        val = self.transform_attribute_value(self.get_attribute_value(observed_value))

        if val >= self.min_ec_level:
            return math.exp(-(5. - val) / self.scale)
        return -1


class TemperatureRangeFilter(RangeFilter):
    """ Filters out observed values with temperatures that fall outside a specified range. """

    def __init__(self, min=float('nan'), max=float('nan')):
        """
        Args:
            min (:obj:`float`, optional): minimum value
            max (:obj:`float`, optional): maximum value
        """
        super(TemperatureRangeFilter, self).__init__(('observation', 'environment', 'temperature', ), min=min, max=max)


class TemperatureNormalFilter(NormalFilter):
    """ Prioritizes observed values with temperatures that are close to `mean`. """

    def __init__(self, mean, std):
        """
        Args:
            mean (:obj:`float`): The mean of the distribution. This indicates the value at which the score will be 1.
            std (:obj:`float`): The standard deviation of the distribution. This determines how quickly the score falls to zero away from the mean.
        """
        super(TemperatureNormalFilter, self).__init__(('observation', 'environment', 'temperature', ), mean, std)


class PhRangeFilter(RangeFilter):
    """ Filters out observed values with pHs that fall outside a specified range. """

    def __init__(self, min=float('nan'), max=float('nan')):
        """
        Args:
            min (:obj:`float`, optional): minimum value
            max (:obj:`float`, optional): maximum value
        """
        super(PhRangeFilter, self).__init__(('observation', 'environment', 'ph', ), min=min, max=max)


class PhNormalFilter(NormalFilter):
    """ Prioritizes observed values with pHs that are close to `mean`. """

    def __init__(self, mean, std):
        """
        Args:
            mean (:obj:`float`): The mean of the distribution. This indicates the value at which the score will be 1.
            std (:obj:`float`): The standard deviation of the distribution. This determines how quickly the score 
                falls to zero away from the mean.
        """
        super(PhNormalFilter, self).__init__(('observation', 'environment', 'ph', ), mean, std)
