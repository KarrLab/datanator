"""
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-05-16
:Copyright: 2017, Karr Lab
:License: MIT
"""

from datetime import datetime
from kinetic_datanator.core import data_model
from kinetic_datanator.util import molecule_util
from kinetic_datanator.util import taxonomy_util
import abc
import getpass
import itertools
import kinetic_datanator.core.data_source
import Levenshtein
import math
import numpy
import scipy.stats
import six
import wc_utils.util.stats
import wc_utils.util.string


class DataQueryGenerator(six.with_metaclass(abc.ABCMeta, object)):
    """ Represents a query of a data source

    #. Find observed values for the exact or similar model components
    #. Filter out observed values from disimilar genetic and environmental conditions and
       rank the remaing observed values by their similarity to the desired genetic and environmental
       conditions

        * Taxonomy
        * Genetic variation (wildtype/mutant)
        * Temperature
        * pH

    #. Calculate a statistical representation of the relevant observed values

    Attributes:
        filters (:obj:`list` of :obj:`Filter`): list of filters
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
            include_variants (:obj:`bool`, optional): if :obj:`True`, also include observed values from mutant taxa
            temperature (:obj:`float`, optional): desired temperature to search for
            temperature_std (:obj:`float`, optional): how much to penalize observed values from other temperatures
            ph (:obj:`float`, optional): desired pH to search for
            ph_std (:obj:`float`, optional): how much to penalize observed values from other pHs
        """
        # todo (enhancement): filter media

        self.filters = filters = []

        if taxon:
            filters.append(TaxonomicDistanceFilter(taxon, max=max_taxon_dist, scale=taxon_dist_scale))

        if not include_variants:
            filters.append(WildtypeFilter())

        if not math.isnan(temperature) and not math.isnan(temperature_std):
            filters.append(TemperatureNormalFilter(temperature, temperature_std))

        if not math.isnan(ph) and not math.isnan(ph_std):
            filters.append(PhNormalFilter(ph, ph_std))

    def run(self, component):
        """

        1. Find observed values for the exact or similar model components and genetic and environmental conditions
        2. Rank the results by their similarity to the model component and the genetic and environmental conditions
        3. Calculate a consensus statistical representation of the relevant observed values

        Args:
            component (:obj:`data_model.EntityInteractionOrProperty`): model component to find data for

        Returns:
            :obj:`list` of :obj:`data_model.Consensus`: statistical consensus of the relevant observed values of
                :obj:`component` and the observed values it was based on
        """
        observed_result = self.get_observed_result(component)
        filter_result = self.filter_observed_results(component, observed_result)
        return filter_result
        # return self.get_consensus(component, filter_result)

    @abc.abstractmethod
    def get_observed_result(self, component):
        """ Find the observed result relevant to :obj:`component`

        Args:
            component (:obj:`model.PhysicalEntity`) or (:obj:`model.PhysicalProperty`): model component to find data for

        Returns:
            :obj:`list` of :obj:`data_model.ObservedValue`: list of relevant observed values
            :obj:`list` of :obj:`data_model.ObservedValue`: list of relevant observed values
        """
        pass

    def filter_observed_results(self, component, observed_results):
        """ Filter out observed values from dissimilar genetic and environmental conditions and
        order the remaining observed values by their similarity to specified genetic and
        environmental conditions.

        Args:
            component (:obj:`data_model.EntityInteractionOrProperty`): model component to find data for
            observed_results (:obj:`list` of :obj:`data_model.ObservedValue`): list of observed values

        Returns:
            :obj:`FilterResult`: filter result
        """
        return FilterRunner(self.filters) \
            .run(component, observed_results, return_info=True)

    def get_consensus(self, component, filter_result):
        """ Calculate a consensus statistical representation of the one or more observed values

        Args:
            component (:obj:`data_model.EntityInteractionOrProperty`): model component to find data for
            filter_result (:obj:`FilterResult`): filter result

        Returns:
            :obj:`list` of :obj:`data_model.Consensus`: statistical consensus of the relevant observed values of
                :obj:`component` and the observed values it was based on
        """
        # group observed values by their subcomponents, attributes
        return consensus

    def metadata_dump(self, component):
        """ Calculate a consensus statistical representation of the one or more observed values

        Args:
            component (:obj:`models.PhysicalEntity` or :obj:`models.PhysicalProperty`): model component dump data for

        Returns:
            :obj:`list` of :obj:`data_model.ObservedResultMetadata`: data model metadata object
        """
        genetics = None
        environment=None
        cross_references= []
        method=None
        synonym= []
        meta = component._metadata


        taxon = meta.taxon[0].name if meta.taxon else None
        variation = meta.cell_line[0].name if meta.cell_line else None
        genetics = data_model.Genetics(taxon=taxon, variation=variation)

        temperature = meta.conditions[0].temperature if meta.conditions else None
        ph = meta.conditions[0].ph if meta.conditions else None
        media = meta.conditions[0].media if meta.conditions else None
        growth_status = meta.conditions[0].growth_status if meta.conditions else None
        growth_system = meta.conditions[0].growth_system if meta.conditions else None
        environment = data_model.Environment(temperature=temperature, ph=ph, media=media, growth_status=growth_status, growth_system=growth_system)

        name = meta.method[0].name if meta.method else None
        description =  meta.method[0].comments if meta.method else None
        performer = meta.method[0].performer if meta.method else None
        hardware =  meta.method[0].hardware if meta.method else None
        software =  meta.method[0].software if meta.method else None
        method = data_model.Method(name=name, description=description, performer=performer, hardware=hardware, software=software)


        if meta.resource:
            for item in meta.resource:
                cross_references.append(data_model.Resource(namespace=item.namespace, id=item._id))
        else:
            cross_references.append(data_model.Resource(namespace=None, id=None))


        if meta.synonym:
            for item in meta.synonym:
                synonym.append(data_model.Synonym(name=item.name))
        else:
            synonym.append(data_model.Synonym(name=None))

        metadata_result = data_model.ObservedResultMetadata(genetics = genetics, environment=environment, cross_references=cross_references, method=method, synonym=synonym)

        return metadata_result


class CachedDataSourceQueryGenerator(DataQueryGenerator):
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
            include_variants (:obj:`bool`, optional): if :obj:`True`, also include observed values from mutant taxa
            temperature (:obj:`float`, optional): desired temperature to search for
            temperature_std (:obj:`float`, optional): how much to penalize observed values from other temperatures
            ph (:obj:`float`, optional): desired pH to search for
            ph_std (:obj:`float`, optional): how much to penalize observed values from other pHs
            data_source (:obj:`kinetic_datanator.core.data_source.CachedDataSource`, optional): cached data source
        """
        super(CachedDataSourceQueryGenerator, self).__init__(
            taxon=taxon, max_taxon_dist=max_taxon_dist, taxon_dist_scale=taxon_dist_scale, include_variants=include_variants,
            temperature=temperature, temperature_std=temperature_std,
            ph=ph, ph_std=ph_std)

        self.data_source = data_source


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

    def run(self, target_component, observed_results, return_info=False):
        """ Filter and order a list of observed values according to a list of filters. Optionally, return additional
        information about the filtering including the scores of the observed values and the indices of the prioritized
        observed values in the input list of observed values.

        1. Calculate the score of each observed value for each filter

            * Scores equal to -1, indicate that the observed value should be discarded
            * Scores between 0 and 1, indicate how much the observed value should be prioritized

        2. Discard any observed value which has at least one score equal to -1
        3. Order the observed values by their mean score

        Args:
            target_component (:obj:`data_model.EntityInteractionOrProperty`): interaction, species, or property to find data about
            observed_results (:obj:`list` of :obj:`data_model.ObservedValue`): list of experimental and/or computational observed values
            return_info (:obj:`bool`, optional): if `True`, also return the scores and indices of the ordered observed values in the input list

        Returns:
            :obj:`list` of :obj:`data_model.ObservedValue` or :obj:`FilterResult`:

                * If `return_info` is `False`: return a list of the observed values which matches the filters, ordered by their mean score
                * If `return_info` is `True`: return a list of the observed values which matches the filters, ordered by their mean score plus additional diagnostic information
        """

        # score observed values against the filters
        all_scores = self.score(target_component, observed_results)

        # filter out observed values that must be discarded (observed values with score = -1)
        obs_vals, scores, i_obs_values = self.filter(observed_results, all_scores)

        # order remaining observed values by their mean score
        obs_vals, scores, i_obs_values = self.order(obs_vals, scores, i_obs_values)

        # return
        if return_info:
            # return ordered list of observed values and additional information
            return FilterResult(obs_vals, scores, i_obs_values, observed_results, all_scores)
        else:
            # return ordered list of observed values
            return obs_vals

    def score(self, target_component, observed_results):
        """ Score observed values against the filters

        Args:
            target_component (:obj:`data_model.EntityInteractionOrProperty`): interaction, species, or property to find data about
            observed_results (:obj:`list` of :obj:`data_model.ObservedValue`): list of experimental and/or computational observed values

        Returns:
            :obj:`list` of :obj:`float`: list of scores
        """
        n_obs = len(observed_results)
        n_filt = len(self.filters)
        scores = numpy.full((n_obs, n_filt, ), numpy.nan)
        for i_filter, filter in enumerate(self.filters):
            for i_obs, obs in enumerate(observed_results):
                scores[i_obs, i_filter] = filter.score(target_component, obs)

        return scores

    def filter(self, observed_results, scores):
        """ Filter out observed values that must be discarded (observed values with score = -1)

        Args:
            observed_results (:obj:`list` of :obj:`data_model.ObservedValue`): list of experimental and/or computational observed values
            scores (:obj:`list` of :obj:`float`): list of scores

        Returns:
            :obj:`tuple`:

                * :obj:`list` of :obj:`data_model.ObservedValue`: list of acceptable observed values (observed values without scores = -1)
                * :obj:`list` of :obj:`float`: list of scores of the acceptable observed values
                * :obj:`list` of :obj:`int`: list of indices of the ordered observed values within the original list of observed values

        """
        ok_observed_results = numpy.extract(numpy.all(scores >= 0, 1).transpose(), observed_results).tolist()
        i_ok_observed_results = numpy.flatnonzero(numpy.all(scores >= 0, 1)).tolist()
        ok_scores = scores[i_ok_observed_results, :]
        return (ok_observed_results, ok_scores, i_ok_observed_results, )

    def order(self, observed_results, scores, i_observations=None):
        """ Order observed values by their mean score

        Args:
            observed_results (:obj:`list` of :obj:`data_model.ObservedValue`): list of observed values
            scores (:obj:`list` of :obj:`float`): list of scores
            i_observations (:obj:`list` of :obj:`int`, optional): list of indices within the original list of observed values

        Returns:
            :obj:`tuple`:

                * :obj:`list` of :obj:`data_model.ObservedValue`: ordered list of observed values
                * :obj:`list` of :obj:`float`: list of scores of the ordered observed values
                * :obj:`list` of :obj:`int`: list of indices of the ordered observed values within the original list of observed values
        """

        if not i_observations:
            i_observations = range(len(observed_results))

        order = numpy.argsort(numpy.mean(scores, 1))[::-1]

        ordered_observed_results = [observed_results[i] for i in order]
        ordered_scores = scores[order, :]
        i_ordered_observed_results = [i_observations[i] for i in order]

        return (ordered_observed_results, ordered_scores, i_ordered_observed_results, )


class FilterResult(object):
    """ Represents the results of applying a list of filters to a dataset

    Attributes:
        observed_results (:obj:`list` of `data_model.ObservedValue`): prioritized list of observed values
        scores (:obj:`numpy.ndarray`): matrix of scores (rows: observed values in same order as in `ordered_observed_results`; columns: filters, in same orders as in `filters`)
        observed_value_indices (:obj:`list` of :obj:`int`): indices of the ordered observed values in the input list of observed values
        all_observed_results (:obj:`list` of `data_model.ObservedValue`): input list of observed values
        all_scores (:obj:`numpy.ndarray`): matrix of scores (rows: observed values in same order as in `observed_results`; columns: filters, in same orders as in `filters`)
    """

    def __init__(self, observed_results, scores, observed_value_indices, all_observed_results, all_scores):
        """
        Args:
            observed_results (:obj:`list` of `data_model.ObservedValue`): prioritized list of observed values
            scores (:obj:`numpy.ndarray`): matrix of scores (rows: observed values in same order as in `ordered_observed_results`; columns: filters, in same orders as in `filters`)
            observed_value_indices (:obj:`list` of :obj:`int`): indices of the ordered observed values in the input list of observed values
            all_observed_results (:obj:`list` of `data_model.ObservedValue`): input list of observed values
            all_scores (:obj:`numpy.ndarray`): matrix of scores (rows: observed values in same order as in `observed_results`; columns: filters, in same orders as in `filters`)
        """
        self.observed_results = observed_results
        self.scores = scores
        self.observed_value_indices = observed_value_indices
        self.all_observed_results = all_observed_results
        self.all_scores = all_scores


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

    def get_attribute_of_observed_value(self, observed_value):
        """ Get the value of the attribute of observed value :obj:`observed_value`

        Args:
            observed_value (:obj:`data_model.ObservedValue`): observed value

        Returns:
            :obj:`object`: value of the attribute of the observed value
        """
        val = observed_value
        for attr in self.attribute:
            if hasattr(val, attr):
                val = getattr(val, attr)
            else:
                return None
        return val

    def compare_observed_value_with_target_component(self, target_component, observed_value):
        """ Compare the observed biological component with the target component

        Args:
            target_component (:obj:`data_model.EntityInteractionOrProperty`): interaction, species, or property to find data about
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`object`: transformed value
        """
        return self.get_attribute_of_observed_value(observed_value)

    def score(self, target_component, observed_value):
        """ Calculate a scaled numeric score betwen 0 and 1 which indicates how well the observed value matches
        one or more criteria. Please see :obj:`FilterRunner` to see how these scores are used to filter and
        order observed values.

        Args:
            target_component (:obj:`data_model.EntityInteractionOrProperty`): interaction, species, or property to find data about
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`float`: score which indicates how well the observed value matches the criteria
        """
        val = self.compare_observed_value_with_target_component(target_component, observed_value)
        if val is None:
            return 0.5
        return val


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

    def score(self, target_component, observed_value):
        """ Calculate a numeric score which indicates how well the observed value matches
        one or more criteria. Please see :obj:`FilterRunner` to see how these scores are used to filter and
        order observed values.

        Args:
            target_component (:obj:`data_model.EntityInteractionOrProperty`): interaction, species, or property to find data about
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`float`: score which indicates how well the observed value matches the criteria
        """
        val = self.compare_observed_value_with_target_component(target_component, observed_value)
        if val is None:
            return 0.5
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

    def score(self, target_component, observed_value):
        """ Calculate a scaled numeric score betwen 0 and 1 which indicates how well the observed value matches
        one or more criteria. Please see :obj:`FilterRunner` to see how these scores are used to filter and
        order observed values.

        Args:
            target_component (:obj:`data_model.EntityInteractionOrProperty`): interaction, species, or property to find data about
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`float`: score which indicates how well the observed value matches the criteria
        """
        val = self.compare_observed_value_with_target_component(target_component, observed_value)
        if val is None:
            return 0.5
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

    def score(self, target_component, observed_value):
        """ Calculate a numeric score which indicates how well the attribute of the observed value matches the specified
        normal distribution (mean, std).

        Args:
            target_component (:obj:`data_model.EntityInteractionOrProperty`): interaction, species, or property to find data about
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`float`: score which indicates how well the observed value matches the normal distribution (mean, std)
        """
        val = self.compare_observed_value_with_target_component(target_component, observed_value)
        if val is None:
            return 0.5
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

    def score(self, target_component, observed_value):
        """ Calculate a numeric score which indicates how well the attribute of the observed value matches the specified
        distribution (center, scale).

        Args:
            target_component (:obj:`data_model.EntityInteractionOrProperty`): interaction, species, or property to find data about
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`float`: score which indicates how well the observed value matches the distribution (center, scale)
        """
        val = self.compare_observed_value_with_target_component(target_component, observed_value)
        if val is None:
            return 0.5
        return math.exp(-(val - self.center) / self.scale)


class SpecieSimilarityFilter(Filter):
    """ Proritize observed species based on their similarity to target species

    Attributes:
        min_similarity (:obj:`float`): minimum acceptable similarity
    """

    def __init__(self, filters, min_similarity=0.5):
        """
        Args:
            filters (:obj:`Filter` or :obj:`list` of :obj:`Filter`): filter or list of filters
            min_similarity (:obj:`float`, optional): minimum acceptable similarity
        """
        super(SpecieSimilarityFilter, self).__init__(filters)
        self.min_similarity = min_similarity

    def score(self, target_component, observed_value):
        """ Calculate a scaled numeric score betwen 0 and 1 which indicates how well the sequence of the
        observed species match that of the target species.

        Args:
            target_component (:obj:`data_model.Specie`): species to find data about
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`float`: similarity between the observed and target structures
        """
        val = self.compare_observed_value_with_target_component(target_component, observed_value)
        if val < self.min_similarity:
            return -1
        return val


class SpecieStructuralSimilarityFilter(SpecieSimilarityFilter):
    """ Proritize observed species based on their structural similarity with target species """

    def __init__(self, min_similarity=0.5):
        """
        Args:
            min_similarity (:obj:`float`, optional): minimum acceptable similarity
        """
        super(SpecieStructuralSimilarityFilter, self).__init__(('observable', 'specie', 'structure', ),
                                                               min_similarity=min_similarity)

    def compare_observed_value_with_target_component(self, target_component, observed_value):
        """ Compare the observed biological component with the target component

        Args:
            target_component (:obj:`data_model.Specie`): species to find data about
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`float`: similarity between the observed and target structures
        """
        target_structure = target_component.structure
        observed_structure = self.get_attribute_of_observed_value(observed_value)
        target_mol = molecule_util.Molecule(structure=target_structure)
        observed_mol = molecule_util.Molecule(structure=observed_structure)
        return target_mol.get_similarity(observed_mol)


class SpecieSequenceSimilarityFilter(SpecieSimilarityFilter):
    """ Proritize observed species based on the Levenshtein distance of their sequences to that of target species """

    def __init__(self, min_similarity=0.5):
        """
        Args:
            min_similarity (:obj:`float`, optional): minimum acceptable similarity
        """
        super(SpecieSequenceSimilarityFilter, self).__init__(('observable', 'specie', 'sequence', ),
                                                             min_similarity=min_similarity)

    def compare_observed_value_with_target_component(self, target_component, observed_value):
        """ Compare the observed biological component with the target component

        Args:
            target_component (:obj:`data_model.Specie`): species to find data about
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`float`: similarity between the observed and target sequences
        """
        target_sequence = target_component.sequence
        observed_sequence = self.get_attribute_of_observed_value(observed_value)

        dist = Levenshtein.distance(target_sequence, observed_sequence)
        return 1. - float(dist) / float(max(len(target_sequence), len(observed_sequence)))


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

    def __init__(self, min_ec_level=3, scale=2./5.):
        """
        Args:
            min_ec_level (:obj:`int`, optional): minimum EC level that must be common to the observed and target reaction
            scale (:obj:`float`, optional): How to exponentially scale of the scores. This determines how quickly the score
                falls to zero.
        """
        super(ReactionSimilarityFilter, self).__init__(('observable', 'interaction'))
        self.min_ec_level = min_ec_level
        self.scale = scale

    def compare_observed_value_with_target_component(self, target_component, observed_value):
        """ Compare the observed biological component with the target component

        Args:
            target_component (:obj:`data_model.Reaction`): reaction to find data about
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`int`: transformed value
        """
        # todo (enhancement): use weights associated with predicted EC numbers

        target_reaction = target_component
        observed_reaction = self.get_attribute_of_observed_value(observed_value)

        def get_participant_species(participants):
            species = []
            for part in participants:
                if not part.specie or not part.specie.structure:
                    return None
                species.append(part.specie)
            return species

        target_reactants = get_participant_species(target_reaction.get_reactants())
        target_products = get_participant_species(target_reaction.get_products())
        observed_reactants = get_participant_species(observed_reaction.get_reactants())
        observed_products = get_participant_species(observed_reaction.get_products())

        # check participants are the same
        def get_formula_connectivities(species):
            if species is None:
                return None

            vals = set()
            for specie in species:
                val = specie.to_inchi(only_formula_and_connectivity=True)
                if val not in [None, '', 'H2O']:
                    vals.add(val)
            return vals

        target_reactant_connectivities = get_formula_connectivities(target_reactants)
        target_product_connectivities = get_formula_connectivities(target_products)
        observed_reactant_connectivities = get_formula_connectivities(observed_reactants)
        observed_product_connectivities = get_formula_connectivities(observed_products)
        if \
                target_reactant_connectivities is not None and \
                target_product_connectivities is not None and \
                observed_reactant_connectivities is not None and \
                observed_product_connectivities is not None and \
                observed_reactant_connectivities == target_reactant_connectivities and \
                observed_product_connectivities == target_product_connectivities:
            return 5

        # check numbers of reactants and products are the same
        if observed_reactants and observed_products and \
                (len(observed_reactant_connectivities) != len(target_reactant_connectivities) or
                    len(observed_product_connectivities) != len(target_product_connectivities)):
            return -1

        # check directions are the same
        def get_side_similarity(target_connectivities, observed_connectivities):
            similarities = numpy.full((len(target_connectivities), len(observed_connectivities)), numpy.nan)
            for t, t_conn in enumerate(target_connectivities):
                for o, o_conn in enumerate(observed_connectivities):
                    similarities[t, o] = molecule_util.Molecule(structure='InChI=1S/' + t_conn) \
                        .get_similarity(molecule_util.Molecule(structure='InChI=1S/' + o_conn))

            # greedy search
            max_vals = []
            for i in range(len(target_connectivities)):
                i_max = numpy.argmax(similarities[i, :])
                max_vals.append(similarities[i, i_max])
                similarities = numpy.delete(similarities, i_max, 1)
            return numpy.mean(max_vals)

        if observed_reactants and observed_products and len(target_reactant_connectivities) == len(target_product_connectivities):
            if \
                    get_side_similarity(target_reactant_connectivities, observed_reactant_connectivities) < \
                    get_side_similarity(target_product_connectivities, observed_reactant_connectivities) or \
                    get_side_similarity(target_product_connectivities, observed_product_connectivities) < \
                    get_side_similarity(target_reactant_connectivities, observed_product_connectivities):
                return -1

        # check membership to same EC class
        target_ecs = target_reaction.get_ec_numbers()
        observed_ecs = observed_reaction.get_ec_numbers()

        for level in range(4, 0, -1):
            part_target_ecs = []
            for target_ec in target_ecs:
                if target_ec.id.count('.') >= level - 1:
                    part_target_ec, _, _ = wc_utils.util.string.partition_nth(target_ec.id, '.', level)
                    if part_target_ec[-1] != '.':
                        part_target_ecs.append(part_target_ec)

            part_observed_ecs = []
            for observed_ec in observed_ecs:
                if observed_ec.id.count('.') >= level - 1:
                    part_observed_ec, _, _ = wc_utils.util.string.partition_nth(observed_ec.id, '.', level)
                    if part_observed_ec[-1] != '.':
                        part_observed_ecs.append(part_observed_ec)

            if part_observed_ecs and set(part_observed_ecs).intersection(part_target_ecs):
                return level

        return 0

    def score(self, target_component, observed_value):
        """ Calculate a scaled numeric score betwen 0 and 1 which indicates how well the observed value matches
        one or more criteria. Please see :obj:`FilterRunner` to see how these scores are used to filter and
        order observed values.

        Args:
            target_component (:obj:`data_model.EntityInteractionOrProperty`): interaction, species, or property to find data about
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`float`: score which indicates how well the observed_value matches the criteria
        """
        val = self.compare_observed_value_with_target_component(target_component, observed_value)
        if val is None:
            return 0.5
        if val >= self.min_ec_level:
            return math.exp(-(5. - val) / self.scale)
        return -1


class ReactionParticipantFilter(Filter):
    """
    Attributes:
        min_similarity (:obj:`float`): mini
    """

    def __init__(self, min_similarity=0.5):
        """
        Args:
            min_similarity (:obj:`float`, optional): mini
        """
        super(ReactionParticipantFilter, self).__init__(('observable', 'specie', 'structure'))
        self.min_similarity = min_similarity

    def compare_observed_value_with_target_component(self, target_component, observed_value):
        """ Compare the observed biological component with the target component

        Args:
            target_component (:obj:`data_model.EntityInteractionOrProperty`): interaction, species, or property to find data about
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`int`: transformed value
        """
        observed_reaction = observed_value.observable.interaction
        observed_specie = observed_value.observable.specie

        # if observed property is not associated with a species, return 1
        if not observed_specie:
            return 1

        # if the observed species has no structure, return None
        observed_structure = observed_specie.structure
        if not observed_structure:
            return None

        # if the observed structure matches one of the target structures, return 1
        if next((True for part in observed_reaction.get_reactants() if part.specie is observed_specie), None):
            target_parts = target_component.get_reactants()
        elif next((True for part in observed_reaction.get_products() if part.specie is observed_specie), None):
            target_parts = target_component.get_products()
        elif next((True for part in observed_reaction.get_modifiers() if part.specie is observed_specie), None):
            target_parts = target_component.get_modifiers()
        else:
            return None

        if not target_parts:
            return None

        target_species = [part.specie for part in target_parts]
        target_structures = [specie.structure for specie in target_species]
        if observed_structure in target_structures:
            return 1

        # if the observed structure is similar to one of the target structures, return 1
        target_formula_connectivities = [specie.to_inchi(only_formula_and_connectivity=True) for specie in target_species]
        observed_formula_connectivities = observed_specie.to_inchi(only_formula_and_connectivity=True)
        if observed_formula_connectivities in target_formula_connectivities:
            return 1

        # return maximal similarity to target species
        return max(target_specie.get_similarity(observed_specie) for target_specie in target_species)

    def score(self, target_component, observed_value):
        """ Calculate a scaled numeric score betwen 0 and 1 which indicates how well the observed value matches
        one or more criteria. Please see :obj:`FilterRunner` to see how these scores are used to filter and
        order observed values.

        Args:
            target_component (:obj:`data_model.EntityInteractionOrProperty`): interaction, species, or property to find data about
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`float`: score which indicates how well the observed_value matches the criteria
        """
        val = self.compare_observed_value_with_target_component(target_component, observed_value)
        if val is None:
            return -1
        if val < self.min_similarity:
            return -1
        return val


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

        super(TaxonomicDistanceFilter, self).__init__(('metadata', 'genetics', 'taxon', ))

        self.taxon = taxon
        self.max = max
        self.scale = scale

    def compare_observed_value_with_target_component(self, target_component, observed_value):
        """ Compare the observed biological component with the target component

        Args:
            target_component (:obj:`data_model.EntityInteractionOrProperty`): target component
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`int`: distance to latest common ancestor with the observed taxon
        """
        observed_taxon = self.get_attribute_of_observed_value(observed_value)
        if not observed_taxon:
            return None

        target_taxon_obj = taxonomy_util.Taxon(name=self.taxon)
        obsserved_taxon_obj = taxonomy_util.Taxon(name=observed_taxon)

        return target_taxon_obj.get_distance_to_common_ancestor(obsserved_taxon_obj)

    def score(self, target_component, observed_value):
        """ Calculate a scaled numeric score betwen 0 and 1 which indicates how well the observed value matches
        one or more criteria. Please see :obj:`FilterRunner` to see how these scores are used to filter and
        order observed values.

        Args:
            target_component (:obj:`data_model.EntityInteractionOrProperty`): interaction, species, or property to find data about
            observed_value (:obj:`data_model.ObservedValue`): experimentally or computationally observed value

        Returns:
            :obj:`float`: score which indicates how well the observed value matches the criteria
        """
        val = self.compare_observed_value_with_target_component(target_component, observed_value)
        if val is None:
            return 0.5

        if not numpy.isnan(self.max) and (numpy.isnan(val) or val > self.max):
            return -1

        return math.exp(-val / self.scale)


class WildtypeFilter(OptionsFilter):
    """ Filter out observed values which were observed for taxa with genetic perturbations """

    #TODO: Need to figure out what the options are for these

    def __init__(self):
        super(WildtypeFilter, self).__init__(('metadata', 'genetics', 'variation', ), [''])


class TemperatureRangeFilter(RangeFilter):
    """ Filters out observed values with temperatures that fall outside a specified range. """

    def __init__(self, min=float('nan'), max=float('nan')):
        """
        Args:
            min (:obj:`float`, optional): minimum value
            max (:obj:`float`, optional): maximum value
        """
        super(TemperatureRangeFilter, self).__init__(('metadata', 'environment', 'temperature', ), min=min, max=max)


class TemperatureNormalFilter(NormalFilter):
    """ Prioritizes observed values with temperatures that are close to `mean`. """

    def __init__(self, mean, std):
        """
        Args:
            mean (:obj:`float`): The mean of the distribution. This indicates the value at which the score will be 1.
            std (:obj:`float`): The standard deviation of the distribution. This determines how quickly the score falls to zero away from the mean.
        """
        super(TemperatureNormalFilter, self).__init__(('metadata', 'environment', 'temperature', ), mean, std)


class PhRangeFilter(RangeFilter):
    """ Filters out observed values with pHs that fall outside a specified range. """

    def __init__(self, min=float('nan'), max=float('nan')):
        """
        Args:
            min (:obj:`float`, optional): minimum value
            max (:obj:`float`, optional): maximum value
        """
        super(PhRangeFilter, self).__init__(('metadata', 'environment', 'ph', ), min=min, max=max)


class PhNormalFilter(NormalFilter):
    """ Prioritizes observed values with pHs that are close to `mean`. """

    def __init__(self, mean, std):
        """
        Args:
            mean (:obj:`float`): The mean of the distribution. This indicates the value at which the score will be 1.
            std (:obj:`float`): The standard deviation of the distribution. This determines how quickly the score
                falls to zero away from the mean.
        """
        super(PhNormalFilter, self).__init__(('metadata', 'environment', 'ph', ), mean, std)


class ConsensusGenerator(object):
    """ """

    def run(self, observed_results, method, weighted=True):
        """

        Args:
            observed_results (:obj:`list` of :obj:`data_model.ObservedValue`): list of
                observed values
            method (:obj:`str`): `mean`, `median`, or `mode`; desired average
                statistic
            weighted (:obj:`bool`, optional): if :obj:`True`, calculate the weighted
                average value

        Returns:
            :obj:`list` of :obj:`data_model.Consensus`: list of consensus values of
                the observed properties

        Raises:
            :obj:`ValueError`: if :obj:`method` is not one of `mean`, `median`,
                or `mode`
        """
        consensuses = []
        for observable, evidence in self.group_observed_results_by_properties(observed_results):
            norm_values, norm_errors, weights, units = self.normalize_observed_results(evidence)
            value, error, method = self.calc_average(norm_values, weights=weights if weighted else None, method=method)
            consensus.append(data_model.Consensus(
                observable=observable,
                value=value,
                error=error,
                units=units,
                evidence=evidence,
                method=method,
                user=getpass.getuser(),
                date=datetime.utcnow(),
            ))

        return consensuses

    def group_observed_results_by_properties(self, observed_results):
        """ Group observed values by their observed properties

        Args:
            observed_results (:obj:`list` of :obj:`data_model.ObservedValue`): list of
                observed values

        Returns:
            :obj:`list` of :obj:`tuple` of (:obj:`str`, :obj:`list` of :obj:`data_model.ObservedValue`):
                list of observed values, grouped by the observed property
        """
        grouped_obs = {}
        for obs, score in zip(filter_result.observed_results, filter_result.scores):
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

    def normalize_observed_results(self, observed_results):
        """ Normalize one or more observed values to SI units

        Args:
            observed_results (:obj:`list` of :obj:`data_model.ObservedValue`): list of
                observed values

        Returns:
            :obj:`tuple`:
                * :obj:`list` of :obj:`float`: normalized observed values
                * :obj:`list` of :obj:`float`: normalized errors of the observed values
                * :obj:`list` of :obj:`float`: weights of the observed values
                * :obj:`str`: units of the normalized observed values
        """
        norm_values = []
        norm_errors = []
        weights = []
        units = ''
        for ov in observed_results:
            value = ob.value * unit_registry(ob.units)
            error = ob.error * unit_registry(ob.units)
            values.append(value.to_base_units().magnitude)
            errors.append(error.to_base_units().magnitude)

            weights.append(ov.relevance)

    def calc_average(self, values, weights=None, method='mean'):
        """ Calculate the weighted or unweighted average of one of more values

        Args:
            values (:obj:`list` of :obj:`float`): list of normalized values
            weights (:obj:`list` of :obj:`float`, optional): weights of :obj:`values`
            method (:obj:`str`, optional): `mean`, `median`, or `mode`; the desired average of
                :obj:`values`

        Returns:
            :obj:`tuple` of :obj:`float`, :obj:`float`, :obj:`data_model.ConsensusMethod`:
                tuple of the average value, its uncertainty, and the method used to calculate
                the average value
        """
        # convert to numpy arrays
        values = 1. * numpy.array(values)
        if weights is not None:
            weights = 1. * numpy.array(weights)

        # ignore nan values
        tfs = numpy.logical_not(numpy.isnan(values))
        if not numpy.any(tfs):
            return (float('nan'), float('nan'), None)
        values = numpy.extract(tfs, values)
        if weights is not None:
            weights = numpy.extract(tfs, weights)

        # ignore weights if any are nan
        if weights is not None:
            if not all(w is not None and not numpy.isnan(w) for w in weights):
                weights = None

        # calculate average statistic
        if method == 'mean':
            if weights is None:
                value = numpy.average(values)
            else:
                value = numpy.average(values, weights=weights)

            if weights is None:
                method = data_model.ConsensusMethod.mean
            else:
                method = data_model.ConsensusMethod.weighted_mean

        elif method == 'median':
            if weights is None:
                value = numpy.median(values)
            else:
                value = wc_utils.util.stats.weighted_median(values, weights)

            if weights is None:
                method = data_model.ConsensusMethod.median
            else:
                method = data_model.ConsensusMethod.weighted_median

        elif method == 'mode':
            if weights is None:
                value = scipy.stats.mode(values).mode[0]
            else:
                value = wc_utils.util.stats.weighted_mode(values, weights)

            if weights is None:
                method = data_model.ConsensusMethod.mode
            else:
                method = data_model.ConsensusMethod.weighted_mode

        else:
            raise ValueError('Unsupported consensus method `{}`'.format(method))

        # calculate error
        # todo: is this the most informative error statistic?
        if values.size == 1:
            error = float('nan')
        else:
            if weights is None:
                error = numpy.std(values)
            else:
                mean = numpy.average(values, weights=weights)
                error = numpy.sqrt(numpy.sum(weights * numpy.power(values - mean, 2) / numpy.sum(weights)))

        # return calculated statistic and error
        return (float(value), float(error), method)
