""" Test of observation filtering and ordering

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-10
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model
from kinetic_datanator.core import filter
from kinetic_datanator.util import warning_util
import copy
import math
import numpy
import scipy.stats
import unittest


warning_util.disable_warnings()


class TestFilters(unittest.TestCase):

    def test_TaxonomicDistanceFilter(self):
        # example 1
        f = filter.TaxonomicDistanceFilter('Mycoplasma pneumoniae M129')
        self.assertEqual(f.max, 8.)

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(taxon='Mycoplasma pneumoniae M129')))
        self.assertEqual(f.score(ov), 1.)

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(taxon='Mycoplasma pneumoniae')))
        self.assertEqual(f.scale, 8./5)
        self.assertEqual(f.score(ov), math.exp(-1/f.scale))

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(taxon='Escherichia coli')))
        self.assertEqual(f.score(ov), math.exp(-5))

        # example 2
        f = filter.TaxonomicDistanceFilter('Mycoplasma pneumoniae M129', max=5)
        self.assertEqual(f.max, 5.)

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(taxon='Mycoplasma pneumoniae M129')))
        self.assertEqual(f.score(ov), 1.)

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(taxon='Mycoplasma pneumoniae')))
        self.assertEqual(f.score(ov), math.exp(-1/f.scale))

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(taxon='Escherichia coli')))
        self.assertEqual(f.score(ov), -1)

        # example 3
        f = filter.TaxonomicDistanceFilter('Mycoplasma genitalium')
        self.assertEqual(f.max, 7.)

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(taxon='Mycoplasma genitalium G37')))
        self.assertEqual(f.score(ov), 1.)

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(taxon='Mycoplasma genitalium')))
        self.assertEqual(f.score(ov), 1)

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(taxon='Mycoplasma')))
        self.assertEqual(f.score(ov), math.exp(-1/f.scale))

    def test_OptionsFilter(self):
        f = filter.OptionsFilter(('observation', 'genetics', 'variation', ), [''])
        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(variation='')))
        self.assertEqual(f.score(ov), 1)
        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(variation='wildtype')))
        self.assertEqual(f.score(ov), -1)
        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(variation='Delta gene-01')))
        self.assertEqual(f.score(ov), -1)

        f = filter.OptionsFilter(('observation', 'genetics', 'variation', ), ['wildtype'])
        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(variation='')))
        self.assertEqual(f.score(ov), -1)
        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(variation='wildtype')))
        self.assertEqual(f.score(ov), 1)
        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(variation='Delta gene-01')))
        self.assertEqual(f.score(ov), -1)

    def test_WildtypeFilter(self):
        f = filter.WildtypeFilter()

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(variation='')))
        self.assertEqual(f.score(ov), 1)

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(variation='Delta gene-01')))
        self.assertEqual(f.score(ov), -1)

    def test_SpecieMolecularSimilarityFilter(self):
        adp = 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'
        atp = 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'
        h2o = 'O'

        f = filter.SpecieMolecularSimilarityFilter(atp)

        ov = data_model.ObservedValue(observable=data_model.Specie(structure=adp))
        numpy.testing.assert_almost_equal(f.score(ov), 0.955, decimal=3)

        ov = data_model.ObservedValue(observable=data_model.Specie(structure=h2o))
        numpy.testing.assert_almost_equal(f.score(ov), 0, decimal=3)

    def test_ReactionSimilarityFilter(self):
        def get_reaction(pi='InChI=1S/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)/p-2', ec='1.1.1.1'):
            atp = data_model.Specie(
                id='atp', structure='InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-4/t4-,6-,7-,10-/m1/s1')
            h2o = data_model.Specie(id='h2o', structure='InChI=1S/H2O/h1H2')
            adp = data_model.Specie(
                id='adp', structure='InChI=1S/C10H15N5O10P2/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(24-10)1-23-27(21,22)25-26(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H2,11,12,13)(H2,18,19,20)/p-3/t4-,6-,7-,10-/m1/s1')
            pi = data_model.Specie(id='pi', structure=pi)
            h = data_model.Specie(id='h', structure='InChI=1S/p+1/i/hH')

            return data_model.Reaction(
                participants=[
                    data_model.ReactionParticipant(coefficient=-1, specie=atp),
                    data_model.ReactionParticipant(coefficient=-1, specie=h2o),
                    data_model.ReactionParticipant(coefficient=1, specie=adp),
                    data_model.ReactionParticipant(coefficient=1, specie=pi),
                    data_model.ReactionParticipant(coefficient=1, specie=h),
                ],
                cross_references=[
                    data_model.Resource(namespace='ec-code', id=ec)
                ])

        f = filter.ReactionSimilarityFilter(get_reaction(), min_ec_level=3, scale=1)

        # same participants
        ov = data_model.ObservedValue(observable=get_reaction())
        numpy.testing.assert_almost_equal(f.score(ov), 1, decimal=3)

        # similiar participants
        ov = data_model.ObservedValue(observable=get_reaction(pi='InChI=1S/H3O4P/c1-5(2,3)4'))
        numpy.testing.assert_almost_equal(f.score(ov), 1, decimal=3)

        # different participants, same 4-digit EC
        ov = data_model.ObservedValue(observable=get_reaction(pi='InChI=1S/H4O4P/c1-5(2,3)4'))
        numpy.testing.assert_almost_equal(f.score(ov), math.exp(-1), decimal=3)

        # different participants, same 3-digit EC
        ov = data_model.ObservedValue(observable=get_reaction(pi='InChI=1S/H4O4P/c1-5(2,3)4', ec='1.1.1.2'))
        numpy.testing.assert_almost_equal(f.score(ov), math.exp(-2), decimal=3)

        ov = data_model.ObservedValue(observable=get_reaction(pi='InChI=1S/H4O4P/c1-5(2,3)4', ec='1.1.1'))
        numpy.testing.assert_almost_equal(f.score(ov), math.exp(-2), decimal=3)

        ov = data_model.ObservedValue(observable=get_reaction(pi='InChI=1S/H4O4P/c1-5(2,3)4', ec='1.1.1.'))
        numpy.testing.assert_almost_equal(f.score(ov), math.exp(-2), decimal=3)

        ov = data_model.ObservedValue(observable=get_reaction(pi='InChI=1S/H4O4P/c1-5(2,3)4', ec='1.1.1.-'))
        numpy.testing.assert_almost_equal(f.score(ov), math.exp(-2), decimal=3)

        # different participants, same 2-digit EC
        ov = data_model.ObservedValue(observable=get_reaction(pi='InChI=1S/H4O4P/c1-5(2,3)4', ec='1.1.2.1'))
        numpy.testing.assert_almost_equal(f.score(ov), -1, decimal=3)

        # target reaction only has 3 digits
        f1 = filter.ReactionSimilarityFilter(get_reaction(ec='1.1.1'), min_ec_level=3, scale=1)
        f2 = filter.ReactionSimilarityFilter(get_reaction(ec='1.1.1.'), min_ec_level=3, scale=1)

        ov = data_model.ObservedValue(observable=get_reaction(pi='InChI=1S/H4O4P/c1-5(2,3)4', ec='1.1.1'))
        numpy.testing.assert_almost_equal(f1.score(ov), math.exp(-2), decimal=3)
        numpy.testing.assert_almost_equal(f2.score(ov), math.exp(-2), decimal=3)

        ov = data_model.ObservedValue(observable=get_reaction(pi='InChI=1S/H4O4P/c1-5(2,3)4', ec='1.1.1.1'))
        numpy.testing.assert_almost_equal(f1.score(ov), math.exp(-2), decimal=3)
        numpy.testing.assert_almost_equal(f2.score(ov), math.exp(-2), decimal=3)

        ov = data_model.ObservedValue(observable=get_reaction(pi='InChI=1S/H4O4P/c1-5(2,3)4', ec='1.1.2.1'))
        numpy.testing.assert_almost_equal(f1.score(ov), -1, decimal=3)
        numpy.testing.assert_almost_equal(f2.score(ov), -1, decimal=3)

    def test_RangeFilter(self):
        def obs_val(temperature):
            return data_model.ObservedValue(
                observation=data_model.Observation(environment=data_model.Environment(temperature=temperature)))

        f = filter.RangeFilter(('observation', 'environment', 'temperature', ), min=15., max=30.)
        self.assertEqual(f.score(obs_val(float('nan'))), -1)
        self.assertEqual(f.score(obs_val(10)), -1)
        self.assertEqual(f.score(obs_val(15.0)), 1)
        self.assertEqual(f.score(obs_val(15.01)), 1)
        self.assertEqual(f.score(obs_val(29.99)), 1)
        self.assertEqual(f.score(obs_val(30.0)), 1)
        self.assertEqual(f.score(obs_val(31)), -1)

        f = filter.RangeFilter(('observation', 'environment', 'temperature', ), min=float('nan'), max=30.)
        self.assertEqual(f.score(obs_val(float('nan'))), -1)
        self.assertEqual(f.score(obs_val(10)), 1)
        self.assertEqual(f.score(obs_val(15.0)), 1)
        self.assertEqual(f.score(obs_val(15.01)), 1)
        self.assertEqual(f.score(obs_val(29.99)), 1)
        self.assertEqual(f.score(obs_val(30.0)), 1)
        self.assertEqual(f.score(obs_val(31)), -1)

        f = filter.RangeFilter(('observation', 'environment', 'temperature', ), min=15., max=float('nan'))
        self.assertEqual(f.score(obs_val(float('nan'))), -1)
        self.assertEqual(f.score(obs_val(10)), -1)
        self.assertEqual(f.score(obs_val(15.0)), 1)
        self.assertEqual(f.score(obs_val(15.01)), 1)
        self.assertEqual(f.score(obs_val(29.99)), 1)
        self.assertEqual(f.score(obs_val(30.0)), 1)
        self.assertEqual(f.score(obs_val(31)), 1)

        f = filter.RangeFilter(('observation', 'environment', 'temperature', ))
        self.assertEqual(f.score(obs_val(float('nan'))), 1)
        self.assertEqual(f.score(obs_val(10)), 1)
        self.assertEqual(f.score(obs_val(15.0)), 1)
        self.assertEqual(f.score(obs_val(15.01)), 1)
        self.assertEqual(f.score(obs_val(29.99)), 1)
        self.assertEqual(f.score(obs_val(30.0)), 1)
        self.assertEqual(f.score(obs_val(31)), 1)

    def test_TemperatureRangeFilter(self):
        def obs_val(temperature):
            return data_model.ObservedValue(
                observation=data_model.Observation(environment=data_model.Environment(temperature=temperature)))

        f = filter.TemperatureRangeFilter(min=15., max=30.)
        self.assertEqual(f.score(obs_val(10)), -1)
        self.assertEqual(f.score(obs_val(15.0)), 1)
        self.assertEqual(f.score(obs_val(15.01)), 1)
        self.assertEqual(f.score(obs_val(29.99)), 1)
        self.assertEqual(f.score(obs_val(30.0)), 1)
        self.assertEqual(f.score(obs_val(31)), -1)

    def test_PhRangeFilter(self):
        def obs_val(ph):
            return data_model.ObservedValue(
                observation=data_model.Observation(environment=data_model.Environment(ph=ph)))

        f = filter.PhRangeFilter(min=5., max=9.)
        self.assertEqual(f.score(obs_val(3)), -1)
        self.assertEqual(f.score(obs_val(5.0)), 1)
        self.assertEqual(f.score(obs_val(5.01)), 1)
        self.assertEqual(f.score(obs_val(8.99)), 1)
        self.assertEqual(f.score(obs_val(9.0)), 1)
        self.assertEqual(f.score(obs_val(10)), -1)

    def test_NormalFilter(self):
        def obs_val(temperature):
            return data_model.ObservedValue(
                observation=data_model.Observation(environment=data_model.Environment(temperature=temperature)))

        f = filter.NormalFilter(('observation', 'environment', 'temperature', ), mean=37, std=1)
        self.assertEqual(f.score(obs_val(37)), 1)
        numpy.testing.assert_almost_equal(f.score(obs_val(36)), 2 * scipy.stats.norm.cdf(-1), decimal=5)
        numpy.testing.assert_almost_equal(f.score(obs_val(38)), 2 * scipy.stats.norm.cdf(-1), decimal=5)
        numpy.testing.assert_almost_equal(f.score(obs_val(35)), 2 * scipy.stats.norm.cdf(-2), decimal=5)
        numpy.testing.assert_almost_equal(f.score(obs_val(39)), 2 * scipy.stats.norm.cdf(-2), decimal=5)

    def test_TemperatureNormalFilter(self):
        def obs_val(temperature):
            return data_model.ObservedValue(
                observation=data_model.Observation(environment=data_model.Environment(temperature=temperature)))

        f = filter.TemperatureNormalFilter(mean=37, std=1)
        self.assertEqual(f.score(obs_val(37)), 1)
        numpy.testing.assert_almost_equal(f.score(obs_val(36)), 2 * scipy.stats.norm.cdf(-1), decimal=5)
        numpy.testing.assert_almost_equal(f.score(obs_val(38)), 2 * scipy.stats.norm.cdf(-1), decimal=5)
        numpy.testing.assert_almost_equal(f.score(obs_val(35)), 2 * scipy.stats.norm.cdf(-2), decimal=5)
        numpy.testing.assert_almost_equal(f.score(obs_val(39)), 2 * scipy.stats.norm.cdf(-2), decimal=5)

    def test_PhNormalFilter(self):
        def obs_val(ph):
            return data_model.ObservedValue(
                observation=data_model.Observation(environment=data_model.Environment(ph=ph)))

        f = filter.PhNormalFilter(mean=7, std=1)
        self.assertEqual(f.score(obs_val(7)), 1)
        numpy.testing.assert_almost_equal(f.score(obs_val(6)), 2 * scipy.stats.norm.cdf(-1), decimal=5)
        numpy.testing.assert_almost_equal(f.score(obs_val(8)), 2 * scipy.stats.norm.cdf(-1), decimal=5)
        numpy.testing.assert_almost_equal(f.score(obs_val(5)), 2 * scipy.stats.norm.cdf(-2), decimal=5)
        numpy.testing.assert_almost_equal(f.score(obs_val(9)), 2 * scipy.stats.norm.cdf(-2), decimal=5)

    def test_ExponentialFilter(self):
        def obs_val(temperature):
            return data_model.ObservedValue(
                observation=data_model.Observation(environment=data_model.Environment(temperature=temperature)))

        f = filter.ExponentialFilter(('observation', 'environment', 'temperature', ), center=1., scale=1.)
        self.assertEqual(f.score(obs_val(1)), 1)
        numpy.testing.assert_almost_equal(f.score(obs_val(100)), 0, decimal=5)
        numpy.testing.assert_almost_equal(f.score(obs_val(2)), math.exp(-1), decimal=5)
        numpy.testing.assert_almost_equal(f.score(obs_val(3)), math.exp(-2), decimal=5)


class TestFilterResult(unittest.TestCase):

    def test(self):
        ov = [
            data_model.ObservedValue(value=1),
            data_model.ObservedValue(value=2),
        ]
        f = [
            filter.RangeFilter(('environment', 'temperature', )),
            filter.RangeFilter(('environment', 'ph', )),
        ]
        s = [0.3, 0.5]

        filter_result = filter.FilterResult(ov, f, s, copy.copy(ov), copy.copy(s), [1, 2])
        self.assertEqual(filter_result.observed_values, ov)
        self.assertEqual(filter_result.filters, f)
        self.assertEqual(filter_result.scores, s)
        self.assertEqual(filter_result.ordered_observed_values, ov)
        self.assertEqual(filter_result.ordered_scores, s)
        self.assertEqual(filter_result.ordered_observed_value_indices, [1, 2])


class TestFilterRunner(unittest.TestCase):

    def test_score(self):
        ov = [
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=7))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=6))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=7))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=6))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=7))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=6))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=7))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=6))),
        ]
        f = [
            filter.TemperatureRangeFilter(min=36.5, max=37.5),
            filter.PhRangeFilter(min=6.5, max=7.5),
            filter.TemperatureNormalFilter(mean=37., std=1),
        ]

        runner = filter.FilterRunner(f[0])
        self.assertEqual(list(runner.score(ov).ravel()), [1., 1., -1., -1., 1., 1., -1., -1.])

        runner = filter.FilterRunner(f[1])
        self.assertEqual(list(runner.score(ov).ravel()), [1., -1., 1., -1., 1., -1., 1., -1.])

        runner = filter.FilterRunner(f[2])
        s = 2 * scipy.stats.norm.cdf(-1)
        numpy.testing.assert_almost_equal(list(runner.score(ov).ravel()), [1., 1., s, s, 1., 1., s, s], decimal=5)

    def test_filter(self):
        ov = [
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=7))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=6))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=7))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=6))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=7))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=6))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=7))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=6))),
        ]

        # one filter
        s = numpy.array([1, -1, 0.5, 1, -1, 0, -1, 0.3], ndmin=2).transpose()
        runner = filter.FilterRunner([])
        filt_o, filt_s, i_filt_o = runner.filter(ov, s)

        self.assertEqual(filt_o, [ov[0], ov[2], ov[3], ov[5], ov[7]])
        numpy.testing.assert_equal(filt_s, numpy.array([1, 0.5, 1, 0, 0.3], ndmin=2).transpose())
        self.assertEqual(i_filt_o, [0, 2, 3, 5, 7])

        # multiple filters
        s = numpy.array([[1, -1, 0.5, 1, -1, 0, -1, 0.3], [0, 0, 0, 0, 0, 0, 0, -1]], ndmin=2).transpose()
        runner = filter.FilterRunner([])
        filt_o, filt_s, i_filt_o = runner.filter(ov, s)

        self.assertEqual(filt_o, [ov[0], ov[2], ov[3], ov[5]])
        numpy.testing.assert_equal(filt_s, numpy.array([[1, 0.5, 1, 0], [0, 0, 0, 0]], ndmin=2).transpose())
        self.assertEqual(i_filt_o, [0, 2, 3, 5])

    def test_order(self):
        ov = [
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=7))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=6))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=7))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=6))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=7))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=6))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=7))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=6))),
        ]

        # one filter
        s = numpy.array([0.9, -1, 0.5, 1, -1.1, 0, -1.3, 0.3], ndmin=2).transpose()
        runner = filter.FilterRunner([])
        order_o, order_s, i_order_o = runner.order(ov, s)

        self.assertEqual(order_o, [ov[3], ov[0], ov[2], ov[7], ov[5], ov[1], ov[4], ov[6]])
        numpy.testing.assert_equal(order_s, numpy.array([1, 0.9, 0.5, 0.3, 0, -1, -1.1, -1.3], ndmin=2).transpose())
        self.assertEqual(i_order_o, [3, 0, 2, 7, 5, 1, 4, 6])

        # multiple filters
        s = numpy.array([[0.9, -1, 0.5, 1, -1.1, 0, -1.3, 0.3], [0, 0, 0, 0, 0, 0, 0, 0]], ndmin=2).transpose()
        runner = filter.FilterRunner([])
        order_o, order_s, i_order_o = runner.order(ov, s)

        self.assertEqual(order_o, [ov[3], ov[0], ov[2], ov[7], ov[5], ov[1], ov[4], ov[6]])
        numpy.testing.assert_equal(order_s, numpy.array(
            [[1, 0.9, 0.5, 0.3, 0, -1, -1.1, -1.3], [0, 0, 0, 0, 0, 0, 0, 0]], ndmin=2).transpose())
        self.assertEqual(i_order_o, [3, 0, 2, 7, 5, 1, 4, 6])

    def test_run(self):
        ov = [
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=7))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=6))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=7))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=6))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37.1, ph=7))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37.1, ph=6))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37.01, ph=7))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37.01, ph=6))),
        ]
        f = [
            filter.TemperatureNormalFilter(mean=37, std=1),
            filter.PhRangeFilter(min=6.5, max=7.5),
        ]
        runner = filter.FilterRunner(f)

        # return_info=True
        result = runner.run(ov, return_info=True)
        self.assertEqual(result.observed_values, ov)
        self.assertEqual(result.filters, f)
        s1 = 2 * scipy.stats.norm.cdf(-1)
        s01 = 2 * scipy.stats.norm.cdf(-0.1)
        s001 = 2 * scipy.stats.norm.cdf(-0.01)
        numpy.testing.assert_almost_equal(result.scores, numpy.array([[1., 1., s1, s1, s01, s01, s001, s001],
                                                                      [1., -1., 1., -1., 1., -1., 1., -1.]], ndmin=2).transpose())

        self.assertEqual(result.ordered_observed_values, [ov[0], ov[6], ov[4], ov[2]])
        numpy.testing.assert_almost_equal(result.ordered_scores, numpy.array([[1., s001, s01, s1], [1., 1., 1., 1.]], ndmin=2).transpose())
        self.assertEqual(result.ordered_observed_value_indices, [0, 6, 4, 2])

        # return_info=False
        ordered_observed_values = runner.run(ov, return_info=False)
        self.assertEqual(ordered_observed_values, [ov[0], ov[6], ov[4], ov[2]])
