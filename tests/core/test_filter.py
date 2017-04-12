""" Test of observation filtering and ordering

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-10
:Copyright: 2017, Karr Lab
:License: MIT
"""

from copy import copy
from kinetic_datanator.core import filter
from kinetic_datanator.core import observation
from numpy import testing as npt
from scipy.stats import norm
import numpy as np
import unittest


class TestFilter(unittest.TestCase):

    @unittest.skip('write me')
    def test_TaxonomicDistanceFilter(self):
        pass  # todo

    def test_OptionsFilter(self):
        f = filter.OptionsFilter(('strain', 'perturbations', ), [''])
        o = observation.Observation(strain=observation.Strain(perturbations=''))
        self.assertEqual(f.score(o), 1)
        o = observation.Observation(strain=observation.Strain(perturbations='wildtype'))
        self.assertEqual(f.score(o), -1)
        o = observation.Observation(strain=observation.Strain(perturbations='Delta gene-01'))
        self.assertEqual(f.score(o), -1)

        f = filter.OptionsFilter(('strain', 'perturbations', ), ['wildtype'])
        o = observation.Observation(strain=observation.Strain(perturbations=''))
        self.assertEqual(f.score(o), -1)
        o = observation.Observation(strain=observation.Strain(perturbations='wildtype'))
        self.assertEqual(f.score(o), 1)
        o = observation.Observation(strain=observation.Strain(perturbations='Delta gene-01'))
        self.assertEqual(f.score(o), -1)

    def test_WildtypeFilter(self):
        f = filter.WildtypeFilter()

        o = observation.Observation(strain=observation.Strain(perturbations=''))
        self.assertEqual(f.score(o), 1)

        o = observation.Observation(strain=observation.Strain(perturbations='Delta gene-01'))
        self.assertEqual(f.score(o), -1)

    def test_RangeFilter(self):
        def obs(temperature):
            return observation.Observation(environment=observation.Environment(temperature=temperature))

        f = filter.RangeFilter(('environment', 'temperature', ), min=15., max=30.)
        self.assertEqual(f.score(obs(float('nan'))), -1)
        self.assertEqual(f.score(obs(10)), -1)
        self.assertEqual(f.score(obs(15.0)), 1)
        self.assertEqual(f.score(obs(15.01)), 1)
        self.assertEqual(f.score(obs(29.99)), 1)
        self.assertEqual(f.score(obs(30.0)), 1)
        self.assertEqual(f.score(obs(31)), -1)

        f = filter.RangeFilter(('environment', 'temperature', ), min=float('nan'), max=30.)
        self.assertEqual(f.score(obs(float('nan'))), -1)
        self.assertEqual(f.score(obs(10)), 1)
        self.assertEqual(f.score(obs(15.0)), 1)
        self.assertEqual(f.score(obs(15.01)), 1)
        self.assertEqual(f.score(obs(29.99)), 1)
        self.assertEqual(f.score(obs(30.0)), 1)
        self.assertEqual(f.score(obs(31)), -1)

        f = filter.RangeFilter(('environment', 'temperature', ), min=15., max=float('nan'))
        self.assertEqual(f.score(obs(float('nan'))), -1)
        self.assertEqual(f.score(obs(10)), -1)
        self.assertEqual(f.score(obs(15.0)), 1)
        self.assertEqual(f.score(obs(15.01)), 1)
        self.assertEqual(f.score(obs(29.99)), 1)
        self.assertEqual(f.score(obs(30.0)), 1)
        self.assertEqual(f.score(obs(31)), 1)

        f = filter.RangeFilter(('environment', 'temperature', ))
        self.assertEqual(f.score(obs(float('nan'))), 1)
        self.assertEqual(f.score(obs(10)), 1)
        self.assertEqual(f.score(obs(15.0)), 1)
        self.assertEqual(f.score(obs(15.01)), 1)
        self.assertEqual(f.score(obs(29.99)), 1)
        self.assertEqual(f.score(obs(30.0)), 1)
        self.assertEqual(f.score(obs(31)), 1)

    def test_TemperatureRangeFilter(self):
        def obs(temperature):
            return observation.Observation(environment=observation.Environment(temperature=temperature))

        f = filter.TemperatureRangeFilter(min=15., max=30.)
        self.assertEqual(f.score(obs(10)), -1)
        self.assertEqual(f.score(obs(15.0)), 1)
        self.assertEqual(f.score(obs(15.01)), 1)
        self.assertEqual(f.score(obs(29.99)), 1)
        self.assertEqual(f.score(obs(30.0)), 1)
        self.assertEqual(f.score(obs(31)), -1)

    def test_PhRangeFilter(self):
        def obs(ph):
            return observation.Observation(environment=observation.Environment(ph=ph))

        f = filter.PhRangeFilter(min=5., max=9.)
        self.assertEqual(f.score(obs(3)), -1)
        self.assertEqual(f.score(obs(5.0)), 1)
        self.assertEqual(f.score(obs(5.01)), 1)
        self.assertEqual(f.score(obs(8.99)), 1)
        self.assertEqual(f.score(obs(9.0)), 1)
        self.assertEqual(f.score(obs(10)), -1)

    def test_NormalFilter(self):
        def obs(temperature):
            return observation.Observation(environment=observation.Environment(temperature=temperature))

        f = filter.NormalFilter(('environment', 'temperature', ), mean=37, std=1)
        self.assertEqual(f.score(obs(37)), 1)
        npt.assert_almost_equal(f.score(obs(36)), 2 * norm.cdf(-1), decimal=5)
        npt.assert_almost_equal(f.score(obs(38)), 2 * norm.cdf(-1), decimal=5)
        npt.assert_almost_equal(f.score(obs(35)), 2 * norm.cdf(-2), decimal=5)
        npt.assert_almost_equal(f.score(obs(39)), 2 * norm.cdf(-2), decimal=5)

    def test_TemperatureNormalFilter(self):
        def obs(temperature):
            return observation.Observation(environment=observation.Environment(temperature=temperature))

        f = filter.TemperatureNormalFilter(mean=37, std=1)
        self.assertEqual(f.score(obs(37)), 1)
        npt.assert_almost_equal(f.score(obs(36)), 2 * norm.cdf(-1), decimal=5)
        npt.assert_almost_equal(f.score(obs(38)), 2 * norm.cdf(-1), decimal=5)
        npt.assert_almost_equal(f.score(obs(35)), 2 * norm.cdf(-2), decimal=5)
        npt.assert_almost_equal(f.score(obs(39)), 2 * norm.cdf(-2), decimal=5)

    def test_PhNormalFilter(self):
        def obs(ph):
            return observation.Observation(environment=observation.Environment(ph=ph))

        f = filter.PhNormalFilter(mean=7, std=1)
        self.assertEqual(f.score(obs(7)), 1)
        npt.assert_almost_equal(f.score(obs(6)), 2 * norm.cdf(-1), decimal=5)
        npt.assert_almost_equal(f.score(obs(8)), 2 * norm.cdf(-1), decimal=5)
        npt.assert_almost_equal(f.score(obs(5)), 2 * norm.cdf(-2), decimal=5)
        npt.assert_almost_equal(f.score(obs(9)), 2 * norm.cdf(-2), decimal=5)

    def test_FilterResult(self):
        o = [
            observation.Observation(value=1),
            observation.Observation(value=2),
        ]
        f = [
            filter.RangeFilter(('environment', 'temperature', )),
            filter.RangeFilter(('environment', 'ph', )),
        ]
        s = [0.3, 0.5]

        filter.FilterResult(o, f, s, copy(o), copy(s), [1, 2])

    def test_FilterRunner_score(self):
        o = [
            observation.Observation(value=1, environment=observation.Environment(temperature=37, ph=7)),
            observation.Observation(value=1, environment=observation.Environment(temperature=37, ph=6)),
            observation.Observation(value=1, environment=observation.Environment(temperature=36, ph=7)),
            observation.Observation(value=1, environment=observation.Environment(temperature=36, ph=6)),
            observation.Observation(value=2, environment=observation.Environment(temperature=37, ph=7)),
            observation.Observation(value=2, environment=observation.Environment(temperature=37, ph=6)),
            observation.Observation(value=2, environment=observation.Environment(temperature=36, ph=7)),
            observation.Observation(value=2, environment=observation.Environment(temperature=36, ph=6)),
        ]
        f = [
            filter.TemperatureRangeFilter(min=36.5, max=37.5),
            filter.PhRangeFilter(min=6.5, max=7.5),
            filter.TemperatureNormalFilter(mean=37., std=1),
        ]

        runner = filter.FilterRunner(f[0])
        self.assertEqual(list(runner.score(o).ravel()), [1., 1., -1., -1., 1., 1., -1., -1.])

        runner = filter.FilterRunner(f[1])
        self.assertEqual(list(runner.score(o).ravel()), [1., -1., 1., -1., 1., -1., 1., -1.])

        runner = filter.FilterRunner(f[2])
        s = 2 * norm.cdf(-1)
        npt.assert_almost_equal(list(runner.score(o).ravel()), [1., 1., s, s, 1., 1., s, s], decimal=5)

    def test_FilterRunner_filter(self):
        o = [
            observation.Observation(value=1, environment=observation.Environment(temperature=37, ph=7)),
            observation.Observation(value=1, environment=observation.Environment(temperature=37, ph=6)),
            observation.Observation(value=1, environment=observation.Environment(temperature=36, ph=7)),
            observation.Observation(value=1, environment=observation.Environment(temperature=36, ph=6)),
            observation.Observation(value=2, environment=observation.Environment(temperature=37, ph=7)),
            observation.Observation(value=2, environment=observation.Environment(temperature=37, ph=6)),
            observation.Observation(value=2, environment=observation.Environment(temperature=36, ph=7)),
            observation.Observation(value=2, environment=observation.Environment(temperature=36, ph=6)),
        ]

        # one filter
        s = np.array([1, -1, 0.5, 1, -1, 0, -1, 0.3], ndmin=2).transpose()
        runner = filter.FilterRunner([])
        filt_o, filt_s, i_filt_o = runner.filter(o, s)

        self.assertEqual(filt_o, [o[0], o[2], o[3], o[5], o[7]])
        npt.assert_equal(filt_s, np.array([1, 0.5, 1, 0, 0.3], ndmin=2).transpose())
        self.assertEqual(i_filt_o, [0, 2, 3, 5, 7])

        # multiple filters
        s = np.array([[1, -1, 0.5, 1, -1, 0, -1, 0.3], [0, 0, 0, 0, 0, 0, 0, -1]], ndmin=2).transpose()
        runner = filter.FilterRunner([])
        filt_o, filt_s, i_filt_o = runner.filter(o, s)

        self.assertEqual(filt_o, [o[0], o[2], o[3], o[5]])
        npt.assert_equal(filt_s, np.array([[1, 0.5, 1, 0], [0, 0, 0, 0]], ndmin=2).transpose())
        self.assertEqual(i_filt_o, [0, 2, 3, 5])

    def test_FilterRunner_order(self):
        o = [
            observation.Observation(value=1, environment=observation.Environment(temperature=37, ph=7)),
            observation.Observation(value=1, environment=observation.Environment(temperature=37, ph=6)),
            observation.Observation(value=1, environment=observation.Environment(temperature=36, ph=7)),
            observation.Observation(value=1, environment=observation.Environment(temperature=36, ph=6)),
            observation.Observation(value=2, environment=observation.Environment(temperature=37, ph=7)),
            observation.Observation(value=2, environment=observation.Environment(temperature=37, ph=6)),
            observation.Observation(value=2, environment=observation.Environment(temperature=36, ph=7)),
            observation.Observation(value=2, environment=observation.Environment(temperature=36, ph=6)),
        ]

        # one filter
        s = np.array([0.9, -1, 0.5, 1, -1.1, 0, -1.3, 0.3], ndmin=2).transpose()
        runner = filter.FilterRunner([])
        order_o, order_s, i_order_o = runner.order(o, s)

        self.assertEqual(order_o, [o[3], o[0], o[2], o[7], o[5], o[1], o[4], o[6]])
        npt.assert_equal(order_s, np.array([1, 0.9, 0.5, 0.3, 0, -1, -1.1, -1.3], ndmin=2).transpose())
        self.assertEqual(i_order_o, [3, 0, 2, 7, 5, 1, 4, 6])

        # multiple filters
        s = np.array([[0.9, -1, 0.5, 1, -1.1, 0, -1.3, 0.3], [0, 0, 0, 0, 0, 0, 0, 0]], ndmin=2).transpose()
        runner = filter.FilterRunner([])
        order_o, order_s, i_order_o = runner.order(o, s)

        self.assertEqual(order_o, [o[3], o[0], o[2], o[7], o[5], o[1], o[4], o[6]])
        npt.assert_equal(order_s, np.array([[1, 0.9, 0.5, 0.3, 0, -1, -1.1, -1.3], [0, 0, 0, 0, 0, 0, 0, 0]], ndmin=2).transpose())
        self.assertEqual(i_order_o, [3, 0, 2, 7, 5, 1, 4, 6])

    def test_FilterRunner_run(self):
        o = [
            observation.Observation(value=1, environment=observation.Environment(temperature=37, ph=7)),
            observation.Observation(value=1, environment=observation.Environment(temperature=37, ph=6)),
            observation.Observation(value=1, environment=observation.Environment(temperature=36, ph=7)),
            observation.Observation(value=1, environment=observation.Environment(temperature=36, ph=6)),
            observation.Observation(value=2, environment=observation.Environment(temperature=37.1, ph=7)),
            observation.Observation(value=2, environment=observation.Environment(temperature=37.1, ph=6)),
            observation.Observation(value=2, environment=observation.Environment(temperature=37.01, ph=7)),
            observation.Observation(value=2, environment=observation.Environment(temperature=37.01, ph=6)),
        ]
        f = [
            filter.TemperatureNormalFilter(mean=37, std=1),
            filter.PhRangeFilter(min=6.5, max=7.5),
        ]
        runner = filter.FilterRunner(f)

        # return_info=True
        result = runner.run(o, return_info=True)
        self.assertEqual(result.observations, o)
        self.assertEqual(result.filters, f)
        s1 = 2 * norm.cdf(-1)
        s01 = 2 * norm.cdf(-0.1)
        s001 = 2 * norm.cdf(-0.01)
        npt.assert_almost_equal(result.scores, np.array([[1., 1., s1, s1, s01, s01, s001, s001],
                                                         [1., -1., 1., -1., 1., -1., 1., -1.]], ndmin=2).transpose())

        self.assertEqual(result.ordered_observations, [o[0], o[6], o[4], o[2]])
        npt.assert_almost_equal(result.ordered_scores, np.array([[1., s001, s01, s1], [1., 1., 1., 1.]], ndmin=2).transpose())
        self.assertEqual(result.ordered_observation_indices, [0, 6, 4, 2])

        # return_info=False
        ordered_observations = runner.run(o, return_info=False)
        self.assertEqual(ordered_observations, [o[0], o[6], o[4], o[2]])
