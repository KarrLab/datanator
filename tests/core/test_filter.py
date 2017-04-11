""" Test of observation filtering and ordering

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-10
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import filter
from kinetic_datanator.core import observation
from numpy.testing import assert_almost_equal
from scipy.stats import norm
import unittest


class TestFilter(unittest.TestCase):

    @unittest.skip('write me')
    def test_TaxonomicDistanceFilter(self):
        pass  # todo

    def test_OptionsFilter(self):
        f = filter.OptionsFilter('strain.perturbations', [''])
        o = observation.Observation(strain=observation.Strain(perturbations=''))
        self.assertEqual(f.score(o), 1)
        o = observation.Observation(strain=observation.Strain(perturbations='wildtype'))
        self.assertEqual(f.score(o), -1)
        o = observation.Observation(strain=observation.Strain(perturbations='Delta gene-01'))
        self.assertEqual(f.score(o), -1)

        f = filter.OptionsFilter('strain.perturbations', ['wildtype'])
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

        f = filter.RangeFilter('environment.temperature', min_val=15., max_val=30.)
        self.assertEqual(f.score(obs(float('nan'))), -1)
        self.assertEqual(f.score(obs(10)), -1)
        self.assertEqual(f.score(obs(15.0)), 1)
        self.assertEqual(f.score(obs(15.01)), 1)
        self.assertEqual(f.score(obs(29.99)), 1)
        self.assertEqual(f.score(obs(30.0)), 1)
        self.assertEqual(f.score(obs(31)), -1)

        f = filter.RangeFilter('environment.temperature', min_val=float('nan'), max_val=30.)
        self.assertEqual(f.score(obs(float('nan'))), -1)
        self.assertEqual(f.score(obs(10)), 1)
        self.assertEqual(f.score(obs(15.0)), 1)
        self.assertEqual(f.score(obs(15.01)), 1)
        self.assertEqual(f.score(obs(29.99)), 1)
        self.assertEqual(f.score(obs(30.0)), 1)
        self.assertEqual(f.score(obs(31)), -1)

        f = filter.RangeFilter('environment.temperature', min_val=15., max_val=float('nan'))
        self.assertEqual(f.score(obs(float('nan'))), -1)
        self.assertEqual(f.score(obs(10)), -1)
        self.assertEqual(f.score(obs(15.0)), 1)
        self.assertEqual(f.score(obs(15.01)), 1)
        self.assertEqual(f.score(obs(29.99)), 1)
        self.assertEqual(f.score(obs(30.0)), 1)
        self.assertEqual(f.score(obs(31)), 1)

        f = filter.RangeFilter('environment.temperature')
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

        f = filter.TemperatureRangeFilter(min_val=15., max_val=30.)
        self.assertEqual(f.score(obs(10)), -1)
        self.assertEqual(f.score(obs(15.0)), 1)
        self.assertEqual(f.score(obs(15.01)), 1)
        self.assertEqual(f.score(obs(29.99)), 1)
        self.assertEqual(f.score(obs(30.0)), 1)
        self.assertEqual(f.score(obs(31)), -1)

    def test_PhRangeFilter(self):
        def obs(ph):
            return observation.Observation(environment=observation.Environment(ph=ph))

        f = filter.PhRangeFilter(min_val=5., max_val=9.)
        self.assertEqual(f.score(obs(3)), -1)
        self.assertEqual(f.score(obs(5.0)), 1)
        self.assertEqual(f.score(obs(5.01)), 1)
        self.assertEqual(f.score(obs(8.99)), 1)
        self.assertEqual(f.score(obs(9.0)), 1)
        self.assertEqual(f.score(obs(10)), -1)

    def test_NormalFilter(self):
        def obs(temperature):
            return observation.Observation(environment=observation.Environment(temperature=temperature))

        f = filter.NormalFilter('environment.temperature', mean=37, std=1)
        self.assertEqual(f.score(obs(37)), 1)
        assert_almost_equal(f.score(obs(36)), 2 * norm.cdf(-1), decimal=5)
        assert_almost_equal(f.score(obs(38)), 2 * norm.cdf(-1), decimal=5)
        assert_almost_equal(f.score(obs(35)), 2 * norm.cdf(-2), decimal=5)
        assert_almost_equal(f.score(obs(39)), 2 * norm.cdf(-2), decimal=5)

    def test_TemperatureNormalFilter(self):
        def obs(temperature):
            return observation.Observation(environment=observation.Environment(temperature=temperature))

        f = filter.TemperatureNormalFilter(mean=37, std=1)
        self.assertEqual(f.score(obs(37)), 1)
        assert_almost_equal(f.score(obs(36)), 2 * norm.cdf(-1), decimal=5)
        assert_almost_equal(f.score(obs(38)), 2 * norm.cdf(-1), decimal=5)
        assert_almost_equal(f.score(obs(35)), 2 * norm.cdf(-2), decimal=5)
        assert_almost_equal(f.score(obs(39)), 2 * norm.cdf(-2), decimal=5)

    def test_PhNormalFilter(self):
        def obs(ph):
            return observation.Observation(environment=observation.Environment(ph=ph))

        f = filter.PhNormalFilter(mean=7, std=1)
        self.assertEqual(f.score(obs(7)), 1)
        assert_almost_equal(f.score(obs(6)), 2 * norm.cdf(-1), decimal=5)
        assert_almost_equal(f.score(obs(8)), 2 * norm.cdf(-1), decimal=5)
        assert_almost_equal(f.score(obs(5)), 2 * norm.cdf(-2), decimal=5)
        assert_almost_equal(f.score(obs(9)), 2 * norm.cdf(-2), decimal=5)

    @unittest.skip('write me')
    def test_FilterResult(self):
        pass  # todo

    @unittest.skip('write me')
    def test_FilterRunner(self):
        pass  # todo
