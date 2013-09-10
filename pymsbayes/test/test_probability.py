#! /usr/bin/env python

import unittest
import os
import math

from pymsbayes.utils.probability import *
from pymsbayes.utils.stats import SampleSummarizer
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class ContinuousUniformDistributionTestCase(PyMsBayesTestCase):

    def test_value_error(self):
        self.assertRaises(ValueError, ContinuousUniformDistribution, 2, 1)

    def test_standard(self):
        d = ContinuousUniformDistribution(0, 1)
        self.assertTrue(isinstance(d, ContinuousUniformDistribution))
        self.assertEqual(d.name, 'uniform')
        self.assertEqual(d.minimum, 0)
        self.assertEqual(d.maximum, 1)
        self.assertEqual(d.mean, 0.5)
        self.assertEqual(d.mean, d.median)
        self.assertEqual(d.variance, float(1)/12) 
        self.assertEqual(d.probability_density(0), 1)
        self.assertEqual(d.probability_density(1), 1)
        self.assertEqual(d.probability_density(0.5), 1)
        self.assertEqual(d.probability_density(1.01), 0)
        self.assertEqual(d.probability_density(-0.01), 0)

    def test_draw(self):
        d = ContinuousUniformDistribution(0, 1)
        s = SampleSummarizer()
        for i in d.draw_iter(100000):
            s.add_sample(i)
        self.assertAlmostEqual(d.mean, s.mean, 2)
        self.assertAlmostEqual(d.variance, s.variance, 2)
        self.assertTrue(s.minimum >= d.minimum)
        self.assertTrue(s.maximum <= s.maximum)
        self.assertAlmostEqual(d.minimum, s.minimum, 4)
        self.assertAlmostEqual(d.maximum, s.maximum, 4)

class DiscreteUniformDistributionTestCase(PyMsBayesTestCase):

    def test_value_error_init(self):
        self.assertRaises(ValueError, DiscreteUniformDistribution, 2, 1)

    def test_value_error_init_float(self):
        self.assertRaises(ValueError, DiscreteUniformDistribution, 1.0, 2.0)

    def test_standard(self):
        d = DiscreteUniformDistribution(1, 3)
        self.assertTrue(isinstance(d, DiscreteUniformDistribution))
        self.assertEqual(d.name, 'discrete-uniform')
        self.assertEqual(d.minimum, 1)
        self.assertEqual(d.maximum, 3)
        self.assertEqual(d.mean, 2)
        self.assertEqual(d.mean, d.median)
        self.assertEqual(d.variance, float(8)/12) 
        self.assertEqual(d.probability_mass(0), 0)
        self.assertEqual(d.probability_mass(1), float(1)/3)
        self.assertEqual(d.probability_mass(2), float(1)/3)
        self.assertEqual(d.probability_mass(3), float(1)/3)
        self.assertEqual(d.probability_mass(4), 0)

    def test_float_mass_error(self):
        d = DiscreteUniformDistribution(1, 3)
        self.assertRaises(ValueError, d.probability_mass, 1.0)

    def test_draw(self):
        d = DiscreteUniformDistribution(1, 3)
        s = SampleSummarizer()
        for i in d.draw_iter(100000):
            s.add_sample(i)
        self.assertAlmostEqual(d.mean, s.mean, 2)
        self.assertAlmostEqual(d.variance, s.variance, 2)
        self.assertTrue(s.minimum >= d.minimum)
        self.assertTrue(s.maximum <= s.maximum)
        self.assertEqual(d.minimum, s.minimum)
        self.assertEqual(d.maximum, s.maximum)

class GammaDistributionTestCase(PyMsBayesTestCase):

    def test_value_error(self):
        self.assertRaises(ValueError, GammaDistribution, 1, 0)
        self.assertRaises(ValueError, GammaDistribution, 0, 1)

    def test_standard(self):
        d = GammaDistribution(2, 5)
        self.assertTrue(isinstance(d, GammaDistribution))
        self.assertEqual(d.name, 'gamma')
        self.assertEqual(d.minimum, 0.0)
        self.assertEqual(d.maximum, float('inf'))
        self.assertEqual(d.mean, 10.0)
        self.assertEqual(d.variance, 50.0)
        self.assertRaises(NotImplementedError, d._get_median)

    def test_draw(self):
        d = GammaDistribution(2, 0.01)
        s = SampleSummarizer()
        for i in d.draw_iter(100000):
            s.add_sample(i)
        self.assertAlmostEqual(d.mean, s.mean, 3)
        self.assertAlmostEqual(d.variance, s.variance, 5)
        self.assertTrue(s.minimum >= d.minimum)
        self.assertTrue(s.maximum <= s.maximum)
        self.assertAlmostEqual(d.minimum, s.minimum, 3)


class BetaDistributionTestCase(PyMsBayesTestCase):

    def test_value_error(self):
        self.assertRaises(ValueError, BetaDistribution, 1, 0)
        self.assertRaises(ValueError, BetaDistribution, 0, 1)

    def test_standard(self):
        d = BetaDistribution(1, 1)
        self.assertTrue(isinstance(d, BetaDistribution))
        self.assertEqual(d.name, 'beta')
        self.assertEqual(d.minimum, 0)
        self.assertEqual(d.maximum, 1)
        self.assertEqual(d.mean, float(1)/2)
        self.assertEqual(d.variance, float(1)/12) 
        self.assertRaises(NotImplementedError, d._get_median)

    def test_draw(self):
        d = BetaDistribution(1, 1)
        s = SampleSummarizer()
        for i in d.draw_iter(100000):
            s.add_sample(i)
        self.assertAlmostEqual(d.mean, s.mean, 2)
        self.assertAlmostEqual(d.variance, s.variance, 2)
        self.assertTrue(s.minimum >= d.minimum)
        self.assertTrue(s.maximum <= s.maximum)
        self.assertAlmostEqual(d.minimum, s.minimum, 4)
        self.assertAlmostEqual(d.maximum, s.maximum, 4)

class MultinomialTestCase(unittest.TestCase):
    def test_init(self):
        p = [0.25, 0.25, 0.25, 0.25]
        m = Multinomial(p)
        self.assertEqual(m.parameters, p)

    def test_parameters_error(self):
        p = [0.25, 0.25, 0.25, 0.24]
        self.assertRaises(ValueError, Multinomial, p)
        m = Multinomial([0.5, 0.5])
        self.assertRaises(ValueError, m._set_parameters, [0.33, 0.33, 0.32])

    def test_count_error(self):
        p = [0.25, 0.25, 0.25, 0.25]
        m = Multinomial(p)
        c = [2, 1, 1]
        self.assertRaises(ValueError, m._check_counts, c)
        self.assertRaises(ValueError, m.coefficient, c)
        self.assertRaises(ValueError, m.log_coefficient, c)
        self.assertRaises(ValueError, m.pmf, c)
        self.assertRaises(ValueError, m.log_pmf, c)

    def test_coeff(self):
        p = [float(1) / 3] * 3
        m = Multinomial(p)
        self.assertAlmostEqual(m.coefficient([3,0,0]), 1.0)
        self.assertAlmostEqual(m.coefficient([2,1,0]), 3.0)
        self.assertAlmostEqual(m.coefficient([1,0,2]), 3.0)
        self.assertAlmostEqual(m.coefficient([1,1,1]), 6.0)
        p = [0.25, 0.25, 0.25, 0.25]
        m = Multinomial(p)
        self.assertAlmostEqual(m.coefficient([1,1,1,1]), 24.0)

    def test_pmf(self):
        p = float(1) / 3
        params = [p] * 3
        m = Multinomial(params)
        self.assertAlmostEqual(m.pmf([3,0,0]), p**3)

        params = [0.5, 0.2, 0.3]
        m = Multinomial(params)
        self.assertAlmostEqual(m.pmf([3,0,0]), 0.5**3)
        self.assertAlmostEqual(m.pmf([0,3,0]), 0.2**3)
        self.assertAlmostEqual(m.pmf([0,0,3]), 0.3**3)
        self.assertAlmostEqual(m.pmf([2,2,2]),
                90 * ((0.5**2) * (0.2**2) * (0.3 ** 2)))

        params = [0.5, 0.2, 0.1, 0.2]
        m = Multinomial(params)
        self.assertAlmostEqual(m.pmf([2,0,1,3]),
                60 *((0.5**2) * (0.2**0) * (0.1**1) * (0.2**3)))

class GetProbabilityFromBayesFactorTestCase(unittest.TestCase):
    def test_bf10_pr1_22(self):
        p = get_probability_from_bayes_factor(10, float(1)/22)
        self.assertAlmostEqual(p, float(10)/31)

    def test_bf20_pr1_22(self):
        p = get_probability_from_bayes_factor(20, float(1)/22)
        self.assertAlmostEqual(p, float(20)/41)

    def test_bf01_pr1_5(self):
        p = get_probability_from_bayes_factor(0.1, float(1)/5)
        self.assertAlmostEqual(p, 0.025/1.025)

if __name__ == '__main__':
    unittest.main()

