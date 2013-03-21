#! /usr/bin/env python

import unittest
import os
import math

from pymsbayes.utils.stats import SampleSummarizer
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class SampleSummarizerTestCase(PyMsBayesTestCase):

    def test_init(self):
        ss = SampleSummarizer('test')
        self.assertEqual(ss.name, 'test')
        self.assertEqual(ss.minimum, None)
        self.assertEqual(ss.maximum, None)
        self.assertEqual(ss.mean, None)
        self.assertEqual(ss.variance, None)
        self.assertEqual(ss.std_deviation, None)
        self.assertEqual(ss.pop_variance, None)

    def test_add_one_sample(self):
        ss = SampleSummarizer('test')
        ss.add_sample(1)
        self.assertEqual(ss.name, 'test')
        self.assertEqual(ss.minimum, 1)
        self.assertEqual(ss.maximum, 1)
        self.assertApproxEqual(ss.mean, 1.0, 1e-9)
        self.assertEqual(ss.variance, float('inf'))
        self.assertEqual(ss.std_deviation, float('inf'))
        self.assertEqual(ss.pop_variance, 0)

        ss = SampleSummarizer('test')
        ss.add_sample(3.45)
        self.assertEqual(ss.name, 'test')
        self.assertEqual(ss.minimum, 3.45)
        self.assertEqual(ss.maximum, 3.45)
        self.assertApproxEqual(ss.mean, 3.45, 1e-9)
        self.assertEqual(ss.variance, float('inf'))
        self.assertEqual(ss.std_deviation, float('inf'))
        self.assertEqual(ss.pop_variance, 0)

    def test_update_samples(self):
        ss = SampleSummarizer('test')
        ss.update_samples([1.0, 2.0, 3.0])
        self.assertEqual(ss.name, 'test')
        self.assertEqual(ss.minimum, 1.0)
        self.assertEqual(ss.maximum, 3.0)
        self.assertApproxEqual(ss.mean, 2.0, 1e-9)
        self.assertApproxEqual(ss.variance, 1.0, 1e-9)
        self.assertEqual(ss.std_deviation, math.sqrt(1.0), 1e-9)
        self.assertApproxEqual(ss.pop_variance, 2/float(3), 1e-9)

if __name__ == '__main__':
    unittest.main()

