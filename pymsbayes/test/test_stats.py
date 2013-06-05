#! /usr/bin/env python

import unittest
import os
import math

from pymsbayes.utils.stats import SampleSummarizer, SampleSummary
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

class SampleSummaryTestCase(PyMsBayesTestCase):

    def test_update(self):
        x = [0.1, -3.2, 3.5, 11.4, -2.3, 3.3, -5.6, 7.8, 2.9, -9.3]
        mn = sum(x) / len(x)
        ss = 0.0
        for i in x:
            ss += ((i - mn)**2)
        v = ss / (len(x) - 1)
        
        summarizer1 = SampleSummarizer()
        summarizer2 = SampleSummarizer()
        summarizer1.update_samples(x[:3])
        summarizer2.update_samples(x[3:])

        summary1 = SampleSummary(
                sample_size = summarizer1.n,
                mean = summarizer1.mean,
                variance = summarizer1.variance)
        summary2 = SampleSummary(
                sample_size = summarizer2.n,
                mean = summarizer2.mean,
                variance = summarizer2.variance)
        summary1.update(summary2)
        self.assertEqual(summary1.n, len(x))
        self.assertAlmostEqual(summary1.mean, mn)
        self.assertAlmostEqual(summary1.variance, v)
        self.assertAlmostEqual(summary1.std_deviation, math.sqrt(v))

    def test_update_default_init(self):
        x = [0.1, -3.2, 3.5, 11.4, -2.3, 3.3, -5.6, 7.8, 2.9, -9.3]
        mn = sum(x) / len(x)
        ss = 0.0
        for i in x:
            ss += ((i - mn)**2)
        v = ss / (len(x) - 1)
        
        summarizer1 = SampleSummarizer()
        summarizer2 = SampleSummarizer()
        summarizer1.update_samples(x[:3])
        summarizer2.update_samples(x[3:])

        summary1 = SampleSummary(
                sample_size = summarizer1.n,
                mean = summarizer1.mean,
                variance = summarizer1.variance)
        summary2 = SampleSummary(
                sample_size = summarizer2.n,
                mean = summarizer2.mean,
                variance = summarizer2.variance)
        s = SampleSummary()
        s.update(summary1)
        s.update(summary2)
        self.assertEqual(s.n, len(x))
        self.assertAlmostEqual(s.mean, mn)
        self.assertAlmostEqual(s.variance, v)
        self.assertAlmostEqual(s.std_deviation, math.sqrt(v))

        s = SampleSummary(sample_size=0, mean=12.1, variance=33.5)
        s.update(summary1)
        s.update(summary2)
        self.assertEqual(s.n, len(x))
        self.assertAlmostEqual(s.mean, mn)
        self.assertAlmostEqual(s.variance, v)
        self.assertAlmostEqual(s.std_deviation, math.sqrt(v))

    def test_update_iter(self):
        x1 = [0.1, -3.2, 3.5, 11.4, -2.3, 3.3, -5.6, 7.8, 2.9, -9.3]
        x2 = [11.1, -13.2, 23.5]
        x3 = [-0.1, 13.2, 3.25, 13.4, 2.3, 3.13, -15.6, -17.8, 2.19, 19.3]
        summarizer = SampleSummarizer()
        summarizer1 = SampleSummarizer()
        summarizer2 = SampleSummarizer()
        summarizer3 = SampleSummarizer()
        summarizer.update_samples(x1 + x1 + x2 + x3)
        summarizer1.update_samples(x1)
        summarizer2.update_samples(x2)
        summarizer3.update_samples(x3)

        summary1 = SampleSummary(
                sample_size = summarizer1.n,
                mean = summarizer1.mean,
                variance = summarizer1.variance)
        summary2 = SampleSummary(
                sample_size = summarizer2.n,
                mean = summarizer2.mean,
                variance = summarizer2.variance)
        summary3 = SampleSummary(
                sample_size = summarizer3.n,
                mean = summarizer3.mean,
                variance = summarizer3.variance)

        summary1.update_iter([summary1, summary2, summary3])
        self.assertEqual(summary1.n, summarizer.n)
        self.assertAlmostEqual(summary1.mean, summarizer.mean)
        self.assertAlmostEqual(summary1.variance, summarizer.variance)
        self.assertAlmostEqual(summary1.std_deviation, summarizer.std_deviation)
        
    def test_update_iter_with_summarizers(self):
        x1 = [0.1, -3.2, 3.5, 11.4, -2.3, 3.3, -5.6, 7.8, 2.9, -9.3]
        x2 = [11.1, -13.2, 23.5]
        x3 = [-0.1, 13.2, 3.25, 13.4, 2.3, 3.13, -15.6, -17.8, 2.19, 19.3]
        summarizer = SampleSummarizer()
        summarizer1 = SampleSummarizer()
        summarizer2 = SampleSummarizer()
        summarizer3 = SampleSummarizer()
        summarizer.update_samples(x1 + x1 + x2 + x3)
        summarizer1.update_samples(x1)
        summarizer2.update_samples(x2)
        summarizer3.update_samples(x3)

        summary1 = SampleSummary(
                sample_size = summarizer1.n,
                mean = summarizer1.mean,
                variance = summarizer1.variance)

        summary1.update_iter([summarizer1, summarizer2, summarizer3])
        self.assertEqual(summary1.n, summarizer.n)
        self.assertAlmostEqual(summary1.mean, summarizer.mean)
        self.assertAlmostEqual(summary1.variance, summarizer.variance)
        self.assertAlmostEqual(summary1.std_deviation, summarizer.std_deviation)

    def test_update_iter_default_init(self):
        x1 = [0.1, -3.2, 3.5, 11.4, -2.3, 3.3, -5.6, 7.8, 2.9, -9.3]
        x2 = [11.1, -13.2, 23.5]
        x3 = [-0.1, 13.2, 3.25, 13.4, 2.3, 3.13, -15.6, -17.8, 2.19, 19.3]
        summarizer = SampleSummarizer()
        summarizer1 = SampleSummarizer()
        summarizer2 = SampleSummarizer()
        summarizer3 = SampleSummarizer()
        summarizer.update_samples(x1 + x2 + x3)
        summarizer1.update_samples(x1)
        summarizer2.update_samples(x2)
        summarizer3.update_samples(x3)

        summary1 = SampleSummary(
                sample_size = summarizer1.n,
                mean = summarizer1.mean,
                variance = summarizer1.variance)
        summary2 = SampleSummary(
                sample_size = summarizer2.n,
                mean = summarizer2.mean,
                variance = summarizer2.variance)
        summary3 = SampleSummary(
                sample_size = summarizer3.n,
                mean = summarizer3.mean,
                variance = summarizer3.variance)

        s = SampleSummary()
        s.update_iter([summary1, summary2, summary3])
        self.assertEqual(s.n, summarizer.n)
        self.assertAlmostEqual(s.mean, summarizer.mean)
        self.assertAlmostEqual(s.variance, summarizer.variance)
        self.assertAlmostEqual(s.std_deviation, summarizer.std_deviation)

        s = SampleSummary(sample_size=0, mean=22.1, variance=33.4)
        s.update_iter([summary1, summary2, summary3])
        self.assertEqual(s.n, summarizer.n)
        self.assertAlmostEqual(s.mean, summarizer.mean)
        self.assertAlmostEqual(s.variance, summarizer.variance)
        self.assertAlmostEqual(s.std_deviation, summarizer.std_deviation)

if __name__ == '__main__':
    unittest.main()

