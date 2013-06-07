#! /usr/bin/env python

import unittest
import os
import math

from pymsbayes.utils.stats import (SampleSummarizer, SampleSummary,
        merge_sample_summary_mappings)
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

class MergeSampleSummaryMappingsTestCase(PyMsBayesTestCase):

    def test_merge(self):
        samples1 = {'s1': [12.5, 64.3, 345.12],
                    's2': [3.53, 45.324, -23.455, 0.77],
                    's3': [2.234, -343.2, 23.21, 34.56, 33.2]}
        samples2 = {'s1': [122.5, 0.3, -345.12, 34.54, 0.001],
                    's2': [-0.5, 46.21, 23.4, 0.7799, 23.4],
                    's3': [23.021, 56.56, 33.2002]}
        samples3 = {'s1': [-0.99, 64.3335, 0.79, -12.3, 67.09, 56.4],
                    's2': [87.5, 99.65],
                    's3': [2.45, -4.5, 3.56, 33.2, 45.6, 67.8, 11.334]}
        samples = [samples1, samples2, samples3]
        summarizers_all = {'s1': SampleSummarizer(),
                           's2': SampleSummarizer(),
                           's3': SampleSummarizer()}
        summarizers1 = {'s1': SampleSummarizer(),
                        's2': SampleSummarizer(),
                        's3': SampleSummarizer()}
        summarizers2 = {'s1': SampleSummarizer(),
                        's2': SampleSummarizer(),
                        's3': SampleSummarizer()}
        summarizers3 = {'s1': SampleSummarizer(),
                        's2': SampleSummarizer(),
                        's3': SampleSummarizer()}
        for k, ss in summarizers_all.iteritems():
            samps = []
            for s in samples:
                samps += s[k]
            ss.update_samples(samps)
        for k, ss in summarizers1.iteritems():
            ss.update_samples(samples1[k])
        for k, ss in summarizers2.iteritems():
            ss.update_samples(samples2[k])
        for k, ss in summarizers3.iteritems():
            ss.update_samples(samples3[k])

        summaries1 = dict(zip(samples1.keys(), [None for i in range(len(samples1.keys()))]))
        for k in summaries1.iterkeys():
            summaries1[k] = SampleSummary(
                    sample_size = summarizers1[k].n,
                    mean = summarizers1[k].mean,
                    variance = summarizers1[k].variance)
        summaries2 = dict(zip(samples2.keys(), [None for i in range(len(samples2.keys()))]))
        for k in summaries2.iterkeys():
            summaries2[k] = SampleSummary(
                    sample_size = summarizers2[k].n,
                    mean = summarizers2[k].mean,
                    variance = summarizers2[k].variance)
        summaries3 = dict(zip(samples3.keys(), [None for i in range(len(samples3.keys()))]))
        for k in summaries3.iterkeys():
            summaries3[k] = SampleSummary(
                    sample_size = summarizers3[k].n,
                    mean = summarizers3[k].mean,
                    variance = summarizers3[k].variance)
        summaries_all = merge_sample_summary_mappings([summaries1,
                summaries2, summaries3])

        for k, ss in summaries_all.iteritems():
            self.assertEqual(ss.n, summarizers_all[k].n)
            self.assertAlmostEqual(ss.mean, summarizers_all[k].mean)
            self.assertAlmostEqual(ss.variance, summarizers_all[k].variance)
            self.assertAlmostEqual(ss.std_deviation, summarizers_all[k].std_deviation)
        self.assertNotEqual(summaries_all['s1'].n, summaries1['s1'].n)
        self.assertNotEqual(summaries_all['s1'].n, summaries2['s1'].n)
        self.assertNotEqual(summaries_all['s1'].n, summaries3['s1'].n)

if __name__ == '__main__':
    unittest.main()

