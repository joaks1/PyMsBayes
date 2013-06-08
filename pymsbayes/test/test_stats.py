#! /usr/bin/env python

import unittest
import os
import math
from cStringIO import StringIO

from pymsbayes.fileio import process_file_arg
from pymsbayes.utils.stats import (SampleSummarizer, SampleSummary,
        SampleSummaryCollection)
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

class SampleSummaryCollectionTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.summary_paths = []
        self.samples = [
                {'s1': [12.5, 64.3, 345.12, 23.42],
                 's2': [3.53, 45.324, -23.455, 0.77],
                 's3': [2.234, -343.2, 23.21, 33.2]},
                {'s1': [122.5, 0.3, -345.12, 34.54, 0.001],
                 's2': [-0.5, 46.21, 23.4, 0.7799, 23.4],
                 's3': [23.021, 56.56, 33.2002, 2.534, 4.22]},
                {'s1': [-0.99, 64.3335],
                 's2': [87.5, 99.65],
                 's3': [2.45, -4.5]},
                ]
        self.header = sorted(self.samples[0].keys())
        self.summarizers = [
                {'s1': SampleSummarizer(),
                 's2': SampleSummarizer(),
                 's3': SampleSummarizer()},
                {'s1': SampleSummarizer(),
                 's2': SampleSummarizer(),
                 's3': SampleSummarizer()},
                {'s1': SampleSummarizer(),
                 's2': SampleSummarizer(),
                 's3': SampleSummarizer()},
                ]
        for i in range(len(self.samples)):
            self.summary_paths.append(self.get_test_path(
                    prefix = 'summary-{0}-'.format(i)))
            for k, ss in self.summarizers[i].iteritems():
                ss.update_samples(self.samples[i][k])
        for i in range(len(self.samples)):
            out, close = process_file_arg(self.summary_paths[i], 'w')
            out.write('{0}\n'.format('\t'.join(self.header)))
            out.write('{0}\n'.format('\t'.join(['{0:.12}'.format(
                    self.summarizers[i][k].mean) for k in self.header])))
            out.write('{0}\n'.format('\t'.join(['{0:.12}'.format(
                    self.summarizers[i][k].std_deviation) for k in self.header])))
            out.write('{0}\n'.format('\t'.join(['{0}'.format(
                    self.summarizers[i][k].n) for k in self.header])))
            if close:
                out.close()

        self.summarizers_all = {'s1': SampleSummarizer(),
                                's2': SampleSummarizer(),
                                's3': SampleSummarizer()}
        for k, ss in self.summarizers_all.iteritems():
            samps = []
            for s in self.samples:
                samps += s[k]
            ss.update_samples(samps)

    def tearDown(self):
        self.tear_down()

    def test_get_from_summary_file(self):
        ssc = SampleSummaryCollection.get_from_summary_file(self.summary_paths[0])
        self.assertEqual(sorted(ssc.keys), sorted(self.summarizers[0].keys()))
        self.assertEqual(sorted(ssc.sample_sums.keys()), sorted(self.summarizers[0].keys()))
        for k, ss in self.summarizers[0].iteritems():
            self.assertEqual(ss.n, ssc.sample_sums[k].n)
            self.assertAlmostEqual(ss.mean, ssc.sample_sums[k].mean)
            self.assertAlmostEqual(ss.std_deviation, ssc.sample_sums[k].std_deviation)

    def test_get_from_summary_stream(self):
        stream, close = process_file_arg(self.summary_paths[0])
        ssc = SampleSummaryCollection.get_from_summary_file(stream)
        if close:
            stream.close()
        self.assertEqual(sorted(ssc.keys), sorted(self.summarizers[0].keys()))
        self.assertEqual(sorted(ssc.sample_sums.keys()), sorted(self.summarizers[0].keys()))
        for k, ss in self.summarizers[0].iteritems():
            self.assertEqual(ss.n, ssc.sample_sums[k].n)
            self.assertAlmostEqual(ss.mean, ssc.sample_sums[k].mean)
            self.assertAlmostEqual(ss.std_deviation, ssc.sample_sums[k].std_deviation)

    def test_get_from_cstring(self):
        stream, close = process_file_arg(self.summary_paths[0])
        cs = StringIO()
        cs.write(stream.read())
        cs.seek(0)
        if close:
            stream.close()
        ssc = SampleSummaryCollection.get_from_summary_file(cs)
        self.assertEqual(sorted(ssc.keys), sorted(self.summarizers[0].keys()))
        self.assertEqual(sorted(ssc.sample_sums.keys()), sorted(self.summarizers[0].keys()))
        for k, ss in self.summarizers[0].iteritems():
            self.assertEqual(ss.n, ssc.sample_sums[k].n)
            self.assertAlmostEqual(ss.mean, ssc.sample_sums[k].mean)
            self.assertAlmostEqual(ss.std_deviation, ssc.sample_sums[k].std_deviation)

    def test_init_update(self):
        ssc = SampleSummaryCollection(self.header)
        ssc2 = SampleSummaryCollection.get_from_summary_file(self.summary_paths[0])
        ssc.update(ssc2)
        self.assertEqual(sorted(ssc.keys), sorted(self.summarizers[0].keys()))
        self.assertEqual(sorted(ssc.sample_sums.keys()), sorted(self.summarizers[0].keys()))
        for k, ss in self.summarizers[0].iteritems():
            self.assertEqual(ss.n, ssc.sample_sums[k].n)
            self.assertAlmostEqual(ss.mean, ssc.sample_sums[k].mean)
            self.assertAlmostEqual(ss.std_deviation, ssc.sample_sums[k].std_deviation)

    def test_init_update_iter(self):
        ssc = SampleSummaryCollection(self.header)
        ssc_list = []
        for i in range(len(self.summary_paths)):
            ssc_list.append(SampleSummaryCollection.get_from_summary_file(
                    self.summary_paths[i]))
        ssc.update_iter(ssc_list)
        self.assertEqual(sorted(ssc.keys), sorted(self.summarizers_all.keys()))
        self.assertEqual(sorted(ssc.sample_sums.keys()), sorted(self.summarizers_all.keys()))
        for k, ss in self.summarizers_all.iteritems():
            self.assertEqual(ss.n, ssc.sample_sums[k].n)
            self.assertAlmostEqual(ss.mean, ssc.sample_sums[k].mean)
            self.assertAlmostEqual(ss.std_deviation, ssc.sample_sums[k].std_deviation)

    def test_merge(self):
        ssc_list = []
        for i in range(len(self.summary_paths)):
            ssc_list.append(SampleSummaryCollection.get_from_summary_file(
                    self.summary_paths[i]))
        ssc = SampleSummaryCollection.merge(ssc_list)
        self.assertEqual(sorted(ssc.keys), sorted(self.summarizers_all.keys()))
        self.assertEqual(sorted(ssc.sample_sums.keys()), sorted(self.summarizers_all.keys()))
        for k, ss in self.summarizers_all.iteritems():
            self.assertEqual(ss.n, ssc.sample_sums[k].n)
            self.assertAlmostEqual(ss.mean, ssc.sample_sums[k].mean)
            self.assertAlmostEqual(ss.std_deviation, ssc.sample_sums[k].std_deviation)
        for ssc2 in ssc_list:
            for k, ss in ssc2.sample_sums.iteritems():
                self.assertNotEqual(ss.n, ssc.sample_sums[k].n)

    def test_write(self):
        ssc_list = []
        for i in range(len(self.summary_paths)):
            ssc_list.append(SampleSummaryCollection.get_from_summary_file(
                    self.summary_paths[i]))
        ssc2 = SampleSummaryCollection.merge(ssc_list)
        ssc2.write(self.summary_paths[0])
        ssc = SampleSummaryCollection.get_from_summary_file(self.summary_paths[0])
        self.assertEqual(sorted(ssc.keys), sorted(self.summarizers_all.keys()))
        self.assertEqual(sorted(ssc.sample_sums.keys()), sorted(self.summarizers_all.keys()))
        for k, ss in self.summarizers_all.iteritems():
            self.assertEqual(ss.n, ssc.sample_sums[k].n)
            self.assertAlmostEqual(ss.mean, ssc.sample_sums[k].mean)
            self.assertAlmostEqual(ss.std_deviation, ssc.sample_sums[k].std_deviation)
        for ssc2 in ssc_list:
            for k, ss in ssc2.sample_sums.iteritems():
                self.assertNotEqual(ss.n, ssc.sample_sums[k].n)

if __name__ == '__main__':
    unittest.main()

