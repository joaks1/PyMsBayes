#! /usr/bin/env python

import unittest
import os
import math
import types
from cStringIO import StringIO

from pymsbayes.fileio import process_file_arg
from pymsbayes.utils import GLOBAL_RNG
from pymsbayes.utils.stats import *
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.test.support import package_paths
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class SampleSummarizerTestCase(PyMsBayesTestCase):

    def test_init(self):
        ss = SampleSummarizer(tag='test')
        self.assertEqual(ss.tag, 'test')
        self.assertEqual(ss.minimum, None)
        self.assertEqual(ss.maximum, None)
        self.assertEqual(ss.mean, None)
        self.assertEqual(ss.variance, None)
        self.assertEqual(ss.std_deviation, None)
        self.assertEqual(ss.pop_variance, None)

    def test_add_one_sample(self):
        ss = SampleSummarizer(tag='test')
        ss.add_sample(1)
        self.assertEqual(ss.tag, 'test')
        self.assertEqual(ss.minimum, 1)
        self.assertEqual(ss.maximum, 1)
        self.assertApproxEqual(ss.mean, 1.0, 1e-9)
        self.assertEqual(ss.variance, float('inf'))
        self.assertEqual(ss.std_deviation, float('inf'))
        self.assertEqual(ss.pop_variance, 0)

        ss = SampleSummarizer(tag='test')
        ss.add_sample(3.45)
        self.assertEqual(ss.tag, 'test')
        self.assertEqual(ss.minimum, 3.45)
        self.assertEqual(ss.maximum, 3.45)
        self.assertApproxEqual(ss.mean, 3.45, 1e-9)
        self.assertEqual(ss.variance, float('inf'))
        self.assertEqual(ss.std_deviation, float('inf'))
        self.assertEqual(ss.pop_variance, 0)

    def test_update_samples(self):
        ss = SampleSummarizer(tag='test')
        ss.update_samples([1.0, 2.0, 3.0])
        self.assertEqual(ss.tag, 'test')
        self.assertEqual(ss.minimum, 1.0)
        self.assertEqual(ss.maximum, 3.0)
        self.assertApproxEqual(ss.mean, 2.0, 1e-9)
        self.assertApproxEqual(ss.variance, 1.0, 1e-9)
        self.assertEqual(ss.std_deviation, math.sqrt(1.0), 1e-9)
        self.assertApproxEqual(ss.pop_variance, 2/float(3), 1e-9)

    def test_init_with_samples(self):
        ss = SampleSummarizer([1.0, 2.0, 3.0])
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

class MedianTestCase(unittest.TestCase):
    def test_empty(self):
        samples = []
        self.assertRaises(ValueError, median, samples)
    
    def test_sample_size_1(self):
        samples = [1.3]
        med = median(samples)
        self.assertEqual(samples[0], med)

    def test_sample_size_even(self):
        samples = [1.1, 1.2, 1.3, 1.4]
        med = median(samples)
        self.assertAlmostEqual(med, 1.25)

    def test_sample_size_odd(self):
        samples = [1.1, 1.2, 1.3, 1.4, 1.5]
        med = median(samples)
        self.assertAlmostEqual(med, 1.3)

class ModeListTestCase(unittest.TestCase):
    def test_empty(self):
        samples = []
        self.assertRaises(ValueError, mode_list, samples)

    def test_ints(self):
        samples = [1,2,3,4,5]
        md = mode_list(samples)
        self.assertEqual(md, samples)

        samples = [1,2,2,3,4,5]
        md = mode_list(samples)
        self.assertEqual(md, [2])
        md = mode_list(samples, bin_width=None)
        self.assertEqual(md, [2])
        md = mode_list(samples, bin_width='a')
        self.assertEqual(md, [2])

        samples = [1,2,2,3,4,5,5]
        md = mode_list(samples)
        self.assertEqual(sorted(md), sorted([2, 5]))

    def test_strings(self):
        samples = ['a', 'b', 'b', 'c', 'd']
        md = mode_list(samples)
        self.assertEqual(md, ['b'])

    def test_floats_no_binning(self):
        samples = [1.1,2.1,2.1,3.1,4.1,5.1]
        md = mode_list(samples, bin_width=None)
        self.assertEqual(md, [2.1])
        md = mode_list(samples, bin_width='auto')
        self.assertNotEqual(md, [2.1])

    def test_floats(self):
        samples = [1.111, 1.112, 1.115, 1.16, 1.121]
        md = mode_list(samples, bin_width = 0.01, zero_value = 'b')
        self.assertEqual(sorted(md), sorted([(1.11, 1.12)]))

class IntervalTestCase(unittest.TestCase):
    def setUp(self):
        self.samples = [GLOBAL_RNG.normalvariate(0, 1) for i in range(100000)]
        self.exp_samples = [GLOBAL_RNG.expovariate(1) for i in range(100000)]

    def test_standard_normal_hpd(self):
        hpdi = get_hpd_interval(self.samples, 0.95)
        self.assertAlmostEqual(hpdi[0], -1.96, places=1)
        self.assertAlmostEqual(hpdi[1], 1.96, places=1)

    def test_standard_normal_quantile(self):
        quants = quantile_95(self.samples)
        q025 = quantile(self.samples, p=0.025)
        q975 = quantile(self.samples, p=0.975)
        self.assertAlmostEqual(q025, quants[0])
        self.assertAlmostEqual(q975, quants[1])
        self.assertAlmostEqual(quants[0], -1.96, places=1)
        self.assertAlmostEqual(quants[1], 1.96, places=1)

    def test_exp_hpd(self):
        hpdi = get_hpd_interval(self.exp_samples, 0.95)
        self.assertAlmostEqual(hpdi[0], 0.0, places=1)
        self.assertAlmostEqual(hpdi[1], 2.9957, places=1)

    def test_exp_quantile(self):
        quants = quantile_95(self.exp_samples)
        q025 = quantile(self.exp_samples, p=0.025)
        q975 = quantile(self.exp_samples, p=0.975)
        self.assertAlmostEqual(q025, quants[0])
        self.assertAlmostEqual(q975, quants[1])
        self.assertAlmostEqual(quants[0], 0.0253, places=1)
        self.assertAlmostEqual(quants[1], 3.6889, places=1)

class GetSummaryTestCase(unittest.TestCase):
    def setUp(self):
        self.samples = [GLOBAL_RNG.normalvariate(0, 1) for i in range(100000)]

    def test_standard_normal(self):
        d = get_summary(self.samples)
        self.assertEqual(d['n'], len(self.samples))
        self.assertEqual(d['range'][0], min(self.samples))
        self.assertEqual(d['range'][1], max(self.samples))
        self.assertAlmostEqual(d['mean'], 0.0, places=1)
        self.assertAlmostEqual(d['median'], 0.0, places=1)
        self.assertEqual(len(d['mode'][0]), 2)
        self.assertAlmostEqual(d['mode'][0][0], 0.0, places=0)
        self.assertAlmostEqual(d['mode'][0][1], 0.0, places=0)
        self.assertAlmostEqual(d['variance'], 1.0, places=1)
        self.assertAlmostEqual(d['qi_95'][0], -1.96, places=1)
        self.assertAlmostEqual(d['qi_95'][1], 1.96, places=1)
        self.assertAlmostEqual(d['hpdi_95'][0], -1.96, places=1)
        self.assertAlmostEqual(d['hpdi_95'][1], 1.96, places=1)
        
class IntegerPartitionTestCase(unittest.TestCase):

    def test_default_init_and_update(self):
        ip = IntegerPartition()
        self.assertFalse(ip._initialized)
        self.assertEqual(ip.n, 0)
        ip.update([0.1, 0.2, 0.2, 0.1, 0.2, 0.3])
        self.assertTrue(ip._initialized)
        self.assertEqual(ip.n, 1)
        self.assertEqual(ip.key, '3,2,1')
        self.assertEqual(ip.integer_partition, [3,2,1])
        self.assertEqual(ip._items, [(3, [0.2]), (2, [0.1]), (1, [0.3])])
        self.assertEqual(list(ip.iteritems()),
                [(3, [0.2]), (2, [0.1]), (1, [0.3])])
        self.assertEqual(list(ip.iter_item_summaries()),
                [(3, get_summary([0.2])),
                 (2, get_summary([0.1])),
                 (1, get_summary([0.3]))])
        self.assertEqual(ip.to_string(),
                ('3:0.2[&age_median=0.2,age_mean=0.2,age_n=1,'
                 'age_range={0.2,0.2},age_hpdi_95={0.2,0.2},'
                 'age_qi_95={0.2,0.2}],'
                 '2:0.1[&age_median=0.1,age_mean=0.1,age_n=1,'
                 'age_range={0.1,0.1},age_hpdi_95={0.1,0.1},'
                 'age_qi_95={0.1,0.1}],'
                 '1:0.3[&age_median=0.3,age_mean=0.3,age_n=1,'
                 'age_range={0.3,0.3},age_hpdi_95={0.3,0.3},'
                 'age_qi_95={0.3,0.3}]'))
        self.assertRaises(Exception, ip._initialize,
                [0.1, 0.2, 0.2, 0.1, 0.2, 0.3])

    def test_init_with_element_vector(self):
        ip = IntegerPartition([0.1, 0.2, 0.2, 0.1, 0.2, 0.3, 0.2])
        self.assertTrue(ip._initialized)
        self.assertEqual(ip.n, 1)
        self.assertEqual(ip.key, '4,2,1')
        self.assertEqual(ip.integer_partition, [4,2,1])
        self.assertEqual(ip._items, [(4, [0.2]), (2, [0.1]), (1, [0.3])])
        self.assertEqual(list(ip.iteritems()),
                [(4, [0.2]), (2, [0.1]), (1, [0.3])])
        self.assertEqual(list(ip.iter_item_summaries()),
                [(4, get_summary([0.2])),
                 (2, get_summary([0.1])),
                 (1, get_summary([0.3]))])
        self.assertEqual(ip.to_string(),
                ('4:0.2[&age_median=0.2,age_mean=0.2,age_n=1,'
                 'age_range={0.2,0.2},age_hpdi_95={0.2,0.2},'
                 'age_qi_95={0.2,0.2}],'
                 '2:0.1[&age_median=0.1,age_mean=0.1,age_n=1,'
                 'age_range={0.1,0.1},age_hpdi_95={0.1,0.1},'
                 'age_qi_95={0.1,0.1}],'
                 '1:0.3[&age_median=0.3,age_mean=0.3,age_n=1,'
                 'age_range={0.3,0.3},age_hpdi_95={0.3,0.3},'
                 'age_qi_95={0.3,0.3}]'))
        self.assertRaises(Exception, ip._initialize,
                [0.1, 0.2, 0.2, 0.1, 0.2, 0.3])

    def test_update_with_element_vectors(self):
        ip = IntegerPartition([0.1, 0.2, 0.2, 0.1, 0.2, 0.3, 0.2])
        self.assertEqual(ip.n, 1)
        ip.update([0.5, 0.4, 0.4, 0.1, 0.1, 0.1, 0.1])
        self.assertEqual(ip.n, 2)
        self.assertEqual(ip.key, '4,2,1')
        self.assertEqual(ip.integer_partition, [4,2,1])
        self.assertEqual(ip._items,
                [(4, [0.2, 0.1]), (2, [0.1, 0.4]), (1, [0.3, 0.5])])
        self.assertEqual(list(ip.iteritems()),
                [(4, [0.2, 0.1]), (2, [0.1, 0.4]), (1, [0.3, 0.5])])
        self.assertEqual(list(ip.iter_item_summaries()),
                [(4, get_summary([0.2, 0.1])),
                 (2, get_summary([0.1, 0.4])),
                 (1, get_summary([0.3, 0.5]))])
        self.assertEqual(ip.to_string(),
                ('4:0.15[&age_median=0.15,age_mean=0.15,age_n=2,'
                 'age_range={0.1,0.2},age_hpdi_95={0.1,0.2},'
                 'age_qi_95={0.1025,0.1975}],'
                 '2:0.25[&age_median=0.25,age_mean=0.25,age_n=2,'
                 'age_range={0.1,0.4},age_hpdi_95={0.1,0.4},'
                 'age_qi_95={0.1075,0.3925}],'
                 '1:0.4[&age_median=0.4,age_mean=0.4,age_n=2,'
                 'age_range={0.3,0.5},age_hpdi_95={0.3,0.5},'
                 'age_qi_95={0.305,0.495}]'))
        self.assertRaises(ValueError, ip.update,
                [0.4, 0.4, 0.4, 0.1, 0.1, 0.1, 0.1])

    def test_update_with_instance(self):
        ip = IntegerPartition([0.1, 0.2, 0.2, 0.1, 0.2, 0.3, 0.2])
        self.assertEqual(ip.n, 1)
        ip2 = IntegerPartition([0.5, 0.4, 0.4, 0.1, 0.1, 0.1, 0.1])
        ip.update(ip2)
        self.assertEqual(ip.n, 2)
        self.assertEqual(ip.key, '4,2,1')
        self.assertEqual(ip.integer_partition, [4,2,1])
        self.assertEqual(ip._items,
                [(4, [0.2, 0.1]), (2, [0.1, 0.4]), (1, [0.3, 0.5])])
        self.assertEqual(list(ip.iteritems()),
                [(4, [0.2, 0.1]), (2, [0.1, 0.4]), (1, [0.3, 0.5])])
        self.assertEqual(list(ip.iter_item_summaries()),
                [(4, get_summary([0.2, 0.1])),
                 (2, get_summary([0.1, 0.4])),
                 (1, get_summary([0.3, 0.5]))])
        self.assertEqual(ip.to_string(),
                ('4:0.15[&age_median=0.15,age_mean=0.15,age_n=2,'
                 'age_range={0.1,0.2},age_hpdi_95={0.1,0.2},'
                 'age_qi_95={0.1025,0.1975}],'
                 '2:0.25[&age_median=0.25,age_mean=0.25,age_n=2,'
                 'age_range={0.1,0.4},age_hpdi_95={0.1,0.4},'
                 'age_qi_95={0.1075,0.3925}],'
                 '1:0.4[&age_median=0.4,age_mean=0.4,age_n=2,'
                 'age_range={0.3,0.5},age_hpdi_95={0.3,0.5},'
                 'age_qi_95={0.305,0.495}]'))
        ip3 = IntegerPartition([0.4, 0.4, 0.4, 0.1, 0.1, 0.1, 0.1])
        self.assertRaises(ValueError, ip.update, ip3)

class IntegerPartitionCollectionTestCase(PyMsBayesTestCase):

    def test_default_init_and_add_element_vector(self):
        ipc = IntegerPartitionCollection()
        self.assertEqual(ipc.n, 0)
        self.assertEqual(ipc.integer_partitions, {})
        ipc.add([0.5, 0.4, 0.4, 0.1, 0.1, 0.1, 0.1])
        self.assertEqual(ipc.n, 1)
        self.assertEqual(ipc.keys(), ['4,2,1'])
        self.assertEqual(list(ipc.iterkeys()), ['4,2,1'])
        self.assertIsInstance(ipc.get('4,2,1'), IntegerPartition)
        self.assertEqual(ipc.get_count('4,2,1'), 1)
        self.assertEqual(ipc.get_frequency('4,2,1'), 1.0)
        s = StringIO()
        ipc.write_summary(s)

    def test_default_init_and_int_part(self):
        ipc = IntegerPartitionCollection()
        self.assertEqual(ipc.n, 0)
        self.assertEqual(ipc.integer_partitions, {})
        ip = IntegerPartition([0.5, 0.4, 0.4, 0.1, 0.1, 0.1, 0.1])
        ipc.add(ip)
        self.assertEqual(ipc.n, 1)
        self.assertEqual(ipc.keys(), ['4,2,1'])
        self.assertTrue(ipc.has_key('4,2,1'))
        self.assertIsInstance(ipc.iterkeys(), types.GeneratorType)
        self.assertEqual(list(ipc.iterkeys()), ['4,2,1'])
        self.assertIsInstance(ipc.get('4,2,1'), IntegerPartition)
        self.assertEqual(ipc.get_count('4,2,1'), 1)
        self.assertEqual(ipc.get_frequency('4,2,1'), 1.0)
        self.assertNotEqual(ipc.get('4,2,1'), ip)
        self.assertSameIntegerPartitions([ipc.get('4,2,1'), ip])
        self.assertSameIntegerPartitions(ipc.values() + [ip])
        self.assertIsInstance(ipc.itervalues(), types.GeneratorType)
        self.assertSameIntegerPartitions(list(ipc.itervalues()) + [ip])
        items = ipc.items()
        self.assertEqual(len(items), 1)
        self.assertEqual(items[0][0], ip.key)
        self.assertSameIntegerPartitions([items[0][1], ip])
        self.assertIsInstance(ipc.iteritems(), types.GeneratorType)
        items = list(ipc.iteritems())
        self.assertEqual(items[0][0], ip.key)
        self.assertSameIntegerPartitions([items[0][1], ip])

    def test_add(self):
        ipc = IntegerPartitionCollection()
        self.assertEqual(ipc.n, 0)
        self.assertEqual(ipc.integer_partitions, {})
        ip1 = IntegerPartition([0.5, 0.4, 0.4, 0.1, 0.1, 0.1, 0.1])
        ip2 = IntegerPartition([0.5, 0.4, 0.4, 0.1, 0.1, 0.1, 0.1])
        ip3 = IntegerPartition([0.1, 0.4, 0.4, 0.1, 0.1, 0.1, 0.1])
        ipc.add(ip1)
        ipc.add(ip2)
        ipc.add(ip3)
        ip1.update(ip2)
        self.assertEqual(ipc.n, 3)
        self.assertEqual(ipc.keys(), ['4,2,1', '5,2'])
        self.assertTrue(ipc.has_key('4,2,1'))
        self.assertTrue(ipc.has_key('5,2'))
        self.assertIsInstance(ipc.iterkeys(), types.GeneratorType)
        self.assertEqual(list(ipc.iterkeys()), ['4,2,1', '5,2'])
        self.assertIsInstance(ipc.get('4,2,1'), IntegerPartition)
        self.assertEqual(ipc.get_count('4,2,1'), 2)
        self.assertEqual(ipc.get_count('5,2'), 1)
        self.assertAlmostEqual(ipc.get_frequency('4,2,1'), 0.666666666)
        self.assertAlmostEqual(ipc.get_frequency('5,2'), 0.333333333)
        self.assertNotEqual(ipc.get('4,2,1'), ip1)
        self.assertNotEqual(ipc.get('4,2,1'), ip2)
        self.assertSameIntegerPartitions([ipc.get('4,2,1'), ip1])
        self.assertNotEqual(ipc.get('5,2'), ip3)
        self.assertSameIntegerPartitions([ipc.get('5,2'), ip3])
        self.assertNotEqual(ipc.values(), [ip1, ip3])
        self.assertSameIntegerPartitions([ipc.values()[0], ip1])
        self.assertSameIntegerPartitions([ipc.values()[1], ip3])
        self.assertIsInstance(ipc.itervalues(), types.GeneratorType)
        self.assertNotEqual(list(ipc.itervalues()), [ip1, ip3])
        self.assertSameIntegerPartitions([list(ipc.itervalues())[0], ip1])
        self.assertSameIntegerPartitions([list(ipc.itervalues())[1], ip3])
        self.assertNotEqual(ipc.items(), [(ip1.key, ip1), (ip3.key, ip3)])
        items = ipc.items()
        self.assertEqual(len(items), 2)
        self.assertEqual(items[0][0], ip1.key)
        self.assertSameIntegerPartitions([items[0][1], ip1])
        self.assertEqual(items[1][0], ip3.key)
        self.assertSameIntegerPartitions([items[1][1], ip3])
        self.assertIsInstance(ipc.iteritems(), types.GeneratorType)
        self.assertNotEqual(list(ipc.iteritems()), [(ip1.key, ip1), (ip3.key, ip3)])
        items = list(ipc.iteritems())
        self.assertEqual(len(items), 2)
        self.assertEqual(items[0][0], ip1.key)
        self.assertSameIntegerPartitions([items[0][1], ip1])
        self.assertEqual(items[1][0], ip3.key)

    def test_add_iter(self):
        ipc = IntegerPartitionCollection()
        self.assertEqual(ipc.n, 0)
        self.assertEqual(ipc.integer_partitions, {})
        ip1 = IntegerPartition([0.5, 0.4, 0.4, 0.1, 0.1, 0.1, 0.1])
        ip2 = IntegerPartition([0.5, 0.4, 0.4, 0.1, 0.1, 0.1, 0.1])
        ip3 = IntegerPartition([0.1, 0.4, 0.4, 0.1, 0.1, 0.1, 0.1])
        ipc.add_iter([ip1, ip2, ip3])
        ip1.update(ip2)
        self.assertEqual(ipc.n, 3)
        self.assertEqual(ipc.keys(), ['4,2,1', '5,2'])
        self.assertTrue(ipc.has_key('4,2,1'))
        self.assertTrue(ipc.has_key('5,2'))
        self.assertIsInstance(ipc.iterkeys(), types.GeneratorType)
        self.assertEqual(list(ipc.iterkeys()), ['4,2,1', '5,2'])
        self.assertIsInstance(ipc.get('4,2,1'), IntegerPartition)
        self.assertEqual(ipc.get_count('4,2,1'), 2)
        self.assertEqual(ipc.get_count('5,2'), 1)
        self.assertAlmostEqual(ipc.get_frequency('4,2,1'), 0.666666666)
        self.assertAlmostEqual(ipc.get_frequency('5,2'), 0.333333333)
        self.assertNotEqual(ipc.get('4,2,1'), ip1)
        self.assertNotEqual(ipc.get('4,2,1'), ip2)
        self.assertSameIntegerPartitions([ipc.get('4,2,1'), ip1])
        self.assertNotEqual(ipc.get('5,2'), ip3)
        self.assertSameIntegerPartitions([ipc.get('5,2'), ip3])
        self.assertNotEqual(ipc.values(), [ip1, ip3])
        self.assertSameIntegerPartitions([ipc.values()[0], ip1])
        self.assertSameIntegerPartitions([ipc.values()[1], ip3])
        self.assertIsInstance(ipc.itervalues(), types.GeneratorType)
        self.assertNotEqual(list(ipc.itervalues()), [ip1, ip3])
        self.assertSameIntegerPartitions([list(ipc.itervalues())[0], ip1])
        self.assertSameIntegerPartitions([list(ipc.itervalues())[1], ip3])
        self.assertNotEqual(ipc.items(), [(ip1.key, ip1), (ip3.key, ip3)])
        items = ipc.items()
        self.assertEqual(len(items), 2)
        self.assertEqual(items[0][0], ip1.key)
        self.assertSameIntegerPartitions([items[0][1], ip1])
        self.assertEqual(items[1][0], ip3.key)
        self.assertSameIntegerPartitions([items[1][1], ip3])
        self.assertIsInstance(ipc.iteritems(), types.GeneratorType)
        self.assertNotEqual(list(ipc.iteritems()), [(ip1.key, ip1), (ip3.key, ip3)])
        items = list(ipc.iteritems())
        self.assertEqual(len(items), 2)
        self.assertEqual(items[0][0], ip1.key)
        self.assertSameIntegerPartitions([items[0][1], ip1])
        self.assertEqual(items[1][0], ip3.key)

class GetCountsTestCase(unittest.TestCase):

    def test_get_counts(self):
        x = [0,0,0,1,1,1,1,2,3,4]
        expected = {0: 3, 1: 4, 2: 1, 3: 1, 4: 1}
        counts = get_counts(x)
        self.assertEqual(counts, expected)

class GetFreqsTestCase(unittest.TestCase):

    def test_get_counts(self):
        x = [0,0,0,1,1,1,1,2,3,4]
        expected = {0: 0.3, 1: 0.4, 2: 0.1, 3: 0.1, 4: 0.1}
        freqs = get_freqs(x)
        self.assertAlmostEqual(sum(freqs.values()), 1.0)
        for k, v in freqs.iteritems():
            self.assertAlmostEqual(v, expected[k])

class FreqLessThanTestCase(unittest.TestCase):

    def test_estimate_prob_zero(self):
        x = [0.0045, 0.00021, 0.00012, 0.009999, 0.001, 0.01, 0.010001, 0.9,
                0.09, 1.3]
        self.assertAlmostEqual(freq_less_than(x, 0.01), 0.5)
        self.assertAlmostEqual(freq_less_than(x, 2.0), 1.0)
        self.assertAlmostEqual(freq_less_than(x, 1.3), 0.9)

class SummarizeDiscreteParametersFromDensitiesTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.pdf_path = package_paths.data_path(
                'abctoolbox_posterior_density_file.txt')

    def tearDown(self):
        self.tear_down()

    def test_discrete_probs(self):
        probs = summarize_discrete_parameters_from_densities(self.pdf_path)
        self.assertEqual(len(probs), 1)
        self.assertEqual(probs.keys(), ['PRI.Psi'])
        self.assertEqual(sorted(probs['PRI.Psi'].keys()), sorted([1, 2]))
        self.assertAlmostEqual(probs['PRI.Psi'][1], 0.9999999950161, places=12)
        self.assertAlmostEqual(probs['PRI.Psi'][2], 4.9838793092e-09, places=12)

if __name__ == '__main__':
    unittest.main()

