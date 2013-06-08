#! /usr/bin/env python

import sys
import os
import math
import copy
from cStringIO import StringIO

from pymsbayes.fileio import process_file_arg
from pymsbayes.utils.parsing import parse_summary_file
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class SampleSummarizer(object):
    def __init__(self, name='sample summarizer'):
        self.name = name
        self._min = None
        self._max = None
        self._n = 0
        self._sum = 0
        self._sum_of_squares = 0
    
    def add_sample(self, x):
        self._n += 1
        self._sum += x
        self._sum_of_squares += x*x
        if not self._min:
            self._min = x
        elif x < self._min:
            self._min = x
        if not self._max:
            self._max = x
        elif x > self._max:
            self._max = x

    def update_samples(self, x_iter):
        for x in x_iter:
            self.add_sample(x)

    def _get_n(self):
        return self._n
    
    n = property(_get_n)

    def _get_min(self):
        return self._min

    def _get_max(self):
        return self._max

    minimum = property(_get_min)
    maximum = property(_get_max)
    
    def _get_mean(self):
        if self._n < 1:
            return None
        return float(self._sum)/self._n
    
    def _get_variance(self):
        if self._n < 1:
            return None
        if self._n == 1:
            return float('inf')
        return (self._sum_of_squares - (self._get_mean()*self._sum))/(self._n-1)

    def _get_std_dev(self):
        if self._n < 1:
            return None
        return math.sqrt(self._get_variance())

    def _get_pop_variance(self):
        if self._n < 1:
            return None
        return (self._sum_of_squares - (self._get_mean()*self._sum))/self._n

    mean = property(_get_mean)
    variance = property(_get_variance)
    pop_variance = property(_get_pop_variance)
    std_deviation = property(_get_std_dev)

    def __str__(self):
        s = StringIO()
        s.write('name = {0}\n'.format(self.name))
        s.write('sample size = {0}\n'.format(self._n))
        s.write('sum = {0}\n'.format(self._sum))
        s.write('sum of squares = {0}\n'.format(self._sum_of_squares))
        s.write('min = {0}\nmax = {1}\n'.format(self._min, self._max))
        s.write('mean = {0}\n'.format(self.mean))
        s.write('variance = {0}\n'.format(self.variance))
        s.write('pop variance = {0}\n'.format(self.pop_variance))
        return s.getvalue()

class SampleSummary(object):
    def __init__(self, sample_size = 0, mean = 0.0, variance = 0.0):
        self.n = sample_size
        self.mean = mean
        self.variance = variance
        if self.n < 1:
            self.n = 0
            self.mean = 0.0
            self.variance = 0.0

    def _get_std_dev(self):
        if self.n < 1:
            return None
        return math.sqrt(self.variance)
    std_deviation = property(_get_std_dev)

    def update(self, sample_summary):
        s1 = self
        s2 = sample_summary
        n = s1.n + s2.n
        mean = ((s1.n / float(n)) * s1.mean) + \
               ((s2.n / float(n)) * s2.mean)
        v = float(((s1.n**2) * s1.variance) + ((s2.n**2) * s2.variance) - \
                (s1.n * s1.variance) - (s1.n * s2.variance) - \
                (s2.n * s1.variance) - (s2.n * s2.variance) + \
                (s1.n * s2.n * s1.variance) + (s1.n * s2.n * s2.variance) + \
                (s1.n * s2.n * ((s1.mean - s2.mean)**2))) / \
                ((s1.n + s2.n - 1) * (s1.n + s2.n))
        self.n = n
        self.mean = mean
        self.variance = v

    def update_iter(self, sample_summaries):
        for ss in sample_summaries:
            self.update(ss)

class SampleSummaryCollection(object):
    def __init__(self, keys):
        self.keys = keys
        self.sample_sums = dict(zip(self.keys,
                [SampleSummary() for i in range(len(self.keys))]))

    @classmethod
    def get_from_summary_file(cls, summary_file_obj):
        stat_moments, header = parse_summary_file(summary_file_obj)
        assert(sorted(header) == sorted(stat_moments.keys()))
        ssc = cls(header)
        for k, v in stat_moments.iteritems():
            ssc.sample_sums[k] = SampleSummary(
                    sample_size = v['n'],
                    mean = v['mean'],
                    variance = v['std_deviation']**2)
        return ssc

    def update(self, sample_sum_collection):
        if sorted(self.keys) != sorted(sample_sum_collection.keys):
            raise Exception('SampleSummaryCollections have different keys')
        for stat_name, s_sum in self.sample_sums.iteritems():
            s_sum.update(sample_sum_collection.sample_sums[stat_name])

    def update_iter(self, sample_sum_collections):
        for ssc in sample_sum_collections:
            self.update(ssc)

    def get_copy(self):
        return copy.deepcopy(self)

    @classmethod
    def merge(cls, sample_sum_collections):
        ssc_iter = iter(sample_sum_collections)
        ssc = ssc_iter.next().get_copy()
        ssc.update_iter(ssc_iter)
        return ssc

    def write(self, file_obj):
        out, close = process_file_arg(file_obj, 'w')
        out.write('{0}\n'.format('\t'.join(self.keys)))
        out.write('{0}\n'.format('\t'.join(
                ['{0:.12f}'.format(self.sample_sums[
                        k].mean) for k in self.keys])))
        out.write('{0}\n'.format('\t'.join(
                ['{0:.12f}'.format(self.sample_sums[
                        k].std_deviation) for k in self.keys])))
        out.write('{0}\n'.format('\t'.join(
                ['{0}'.format(self.sample_sums[
                        k].n) for k in self.keys])))
        if close:
            out.close()

