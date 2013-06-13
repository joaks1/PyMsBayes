#! /usr/bin/env python

import sys
import os
import math
import copy
import operator
from cStringIO import StringIO

from pymsbayes.fileio import process_file_arg
from pymsbayes.utils.parsing import parse_summary_file
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

def get_counts(elements):
    keys = sorted(set(elements))
    counts = dict(zip(keys, [0 for i in range(len(keys))]))
    for el in elements:
        counts[el] += 1
    return counts

def median(samples):
    """
    Return the median value from a list of samples.
    """
    s = sorted(samples)
    n = len(s)
    if n < 1:
        raise ValueError('empty samples')
    mdn = None
    if n % 2 == 0:
        mdn = (s[((n / 2) - 1)] + s[(n / 2)]) / 2
    else:
        mdn = s[((n - 1) / 2)]
    return mdn

def mode_list(samples, bin_width = None):
    """
    Return a list of modes from a list of values.

    If `bin_width` is None or zero, the mode will be calculated as if the
    samples are discrete (e.g., ints or strings). If a `bin_width` is provided,
    the samples are treated like floats and are binned into categories of width
    `bin_width`.

    Modified from DendroPy (Copyright Jeet Sukumaran and Mark T. Holder;
    licensed under BSD License; http://pythonhosted.org/DendroPy/):

    Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
    for phylogenetic computing. Bioinformatics 26: 1569-1571.
    """
    if not samples:
        raise ValueError('empty samples')
    counts = {}
    for s in samples:
        if bin_width:
            index = int(round(s / bin_width))
        else:
            index = s
        counts[index] = counts.get(index, 0) + 1
    count_tups = sorted(counts.iteritems(), key = operator.itemgetter(1),
            reverse = True)
    max_count = count_tups[0][1]
    if bin_width:
        return [(val * bin_width) for val, cnt in count_tups if cnt >= max_count]
    return [val for val, cnt in count_tups if cnt >= max_count]

def get_hpd_interval(samples, interval_prob = 0.95):
    """
    Return tuple of highest posterior density (HPD).

    The interval is estimated via a sliding window to find the narrowest
    interval that contains the specified proportion of samples.
    """
    if not samples:
        raise ValueError('empty samples')
    if interval_prob <= 0.0:
        raise ValueError('hpd probability interval is zero')
    samples = sorted(samples)
    tail_prob = 1 - interval_prob
    n = len(samples)
    num_exclude = int(round((n * tail_prob) - 0.5))
    if num_exclude < 1:
        num_exclude = 0
    widths = []
    # sliding window to find possible interval widths
    for i in range(num_exclude+1):
        lower = samples[i]
        upper = samples[(n - 1) - num_exclude + i]
        widths.append(upper - lower)
    min_width = min(widths)
    min_index = widths.index(min_width)
    return(samples[min_index], samples[(n - 1) - num_exclude + min_index])

def quantile(samples, p): 
    """
    Return quantile associated with probability `p`.

    Modified from code by Wai Yip Tung, licensed under PSF license and available
    here:
    http://code.activestate.com/recipes/511478-finding-the-percentile-of-the-values/
    """
    if not samples:
        raise ValueError('empty samples')
    s = sorted(samples)
    k = (len(s) - 1) * p
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return s[int(k)]
    d0 = s[int(f)] * (c - k)
    d1 = s[int(c)] * (k - f)
    return d0 + d1

def quantile_95(samples):
    """
    Return tuple of interval of 2.5% and 97.5% quantiles.
    """
    return (quantile(samples, 0.025), quantile(samples, 0.975))

def get_summary(samples, bin_width = None):
    """
    Return a dict of summaries calculated from the samples.

    The dict has the following items:
        'n': sample_size
        'mean': mean
        'median': median
        'variance': variance
        'range': range
        'hpdi_95': tuple of 95% highest posterior density interval
        'qi_95': tuple of 2.5% to 97.5% quantile interval
    """
    ss = SampleSummarizer()
    ss.update_samples(samples)
    return {'n': ss.n,
            'mean': ss.mean,
            'median': median(samples),
            'variance': ss.variance,
            'range': (min(samples), max(samples)),
            'hpdi_95': get_hpd_interval(samples, 0.95),
            'qi_95': quantile_95(samples)}
         

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

class IntegerPartition(object):
    def __init__(self, element_vector = None):
        self._initialized = False
        self._items = []
        self.integer_partition = []
        self.key = None
        self.n = 0
        if element_vector:
            self._initialize(element_vector)

    def iteritems(self):
        return ((k, v) for k, v in self._items)

    def iter_item_summaries(self):
        return ((k, get_summary(v)) for k, v in self._items)

    def to_string(self):
        elements = []
        for k, v in self.iter_item_summaries():
            int_age = ':'.join([str(k), str(v['median'])])
            comments = [('age_median', v['median']),
                        ('age_mean', v['mean']),
                        ('age_n', v['n']),
                        ('age_range', '{{{0},{1}}}'.format(*v['range'])),
                        ('age_hpdi_95', '{{{0},{1}}}'.format(*v['hpdi_95'])),
                        ('age_qi_95', '{{{0},{1}}}'.format(*v['qi_95']))]
            comment_str = ','.join(['='.join([x, str(y)]) for x, y in comments])
            elements.append('{0}[&{1}]'.format(int_age, comment_str))
        return ','.join(elements)

    def __str__(self):
        return self.to_string()

    def _initialize(self, element_vector):
        if self._initialized:
            raise Exception('cannot re-initialize IntegerPartition instance')
        if isinstance(element_vector, IntegerPartition):
            self._items = copy.deepcopy(element_vector._items)
            self.integer_partition = copy.deepcopy(
                    element_vector.integer_partition)
            self.key = copy.deepcopy(element_vector.key)
            self.n = element_vector.n
        else:
            counts = get_counts(element_vector)
            self._items = sorted([(v, [k]) for k, v in counts.iteritems()],
                    reverse = True)
            self.integer_partition = [k for k, v in self._items]
            self.key = ','.join([str(x) for x in self.integer_partition])
            self.n = 1
        self._initialized = True
    
    def update(self, integer_partition):
        if not self._initialized:
            self._initialize(integer_partition)
            return
        if isinstance(integer_partition, IntegerPartition):
            ip = integer_partition
        else:
            ip = IntegerPartition(integer_partition)
        if ip.key != self.key:
            raise ValueError('integer partition passed to update method '
                    '({0}) does not match partition of current instance '
                    '({1})'.format(ip.key, self.key))
        for i in range(len(ip._items)):
            assert ip.n == len(ip._items[i][1])
            self._items[i][1].extend(ip._items[i][1])
        self.n += ip.n

class IntegerPartitionCollection(object):
    def __init__(self, integer_partitions = None):
        self.integer_partitions = {}
        self.n = 0
        if integer_partitions:
            self.add_iter(integer_partitions)

    def add(self, integer_partition):
        if isinstance(integer_partition, IntegerPartition):
            ip = integer_partition
        else:
            ip = IntegerPartition(integer_partition)
        if self.integer_partitions.has_key(ip.key):
            self.integer_partitions[ip.key].update(ip)
        else:
            self.integer_partitions[ip.key] = copy.deepcopy(ip)
        self.n += ip.n

    def add_iter(self, integer_partitions):
        for ip in integer_partitions:
            self.add(ip)

    def keys(self):
        return [k for k, v in self.iteritems()]

    def iterkeys(self):
        return (k for k in self.keys())

    def values(self):
        return sorted(self.integer_partitions.values(),
                key = lambda x : x.n,
                reverse = True)

    def itervalues(self):
        return (v for v in self.values())

    def items(self):
        return sorted(self.integer_partitions.iteritems(),
                key = lambda x : x[1].n,
                reverse = True)

    def iteritems(self):
        return ((k, v) for k, v in self.items())

    def has_key(self, k):
        return self.integer_partitions.has_key(k)

    def get(self, k, d = None):
        return self.integer_partitions.get(k, d)

    def get_count(self, key):
        if self.has_key(key):
            return self.get(key).n
        else:
            return 0

    def get_frequency(self, key):
        return self.get_count(key) / float(self.n)

    def write_summary(self, file_obj):
        out, close = process_file_arg(file_obj)
        out.write('count\tfrequency\tdivergence_model\tdiv_model_with_age_info\n')
        for k, v in self.iteritems():
            out.write('{0}\t{1}\t{2}\t{3}\n'.format(
                    v.n,
                    self.get_frequency(v.key),
                    v.key,
                    v.to_string()))
        if close:
            out.close()

