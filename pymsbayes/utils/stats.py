#! /usr/bin/env python

import sys
import os
import math
import copy
import operator
import decimal
import fractions
import re
from cStringIO import StringIO

from pymsbayes.utils import GLOBAL_RNG
from pymsbayes.fileio import process_file_arg
from pymsbayes.utils import parsing
from pymsbayes.utils import functions
from pymsbayes.utils import probability
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

def freq_less_than(values, zero_threshold = 0.01):
    count = 0
    total = 0
    for v in values:
        if v < zero_threshold:
            count += 1
        total += 1
    if total < 1:
        return None
    return count / float(total)

def get_counts(elements):
    keys = sorted(set(elements))
    counts = dict(zip(keys, [0 for i in range(len(keys))]))
    for el in elements:
        counts[el] += 1
    return counts

def get_freqs(elements):
    counts = get_counts(elements)
    total = float(sum(counts.itervalues()))
    freqs = {}
    for k, v in counts.iteritems():
        freqs[k] = v / total
    return freqs

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

def spans_zero(samples):
    neg, pos = False, False
    for s in samples:
        if s < 0.0 and not neg:
            neg = True
        if s > 0.0 and not pos:
            pos = True
        if neg and pos:
            return True
    return False

def has_floats(samples):
    for s in samples:
        if isinstance(s, float) or isinstance(s, fractions.Fraction) or \
                isinstance(s, decimal.Decimal):
            return True
    return False

def get_bin_width(samples, algorithm = 'freedman-diaconis'):
    """
    Return "best" histogram bin width for samples using specified `algorithm`.
    Options for `algorithm` argument include:

    `algorithm = 'f'|'fd'|'freedman-diaconis'`
        Use the Freedman-Diaconis rule to calculate bin widths as
        2*IQR(samples) / n^(1/3), where n is the number of samples.

    `algorithm = 'c'|'custom'`
        Use custom (hacked) rescaled version of the Freedman-Diaconis rule to
        calculate. It returns the Freedman-Diaconis bin width multiplied by
        (n^(8/7) / (2 * n).
        This is very similar to the vanilla Freedman-Diaconis bin width at
        small sample sizes, but returns a more conservative (wider) bin width
        as sample sizes get large. I have found the F-D bin widths to be too
        narrow at large samples sizes (e.g., n > 10000) and this adjustment can
        allow more consistent estimation of the mode, yet be more precise than
        Sturges` and Doane's formulae.

    `algorithm = 's'|'sturges'`
        Use Sturges' formula to calculate the number of bins as k = LOG2(n) + 1
        where n is the number of samples, then return sample_range / k as the
        bin width.

    `algorithm = 'r'|'rice'`
        Use Rice Rule to estimate the number of bins as k = 2n^(1/3), where n
        is the number of samples, then return sample_range / k as the bin
        width.

    `algorithm = 'd'|'doane'`
        Use Doane's formula to estimate k:
            k = 1 + Log2(n) + Log2(1 + (|skewness|/sigma))
        where n is the number of samples, and
            sigma = sqrt((6(n-2)/((n+1)(n+3))))
        Then return sample_range / k as the bin width.
    """
    if not samples:
        return None
    if len(samples) < 2:
        return math.fabs(samples[0])
    a = algorithm.strip().lower()
    n = len(list(samples))
    if a in ['c', 'custom']:
        scaler = (n ** (float(8)/7)) / (2 * n)
        return scaler * get_bin_width(samples, 'freedman-diaconis')
    if a in ['f', 'fd', 'freedman-diaconis']:
        iqr = quantile(samples, 0.75) - quantile(samples, 0.25)
        return 2 * (float(iqr) / (n ** (float(1)/3)))
    elif a in ['s', 'sturges']:
        k = math.ceil(math.log(n, 2) + 1)
        return (max(samples) - min(samples)) / k
    elif a in ['d', 'doane']:
        if samples < 3:
            return get_bin_width(samples, 'freedman-diaconis')
        sigma = math.sqrt((6 * (n - 2)) / float((n + 1) * (n + 3)))
        ss = SampleSummarizer(samples)
        k = 1 + math.log(n, 2) + math.log((1 + (math.fabs(
                ss.skewness) / sigma)), 2)
        return (max(samples) - min(samples)) / k
    elif a in ['r', 'rice']:
        k = math.ceil(2 * (n ** (float(1)/3)))
        return (max(samples) - min(samples)) / k
    raise ValueError('unsupported `a` argument: {0!r}'.format(a))

def mode_list(samples, bin_width = 'auto', zero_value = 'boundary'):
    """
    Return a list of modes, or mode bins, from a list of values.

    Arguments include:

        `samples` is an iterable set of values, which can be integers, strings
        or floats.

        `bin_width` controls the behavior of the mode estimation, with the
        following options:
        
            `bin_width = 'a'|'auto'` - The default. The function automatically
            determines whether to treat the values as discrete or continuous by
            checking for floating point numbers in the sample. If there are no
            floats the samples are treated as discrete and a list of the most
            common values is returned. If there are floats, a bin width is
            determined by calling `get_bin_width(samples, algorithm='custom')`.
            The values are then binned using this bin width and the
            `zero_value` argument, and a list of tuples is returned, each tuple
            represents the lower and upper bounds of the most common bins.

            `bin_width = None|0` - The samples are treated as
            discrete and a list of the most common values is returned.

            `bin_width = <NUMBER OTHER THAN ZERO>` - The samples are treated as
            floats and are binned into categories of width `bin_width` to
            determine the mode.

            `bin_width =
                'c'|'custom'
                'f'|'fd'|'freedman-diaconis'|
                's'|'sturges'|
                'r'|'rice'|
                'd'|'doane'`
            The 'best' bin width is determined using the specified algorithm
            (see `get_bin_width` function for details regarding the algorithm
            options).
 
        `zero_value` zero always corresponds to a bin, and this option controls
        whether zero is at the center of a bin or at the edge. Options include:

            `zero_value = 'b'|'boundary'` - zero is a boundary between bins.
            In most cases choosing between 'b' and 'c' will be arbitrary, but
            if the samples are bounded by zero (i.e., the parameter is either
            strictly positive or negative, zero should be set as a boundary.

            `zero_value = 'c'|'center'` - zero is at the center of a bin.  If
            the samples are bounded by zero, use 'boundary' instead.  However,
            this option can be useful if the samples span zero and are
            suspected to be centered at zero.

    The function returns:

        If values are treated as discrete (i.e., `bin_width = None`), a list
        containing the most common values is returned. The list will contain
        multiple values if there is a tie.

        If values are treated as floats (i.e. `bin_width != None`), a list of
        tuples containing the lower and upper boundaries of the most common
        bins is returned. The list will contain multiple tuples each
        representing a most common bin, if there is a tie.

    Some examples:
        >>> from pymsbayes.utils.stats import mode_list
        >>> x = range(10) + [2]
        >>> mode_list(x)  # treat values as discrete by default
        [2]
        >>> x += [6]
        >>> mode_list(x)  # a tie!
        [2, 6]
        >>> x = ['a', 'b', 'b', 'c', 'c', 'b']
        >>> # strings work too when treated as discrete
        >>> mode_list(x)
        ['b']
        >>> import random
        >>> x = [random.Random().expovariate(1) for i in range(10000)]
        >>> # specify bin width for continuous values
        >>> mode_list(x, bin_width='auto')
        [(0.0, 0.10405355148832289)]
        >>> x = [random.Random().normalvariate(1, 1) for i in range(10000)]
        >>> mode_list(x, bin_width='auto')
        [(0.8910191831744725, 1.0183076379136828)]
        >>> x = [random.Random().normalvariate(0, 1) for i in range(10000)]
        >>> # zero is a bin boundary by default
        >>> mode_list(x, bin_width='auto') 
        [(-0.1263661814981197, 0.0)]
        >>> # specify zero_value as bin center to get mode that spans zero
        >>> mode_list(x, bin_width='auto', zero_value='center')
        [(-0.06318309074905985, 0.06318309074905985)]

    The beginnings of this function were based on the mode function in DendroPy
    (Copyright Jeet Sukumaran and Mark T. Holder; licensed under BSD License;
    http://pythonhosted.org/DendroPy/):

    Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
    for phylogenetic computing. Bioinformatics 26: 1569-1571.
    """
    if not samples:
        raise ValueError('empty samples')
    if len(list(samples)) == 1:
        return list(samples)
    zero_value = zero_value.strip().lower()
    discrete = False
    if not bin_width:
        discrete = True
    elif hasattr(bin_width, 'lower'):
        bin_width = bin_width.strip().lower()
        if bin_width in ['a', 'auto']:
            discrete = not has_floats(samples)
            if discrete:
                bin_width = None
            else:
                bin_width = 'c'
        else:
            discrete = False
        if not discrete:
            bin_width = get_bin_width(samples, bin_width)
            if bin_width == 0.0:
                bin_width = (max(samples) - min(samples)) / float(10)
                if bin_width == 0.0:
                    bin_width = 0.001
    if not discrete:
        bw = float(bin_width)
    counts = {}
    for s in samples:
        if discrete:
            index = s
        else:
            if zero_value in ['b', 'boundary']:
                index = int(math.floor(s / bw))
                bounds = (0.0, bw)
            elif zero_value in ['c', 'center']:
                index = int(math.floor((s / bw) + 0.5))
                bounds = ((bw / 2), (bw / 2))
            else:
                raise ValueError('unsupported `zero_value` argument: '
                        '{0!r}'.format(zero_value))
        counts[index] = counts.get(index, 0) + 1
    count_tups = sorted(counts.iteritems(), key = operator.itemgetter(1),
            reverse = True)
    max_count = count_tups[0][1]
    if discrete:
        return [val for val, cnt in count_tups if cnt >= max_count]
    return [((val * bin_width) - bounds[0], (val * bin_width) + bounds[1]) \
            for val, cnt in count_tups if cnt >= max_count]

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

def get_summary(samples, bin_width = 'auto'):
    """
    Return a dict of summaries calculated from the samples.

    The dict has the following items:
        'n': sample_size
        'mean': mean
        'median': median
        'modes': mode (tuple if binning)
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
            'modes': mode_list(samples, bin_width),
            'variance': ss.variance,
            'range': (min(samples), max(samples)),
            'hpdi_95': get_hpd_interval(samples, 0.95),
            'qi_95': quantile_95(samples)}
         
def summarize_discrete_parameters_from_densities(
        parameter_density_file,
        discrete_parameter_patterns = (parsing.MODEL_PATTERNS +
                parsing.PSI_PATTERNS +
                parsing.DIV_MODEL_PATTERNS),
        include_omega_summary = False,
        omega_threshold = 0.01,
        include_cv_summary = False,
        cv_threshold = 0.01):
    patterns = list(discrete_parameter_patterns)
    if include_omega_summary:
        patterns += parsing.OMEGA_PATTERNS
    if include_cv_summary:
        patterns += parsing.CV_PATTERNS
    densities = None
    for i, pd in enumerate(parsing.parameter_density_iter(
            parameter_density_file,
            patterns)):
        if not densities:
            densities = dict(zip(pd.keys(), [{} for i in range(len(pd))]))
        for k, val_dens_tup in pd.iteritems():
            if k == 'PRI.omega':
                if val_dens_tup[0] < omega_threshold:
                    densities[k][0] = densities[k].get(0, 0.0) + val_dens_tup[1]
                else:
                    densities[k][1] = densities[k].get(1, 0.0) + val_dens_tup[1]
            elif k == 'PRI.cv':
                if val_dens_tup[0] < cv_threshold:
                    densities[k][0] = densities[k].get(0, 0.0) + val_dens_tup[1]
                else:
                    densities[k][1] = densities[k].get(1, 0.0) + val_dens_tup[1]
            else:
                idx = int(round(val_dens_tup[0]))
                densities[k][idx] = densities[k].get(idx, 0.0) + val_dens_tup[1]
    probs = dict(zip(densities.keys(),
            [{} for i in range(len(densities))]))
    for pname, pdensities in densities.iteritems():
        total_dens = sum(pdensities.itervalues())
        for idx, density in pdensities.iteritems():
            probs[pname][idx] = density / total_dens
    return probs


class SampleSummarizer(object):
    count = 0
    def __init__(self, samples = None, tag = None):
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.tag = tag
        self._min = None
        self._max = None
        self._n = 0
        self._mean = 0.0
        self._sum_devs_2 = 0.0
        self._sum_devs_3 = 0.0
        self._sum_devs_4 = 0.0
        if samples:
            self.update_samples(samples)
    
    def add_sample(self, x):
        n = self._n + 1
        d = x - self._mean
        d_n = d / n
        d_n2 = d_n * d_n
        self._mean = self._mean + d_n
        first_term =  d * d_n * self._n
        self._sum_devs_4 += (first_term * d_n2 * ((n * n) - (3 * n) + 3)) + \
                (6 * d_n2 * self._sum_devs_2) - (4 * d_n * self._sum_devs_3)
        self._sum_devs_3 += (first_term * d_n * (n - 2)) - \
                (3 * d_n * self._sum_devs_2)
        self._sum_devs_2 += first_term
        self._n = n
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
        return self._mean
    
    def _get_variance(self):
        if self._n < 1:
            return None
        if self._n == 1:
            return float('inf')
        return (self._sum_devs_2 / (self._n - 1))

    def _get_std_dev(self):
        if self._n < 1:
            return None
        return math.sqrt(self._get_variance())

    def _get_pop_variance(self):
        if self._n < 1:
            return None
        return (self._sum_devs_2 / self._n)

    mean = property(_get_mean)
    variance = property(_get_variance)
    pop_variance = property(_get_pop_variance)
    std_deviation = property(_get_std_dev)

    def _get_skewness(self):
        return ((self._sum_devs_3 * math.sqrt(self._n)) / \
                (self._sum_devs_2 ** (float(3)/2)))
    def _get_kurtosis(self):
        return (((self._n * self._sum_devs_4) / (self._sum_devs_2 ** 2)) - 3)

    skewness = property(_get_skewness)
    kurtosis = property(_get_kurtosis)

    def __str__(self):
        s = StringIO()
        s.write('name = {0}\n'.format(self.name))
        s.write('sample size = {0}\n'.format(self._n))
        s.write('min = {0}\nmax = {1}\n'.format(self._min, self._max))
        s.write('mean = {0}\n'.format(self.mean))
        s.write('variance = {0}\n'.format(self.variance))
        s.write('pop variance = {0}\n'.format(self.pop_variance))
        s.write('skewness = {0}\n'.format(self.skewness))
        s.write('kurtosis = {0}\n'.format(self.kurtosis))
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
        stat_moments, header = parsing.parse_summary_file(summary_file_obj)
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

class ListConditionEvaluator(object):
    index_pattern = re.compile(r'(?P<index>\d+)')
    valid_relations = ['==', '!=', '<', '>', '<=', '>=', 'and', 'or', '&',
            '|', 'not']
    def __init__(self, expression_str,
            index_labels = None):
        relations = [s for s in re.split(r'[0-9()\s\t\n]',
                expression_str) if s != '']
        for r in relations:
            if not r in self.valid_relations:
                raise SyntaxError('expression {0!r} contains invalid relation '
                        '{1!r}'.format(expression_str, r))
        self.index_labels = None
        if index_labels:
            self.index_labels = list(index_labels)
        exp_str = ''.join(expression_str.split())
        self.indices = [int(i) for i in self.index_pattern.findall(exp_str)]
        exp = self.index_pattern.sub(' l[\g<index>] ', exp_str)
        self.expression = exp.strip()
        self.pretty_expression = self.expression
        pretty_exp = None
        if self.index_labels:
            pretty_exp = exp_str
            for i in sorted(self.indices, reverse = True):
                try:
                    pretty_exp = pretty_exp.replace(str(i), ' {0} '.format(
                            self.index_labels[i]))
                except IndexError as ie:
                    raise IndexError('expression {0!r} uses indices out of '
                            'range for labels {1}'.format(self.expression,
                                    self.index_labels))
        if pretty_exp:
            self.pretty_expression = pretty_exp.strip()

    def __str__(self):
        return self.expression

    def evaluate(self, values):
        l = list(values)
        try:
            ret = eval(self.expression)
        except IndexError as ie:
            raise IndexError('expression {0!r} uses indices out of range for '
                    '{1}'.format(self.expression, l))
        except SyntaxError as se:
            _LOG.error('ERROR: expression {0!r} has invalid syntax'.format(
                    self.expression))
            raise se
        return ret

class Partition(object):
    def __init__(self, element_vector = None):
        self._initialized = False
        self.values = {}
        self.partition = []
        self.key = None
        self.n = 0
        if element_vector:
            self._initialize(element_vector)

    def _initialize(self, element_vector):
        if self._initialized:
            raise Exception('cannot re-initialize Partition instance')
        if isinstance(element_vector, Partition):
            self.values = copy.deepcopy(element_vector.values)
            self.partition = copy.deepcopy(
                    element_vector.partition)
            self.key = copy.deepcopy(element_vector.key)
            self.n = element_vector.n
        else:
            self.partition, self.values = self.standardize(element_vector)
            self.key = ','.join([str(x) for x in self.partition])
            self.n = 1
        self._initialized = True

    def standardize(self, element_vector):
        el_map = {}
        next_idx = 0
        partition = []
        values = {}
        for i, el in enumerate(element_vector):
            if not el_map.has_key(el):
                el_map[el] = next_idx
                values[next_idx] = [el]
                next_idx += 1
            partition.append(el_map[el])
        return partition, values

    def itervalues(self):
        return ((k, self.values[k]) for k in sorted(self.values.iterkeys()))

    def iter_value_summaries(self):
        return ((k, get_summary(self.values[k])) for k in sorted(
                self.values.iterkeys()))

    def value_summary_string(self):
        elements = []
        for k, v in self.iter_value_summaries():
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

    def update(self, partition):
        if not self._initialized:
            self._initialize(partition)
            return
        if isinstance(partition, Partition):
            p = partition
        else:
            p = Partition(partition)
        if p.key != self.key:
            raise ValueError('partition passed to update method '
                    '({0}) does not match partition of current instance '
                    '({1})'.format(p.key, self.key))
        for k, values in p.values.iteritems():
            assert p.n == len(values)
            self.values[k].extend(values)
        self.n += p.n

    def get_element_vector(self, index = 0):
        if self.n < 1:
            return None
        return [self.values[i][index] for i in self.partition]

    def element_vector_iter(self):
        for i in range(self.n):
            yield self.get_element_vector(i)

    def get_condition_count(self, list_condition_evaluator_instance):
        count = 0
        for values in self.element_vector_iter():
            if list_condition_evaluator_instance.evaluate(values):
                count += 1
        return count

    def dirichlet_process_prior_probability(self, alpha, log = False):
        alpha = float(alpha)
        assert alpha > 0.0
        log_alpha = math.log(alpha)
        log_probs = [0.0]
        num_subsets = 1
        subset_counts = [1]
        elements = [0]
        for i in range(1, len(self.partition)):
            el = self.partition[i]
            if el not in elements:
                log_probs.append(log_alpha - math.log(alpha + i))
                elements.append(num_subsets)
                subset_counts.append(1)
                num_subsets += 1
                continue
            log_probs.append(math.log(subset_counts[el]) - math.log(alpha + i))
            subset_counts[el] += 1
        log_p = sum(log_probs)
        if log:
            return log_p
        return math.exp(log_p)

    def get_integer_partition(self):
        ip = IntegerPartition(self.get_element_vector(0))
        for i in range(1, self.n):
            ip.update(IntegerPartition(self.get_element_vector(i)))
        return ip

    def instance_from_base_distribution(self, base_distribution,
            partition = None,
            rng = None):
        if not rng:
            rng = GLOBAL_RNG
        elements = partition
        if partition == None:
            elements = self.partition
        elif hasattr(partition, 'partition'):
            elements = partition.partition
        value_map = {}
        for c in set(elements):
            value_map[c] = base_distribution.draw(rng)
        return Partition([value_map[el] for el in elements])

    def dirichlet_process_draw(self, alpha,
            num_elements = None,
            base_distribution = None,
            rng = None):
        if not num_elements:
            num_elements = len(self.partition)
        assert num_elements > 0
        if not rng:
            rng = GLOBAL_RNG
        if hasattr(alpha, 'draw'):
            alpha = alpha.draw(rng)
        alpha = float(alpha)
        assert alpha > 0.0
        elements = [0]
        subset_counts = [1]
        num_subsets = 1
        for i in range(1, num_elements):
            new_subset_prob = (alpha / (alpha + i))
            u = rng.random()
            u -= new_subset_prob
            if u < 0.0:
                elements.append(num_subsets)
                subset_counts.append(1)
                num_subsets += 1
                continue
            for j in range(num_subsets):
                subset_prob = (subset_counts[j] / (alpha + i))
                u -= subset_prob
                if u < 0.0:
                    elements.append(j)
                    subset_counts[j] += 1
                    break
            if u > 0.0:
                elements.append(num_subsets - 1)
                subset_counts[num_subsets - 1] += 1
        if base_distribution:
            return self.instance_from_base_distribution(
                    base_distribution = base_distribution,
                    partition = elements,
                    rng = rng)
        return Partition(elements)

    def dirichlet_process_draw_iter(self, alpha,
            num_samples = 1000,
            num_elements = None,
            base_distribution = None,
            rng = None):
        for i in range(num_samples):
            yield self.dirichlet_process_draw(
                    alpha = alpha,
                    num_elements = num_elements,
                    base_distribution = base_distribution,
                    rng = rng)

    def psi_multinomial_draw(self,
            num_elements = None,
            base_distribution = None,
            rng = None):
        if not num_elements:
            num_elements = len(self.partition)
        assert num_elements > 0
        if not rng:
            rng = GLOBAL_RNG
        psi_prior = probability.DiscreteUniformDistribution(1, num_elements)
        psi = psi_prior.draw(rng)
        subsets = list(range(psi))
        elements = [x for x in subsets]
        while len(elements) <= num_elements:
            elements.append(rng.choice(subsets))
        rng.shuffle(elements)
        if base_distribution:
            return self.instance_from_base_distribution(
                    base_distribution = base_distribution,
                    partition = elements,
                    rng = rng)
        return Partition(elements)

    def psi_multinomial_draw_iter(self,
            num_samples = 1000,
            num_elements = None,
            base_distribution = None,
            rng = None):
        for i in range(num_samples):
            yield self.psi_multinomial_draw(
                    num_elements = num_elements,
                    base_distribution = base_distribution,
                    rng = rng)

    def uniform_integer_partition_draw_iter(self,
            num_samples = 1000,
            num_elements = None,
            base_distribution = None,
            rng = None):
        ips = IntegerPartition.get_all_integer_partitions(
                num_elements)
        for i in range(num_samples):
            ip = ips.draw(rng)
            elements = ip.get_random_element_vector(index = 0,
                    rng = rng)
            if base_distribution:
                yield self.instance_from_base_distribution(
                        base_distribution = base_distribution,
                        partition = elements)
            else:
                yield Partition(elements)

    def number_of_partitions(self):
        """
        Total number of unique partitions (I.e., the Bell number).
        """
        n = len(self.partition)
        bell_num = 0
        for i in range(n):
            bell_num += stirling2(n, (i+1))
        return bell_num

    def number_of_partitions_into_k_subsets(self, k):
        """
        The number of unique ways to partition into `k` categories. (I.e., the
        Stirling number of 2nd kind).
        """
        return stirling2(len(self.partition), k)

    def get_dpp_expected_num_cats(self, concentration):
        """
        Calculates and returns the expectation of the number of categories
        under a Dirichlet process controlled by concentration parameter
        `concentration` for the number of elements in this `Partition`
        instance.

        Modified from `expNumTables` function of `util.h` from
        [`DPPDiv`](http://phylo.bio.ku.edu/content/tracy-heath-dppdiv) version
        1.0b (Copyright Tracy Heath, Mark Holder, and John Huelsenback;
        licensed under GPL v3;
        <http://phylo.bio.ku.edu/content/tracy-heath-dppdiv>).
        """
        expected_ncats = 0.0
        for i in range(1, len(self.partition) + 1):
            expected_ncats += (1.0 / (i - 1.0 + concentration))
        return (expected_ncats * concentration)

    def get_hyper_gamma_scale_from_shape_and_dpp_expected_ncats(self, shape, ncats):
        """
        Calculates and returns the scale parameter (given the shape parameter
        `shape`) of a gamma hyper prior on the concentration parameter of the
        Dirichlet process such that the prior expectation for the number of
        categories is `ncats`. The calculation is based on the number of
        elements in this `Partition` instance.
        """
        alpha = self.get_dpp_concentration(ncats)
        return (alpha / shape)

    def get_dpp_concentration(self, expected_num_cats,
            increment = 0.1,
            precision = 0.000001):
        """
        Calculates and returns the Dirichlet-process concentration parameter
        that has an expectation for the number of categories equal to
        `expected_num_cats`.

        Modified from `calculateFromPriorMean` function of `util.h` from
        [`DPPDiv`](http://phylo.bio.ku.edu/content/tracy-heath-dppdiv) version
        1.0b (Copyright Tracy Heath, Mark Holder, and John Huelsenback;
        licensed under GPL v3;
        <http://phylo.bio.ku.edu/content/tracy-heath-dppdiv>).
        """
        alpha = precision
        n = len(self.partition)
        current_exp_ncats = self.get_dpp_expected_num_cats(alpha)
        if expected_num_cats <= 1.0:
            expected_num_cats = 1.01
        increase = False
        if current_exp_ncats < expected_num_cats:
            increase = True
        while math.fabs(current_exp_ncats - expected_num_cats) > precision:
            if (current_exp_ncats < expected_num_cats) and increase:
                alpha += increment
            elif (current_exp_ncats > expected_num_cats) and (not increase):
                alpha -= increment
            elif (current_exp_ncats < expected_num_cats) and (not increase):
                increment /= 2.0
                increase = True
                alpha += increment
            else:
                increment /= 2.0
                increase = False
                alpha -= increment
            current_exp_ncats = self.get_dpp_expected_num_cats(alpha)
        return alpha


class PartitionCollection(object):
    def __init__(self, partitions = None):
        self.partitions = {}
        self.n = 0
        if partitions:
            self.add_iter(partitions)

    def add(self, partition):
        if isinstance(partition, Partition):
            p = partition
        else:
            p = Partition(partition)
        if self.partitions.has_key(p.key):
            self.partitions[p.key].update(p)
        else:
            self.partitions[p.key] = copy.deepcopy(p)
        self.n += p.n

    def add_iter(self, partitions):
        for p in partitions:
            self.add(p)

    def keys(self):
        return [k for k, v in self.iteritems()]

    def iterkeys(self):
        return (k for k in self.keys())

    def values(self):
        return sorted(self.partitions.values(),
                key = lambda x : x.n,
                reverse = True)

    def itervalues(self):
        return (v for v in self.values())

    def items(self):
        return sorted(self.partitions.iteritems(),
                key = lambda x : x[1].n,
                reverse = True)

    def iteritems(self):
        return ((k, v) for k, v in self.items())

    def has_key(self, k):
        return self.partitions.has_key(k)

    def get(self, k, d = None):
        return self.partitions.get(k, d)

    def get_count(self, key):
        if self.has_key(key):
            return self.get(key).n
        else:
            return 0

    def get_counts(self):
        counts = {}
        for k in self.iterkeys():
            counts[k] = self.get_count[k]
        return counts

    def get_frequency(self, key):
        return self.get_count(key) / float(self.n)

    def get_frequencies(self):
        freqs = {}
        for k in self.iterkeys():
            freqs[k] = self.get_frequency(k)
        return freqs

    def get_condition_count(self, list_condition_evaluator_instance):
        count = 0
        for p in self.partitions.itervalues():
            count += p.get_condition_count(list_condition_evaluator_instance)
        return count

    def get_condition_frequency(self, list_condition_evaluator_instance):
        return (self.get_condition_count(list_condition_evaluator_instance) /
                float(self.n))

    def prob_clustered(self, element_indices):
        element_indices = list(element_indices)
        prob_shared = 0.0
        for k, p in self.iteritems():
            div_indices = [d for i, d in enumerate(
                    p.partition) if i in element_indices]
            if len(set(div_indices)) == 1:
                prob_shared += self.get_frequency(k)
        return prob_shared
    
    def get_summary(self):
        stats = []
        for k in self.iterkeys():
            stats.append((k, {
                    'count': self.get_count(k),
                    'frequency': self.get_frequency(k),
                    'string': self.get(k).value_summary_string()}))
        return stats

    def write_summary(self, file_obj):
        out, close = process_file_arg(file_obj)
        out.write('count\tfrequency\tdivergence_model\tdiv_age_info\n')
        for k, v in self.iteritems():
            out.write('{0}\t{1}\t{2}\t{3}\n'.format(
                    v.n,
                    self.get_frequency(v.key),
                    v.key,
                    v.value_summary_string()))
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

    def uniform_prior_probability(self, num_elements = None):
        if not num_elements:
            num_elements = sum(self.integer_partition)
        assert num_elements > 0
        return float(1) / IntegerPartition.number_of_int_partitions(
                num_elements)

    def psi_uniform_prior_probability(self):
        psi = len(self.integer_partition)
        num_elements = sum(self.integer_partition)
        ip = IntegerPartition.number_of_int_partitions_by_k(num_elements)
        num_parts = ip[psi - 1]
        return float(1) / (num_elements * num_parts)

    def psi_multinomial_prior_probability(self):
        psi = len(self.integer_partition)
        num_elements = sum(self.integer_partition)
        psi_prob = float(1) / num_elements
        num_trials = num_elements - psi
        counts = [self.integer_partition[i] - 1 for i in range(psi)]
        assert num_trials == sum(counts)
        params = [float(1) / psi] * len(counts)
        m = proability.Multinomial(params)
        return psi_prob * m.pmf(counts)

    def get_random_element_vector(self, index = 0, rng = None):
        if self.n < 1:
            return None
        if not rng:
            rng = GLOBAL_RNG
        elements = []
        for k, v in self.iteritems():
            elements.extend([v[index]] * k)
        rng.shuffle(elements)
        return elements

    @classmethod
    def cumulative_number_of_int_partitions_by_k(cls, num_elements):
        assert num_elements > 0
        table = [[0] * (num_elements + 1) for i in range(num_elements + 1)]
        for i in range(1, num_elements + 1):
            table[0][i] = 1
        for i in range(1, num_elements + 1):
            for j in range(1, num_elements + 1):
                if i - j < 0:
                    table[i][j] = table[i][j-1]
                    continue
                table[i][j] = table[i][j-1] + table[i-j][j]
        return table[num_elements][1:]
    
    @classmethod
    def number_of_int_partitions_by_k(cls, num_elements):
        c = cls.cumulative_number_of_int_partitions_by_k(num_elements)
        return [c[0]] + [c[i] - c[i-1] for i in range(1, num_elements)]
    
    @classmethod
    def number_of_int_partitions(cls, num_elements):
        return cls.cumulative_number_of_int_partitions_by_k(num_elements)[-1]
    
    @classmethod
    def cumulative_frequency_of_int_partitions_by_k(cls, num_elements):
        c = cls.cumulative_number_of_int_partitions_by_k(num_elements)
        t = float(c[-1])
        return [x / t for x in c] 
    
    @classmethod
    def frequency_of_int_partitions_by_k(cls, num_elements):
        n = cls.number_of_int_partitions_by_k(num_elements)
        t = float(sum(n))
        return [x / t for x in n]

    @classmethod
    def get_all_integer_partitions(cls, num_elements):
        """
        A method for generating all partitions of an integer.
        
        This method generates all of the integer partitions of an integer in
        anti-lexicographic order. The code is based on Algorithm ZS1 from:
        
        Zoghbi, A., and I. Stojmenovic. 1998. Fast algorithms for generating
            generating integer partitions. Intern. J. Computer Math.
            70:319--332.
        
        @param num_elements
          The integer to partition.
        @return
          The method returns an IntegerPartitionCollection instance.
        """

        assert num_elements > 0
        partitions = IntegerPartitionCollection()

        ip = cls.number_of_int_partitions(num_elements)
        x = [1] * num_elements
        x[0] = num_elements
        m = 0
        h = 0
        partitions.add(cls([0] * num_elements))
        part_index = 0
        while x[0] != 1:
            part_index += 1
            if x[h] == 2:
                m += 1
                x[h] = 1
                h -= 1
            else:
                r = x[h] - 1
                t = m - h + 1
                x[h] = r
                while t >= r:
                    h += 1
                    x[h] = r
                    t -= r
                if t == 0:
                    m = h
                else:
                    m = h + 1
                    if t > 1:
                        h += 1
                        x[h] = t
            s = 0;
            part = []
            for i in range(num_elements):
                part.extend([i] * x[i])
                s += x[i];
                if s >= num_elements:
                    break
            partitions.add(cls(part))
        return partitions;
    

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

    def get_counts(self):
        counts = {}
        for k in self.iterkeys():
            counts[k] = self.get_count[k]
        return counts

    def get_frequency(self, key):
        return self.get_count(key) / float(self.n)

    def get_frequencies(self):
        freqs = {}
        for k in self.iterkeys():
            freqs[k] = self.get_frequency(k)
        return freqs
    
    def draw(self, rng = None):
        if not rng:
            rng = GLOBAL_RNG
        return rng.choice(self.values())
    
    def get_summary(self):
        stats = []
        for k in self.iterkeys():
            stats.append((k, {
                    'count': self.get_count(k),
                    'frequency': self.get_frequency(k),
                    'string': self.get(k).to_string()}))
        return stats

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

class ValidationProbabilities(object):
    def __init__(self, correct_to_est_prob_tups = None):
        self._bin_upper_limits = list(functions.frange(0.0, 1.0, 20,
                include_end_point = True))[1:]
        self.true_probs = None
        self.estimated_probs = None
        if correct_to_est_prob_tups:
            self.set_probs_from_correct_to_est_prob_tups(
                    correct_to_est_prob_tups)

    def _bin_correct_to_est_prob_tups(self, correct_prob_tups):
        bins = [[] for i in range(len(self._bin_upper_limits))]
        n = 0
        for (c, p) in correct_prob_tups:
            n += 1
            binned = False
            for i, l in enumerate(self._bin_upper_limits):
                if p < l:
                    bins[i].append((c, p))
                    binned = True
                    break
            if not binned:
                bins[i].append((c, p))
        length = 0
        for b in bins:
            length += len(b)
        assert length == n, '{0} tuples passed, but only {1} in bins'.format(
                n, length)
        assert len(bins) == len(self._bin_upper_limits)
        return bins

    def set_probs_from_correct_to_est_prob_tups(self, correct_prob_tups):
        bins = self._bin_correct_to_est_prob_tups(correct_prob_tups)
        half_bin_width = (1.0 / len(bins)) / 2.0
        self.estimated_probs = []
        self.true_probs = []
        for i, l in enumerate(self._bin_upper_limits):
            correct_vals = [c for (c, p) in bins[i]]
            if len(correct_vals) < 1:
                continue
            correct_prob = sum(correct_vals) / float(len(correct_vals))
            self.true_probs.append(correct_prob)
            self.estimated_probs.append(l - half_bin_width)

def mean_squared_error(x, y):
    if not len(x) == len(y):
        raise ValueError('x and y must be the same length')
    sse = 0.0
    for i in range(len(x)):
        sse += (x[i] - y[i]) ** 2
    return sse / len(x)

def root_mean_square_error(x, y):
    return mean_squared_error(x, y) ** 0.5

def stirling2_row():
    """
    Generate next row of Stirling numbers of the second kind
    """
    
    row = [1]
    yield row
    while True:
        mrow = [i*j for (i,j) in enumerate(row)]
        row = [i+j for (i,j) in zip(mrow + [0], [0] + row)]
        yield row

def stirling1_row():
    """
    Generate next row of Stirling numbers of the first kind
    """
    
    row = [1]
    yield row
    while True:
        n = len(row)
        mrow = [i*(n-1) for i in row]
        row = [i+j for (i,j) in zip(mrow + [0], [0] + row)]
        yield row
        
def stirling2(n, k):
    """
    Number of ways to partition n objects into k subsets
    """
    
    if n == 0:
        return 1
    if k == 1:
        return 1
    s = stirling2_row()
    for i in range(n + 1):
        r = s.next()
    return r[k]

def stirling1(n, k):
    """
    Number of ways to arrange n objects into k cycles
    """
    
    if n == 0:
        return 1
    if k == 1:
        return 1
    s = stirling1_row()
    for i in range(n + 1):
        r = s.next()
    return r[k]

