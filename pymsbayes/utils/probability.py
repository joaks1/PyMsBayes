#! /usr/bin/env python

import sys
import os
import math

from pymsbayes.utils import GLOBAL_RNG
from pymsbayes.utils.oproperty import OProperty as property
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class Distribution(object):
    name = 'Abstract distribution object'

    def _get_min(self):
        raise NotImplementedError

    def _get_max(self):
        raise NotImplementedError

    minimum = property(_get_min)
    maximum = property(_get_max)

    def _get_mean(self):
        raise NotImplementedError

    mean = property(_get_mean)

    def _get_median(self):
        raise NotImplementedError

    median = property(_get_median)

    def _get_mode(self):
        raise NotImplementedError

    mode = property(_get_mode)

    def _get_variance(self):
        raise NotImplementedError

    variance = property(_get_variance)

    def _get_std_dev(self):
        try:
            sd = math.sqrt(self.variance)
        except:
            raise NotImplementedError
        return sd

    std_deviation = property(_get_std_dev)

class ContinuousUniformDistribution(Distribution):
    name = 'uniform'
    def __init__(self, minimum, maximum):
        if not maximum >= minimum:
            raise ValueError('max must be greater than min')
        self._min = float(minimum)
        self._max = float(maximum)

    def _get_min(self):
        return self._min

    def _get_max(self):
        return self._max

    def _get_mean(self):
        return float(self._min + self._max) / 2

    def _get_median(self):
        return self.mean

    def _get_variance(self):
        return (self._max - self._min)**2 / float(12)

    def probability_density(self, x):
        if self._min <= x <= self._max:
            try:
                return 1 / float(self._max - self._min)
            except ZeroDivisionError:
                return float('inf')
        else:
            return 0
    
    def draw(self, rng=None):
        if not rng:
            rng = GLOBAL_RNG
        return rng.uniform(self._min, self._max)

    def draw_iter(self, n, rng=None):
        if not rng:
            rng = GLOBAL_RNG
        for i in range(n):
            yield self.draw(rng)

    def __str__(self):
        return 'U({0},{1})'.format(self._min, self._max)
 
class BetaDistribution(Distribution):
    name = 'beta'
    def __init__(self, alpha, beta, multiplier=1):
        if not alpha > 0 or not beta > 0:
            raise ValueError('alpha and beta must be positive')
        self._a = alpha
        self._b = beta
        self._m = multiplier

    def _get_alpha(self):
        return self._a

    def _get_beta(self):
        return self._b

    alpha = property(_get_alpha)
    beta = property(_get_beta)

    def _get_min(self):
        return 0

    def _get_max(self):
        return self._m

    def _get_mean(self):
        return (self._a / float(self._a + self._b)) * self._m

    def _get_mode(self):
        m = (self._a - 1) / float(self._a + self._b - 2)
        return m * self._m

    def _get_variance(self):
        n = self._a * self._b
        d = ((self._a + self._b)**2) * (self._a + self._b + 1)
        return (n/float(d)) * (self._m**2)

    def draw(self, rng=None):
        if not rng:
            rng = GLOBAL_RNG
        return rng.betavariate(self._a, self._b)
        
    def draw_iter(self, n, rng=None):
        if not rng:
            rng = GLOBAL_RNG
        for i in range(n):
            yield self.draw(rng)

    def __str__(self):
        return 'Beta({0},{1})*{2}'.format(self._a, self._b, self._m)

class DiscreteUniformDistribution(Distribution):
    name = 'discrete-uniform'
    def __init__(self, minimum, maximum):
        if not maximum >= minimum:
            raise ValueError('max must be greater than min')
        if not isinstance(minimum, int) or not isinstance(maximum, int):
            raise ValueError('max and min must be integers')
        self._min = minimum
        self._max = maximum
        self._n = maximum - minimum + 1

    def _get_min(self):
        return self._min

    def _get_max(self):
        return self._max

    def _get_mean(self):
        return float(self._min + self._max) / 2

    def _get_median(self):
        return self.mean

    def _get_variance(self):
        return ((self._n**2) - 1) / float(12)

    def probability_mass(self, x):
        if not isinstance(x, int):
            raise ValueError('x must be an integer')
        if self._min <= x <= self._max:
            return 1 / float(self._n)
        else:
            return 0
    
    def draw(self, rng=None):
        if not rng:
            rng = GLOBAL_RNG
        return rng.randint(self._min, self._max)

    def draw_iter(self, n, rng=None):
        if not rng:
            rng = GLOBAL_RNG
        for i in range(n):
            yield self.draw(rng)

    def __str__(self):
        return 'U({0},{1},...,{2})'.format(self._min, self._min+1, self._max)

