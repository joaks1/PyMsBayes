#! /usr/bin/env python

import unittest
import os
import sys
import math
import multiprocessing

from pymsbayes.workers import *
from pymsbayes.manager import *
from pymsbayes.utils.probability import *
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.test import TestLevel, test_enabled
from pymsbayes.utils.parsing import parse_parameters
from pymsbayes.utils.functions import long_division
from pymsbayes.utils.stats import mode_list
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

SAMPLE_TABLE = """BEGIN SAMPLE_TBL
pair0	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair0.locus0.fasta
pair0	locus1	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair0.locus1.fasta
pair1	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair2.locus0.fasta
pair2	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair3.locus0.fasta
pair2	locus1	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair3.locus1.fasta
pair3	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair3.locus0.fasta
END SAMPLE_TBL
"""

class PriorTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = self.get_test_path(prefix='test-config-')
        self.prior_path = self.get_test_path(prefix='test-prior-')
        self.preamble = ''
        self.sample_table = ''
        self.samples = None

    def tearDown(self):
        self.tear_down()

    def write_cfg(self):
        with open(self.cfg_path, 'w') as out:
            if self.sample_table:
                out.write('{0}\n{1}\n'.format(self.preamble, self.sample_table))
            else:
                out.write('{0}\n{1}\n'.format(self.preamble, SAMPLE_TABLE))

    def generate_prior(self, sample_size=1000, batch_size=250, np=4):
        self.write_cfg()
        num_batches, remainder = long_division(sample_size, batch_size)
        workers = []
        for i in range(num_batches):
            w = MsBayesWorker(
                    temp_fs = self.temp_fs,
                    sample_size = batch_size,
                    config_path = self.cfg_path,
                    schema = 'abctoolbox')
            WORK_FORCE.put(w)
            workers.append(w)
        if remainder > 0:
            w = MsBayesWorker(
                    temp_fs = self.temp_fs,
                    sample_size = remainder,
                    config_path = self.cfg_path,
                    schema = 'abctoolbox')
            WORK_FORCE.put(w)
            workers.append(w)

        result_q = multiprocessing.Queue()
        managers = []
        for i in range(np):
            m = Manager(work_queue = WORK_FORCE,
                    result_queue = result_q)
            m.start()
            managers.append(m)
        for i in range(len(workers)):
            workers[i] = result_q.get()
        for m in managers:
            m.join()
        assert WORK_FORCE.empty()
        assert result_q.empty()
        
        merge_prior_files([w.prior_path for w in workers], self.prior_path)
        self.samples = parse_parameters(self.prior_path, include_thetas=True)
        self.samples['unique_taus'] = []
        self.samples['all_a_thetas'] = []
        self.samples['all_d1_thetas'] = []
        self.samples['all_d2_thetas'] = []
        self.samples['mean_d_thetas'] = []
        for tau_list in self.samples['taus']:
            self.samples['unique_taus'].extend(list(set(tau_list)))
        for lst in self.samples['a_thetas']:
            self.samples['all_a_thetas'].extend(lst)
        for i in range(len(self.samples['d1_thetas'])):
            self.samples['all_d1_thetas'].extend(self.samples['d1_thetas'][i])
            self.samples['all_d2_thetas'].extend(self.samples['d2_thetas'][i])
            for j in range(len(self.samples['d1_thetas'][i])):
                mean_d_theta = (self.samples['d1_thetas'][i][j] + \
                        self.samples['d2_thetas'][i][j]) / 2
                self.samples['mean_d_thetas'].append(mean_d_theta)

    def test_old_unconstrained(self):
        if not test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        self.preamble = """
upperTheta = 0.1
lowerTheta = 0.0001
upperTau = 10.0
lowerTau = 5.0
numTauClasses = 0
upperMig = 0
upperRec = 0
upperAncPopSize = 1.0
constrain = 0
subParamConstrain = 111111111
"""

        self.generate_prior(sample_size=1000, batch_size=250, np=4)
        self.assertFalse(self.samples == None)
        theta = ContinuousUniformDistribution(0.0001, 0.1)
        tau = ContinuousUniformDistribution(5.0, 10.0)
        a_theta_ss = SampleSummarizer(self.samples['all_a_thetas'])
        d_theta_ss = SampleSummarizer(self.samples['mean_d_thetas'])
        tau_ss = SampleSummarizer(self.samples['unique_taus'])
        self.assertAlmostEqual(tau_ss.mean, tau.mean, places=1)
        self.assertAlmostEqual(tau_ss.variance, tau.variance, places=0)
        self.assertAlmostEqual(a_theta_ss.mean, theta.mean, places=2)
        self.assertAlmostEqual(d_theta_ss.mean, theta.mean, places=2)
        self.assertAlmostEqual(a_theta_ss.variance, theta.variance, places=4)
        self.assertAlmostEqual(d_theta_ss.variance, theta.variance, places=4)
        psi_freqs = get_freqs(self.samples['psi'])
        for psi, freq in psi_freqs.iteritems():
            self.assertAlmostEqual(freq, 0.25, places=1)

    def test_old_constrained(self):
        if not test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        self.preamble = """
upperTheta = 0.1
lowerTheta = 0.0001
upperTau = 10.0
lowerTau = 5.0
numTauClasses = 4
upperMig = 0
upperRec = 0
upperAncPopSize = 1.0
constrain = 0
subParamConstrain = 111111111
"""

        self.generate_prior(sample_size=1000, batch_size=250, np=4)
        self.assertFalse(self.samples == None)
        theta = ContinuousUniformDistribution(0.0001, 0.1)
        tau = ContinuousUniformDistribution(5.0, 10.0)
        a_theta_ss = SampleSummarizer(self.samples['all_a_thetas'])
        d_theta_ss = SampleSummarizer(self.samples['mean_d_thetas'])
        tau_ss = SampleSummarizer(self.samples['unique_taus'])
        self.assertAlmostEqual(tau_ss.mean, tau.mean, places=1)
        self.assertAlmostEqual(tau_ss.variance, tau.variance, places=0)
        self.assertAlmostEqual(a_theta_ss.mean, theta.mean, places=2)
        self.assertAlmostEqual(d_theta_ss.mean, theta.mean, places=2)
        self.assertAlmostEqual(a_theta_ss.variance, theta.variance, places=4)
        self.assertAlmostEqual(d_theta_ss.variance, theta.variance, places=4)
        self.assertEqual(self.samples['psi'], [4]*1000)
        for t in self.samples['taus']:
            self.assertEqual(len(set(t)), 4)

    def test_new_uniform_theta000(self):
        if not test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        self.preamble = """
concentrationShape = 0
concentrationScale = 0
thetaShape = 1
thetaScale = 0.001
ancestralThetaShape = 0
ancestralThetaScale = 0
thetaParameters = 000
tauShape = 1.0
tauScale = 2.0
bottleProportionShapeA = 0
bottleProportionShapeB = 0
bottleProportionShared = 1
numTauClasses = 0
migrationShape = 0
migrationScale = 0
recombinationShape = 0
recombinationScale = 0
constrain = 0
subParamConstrain = 111111111
"""

        self.generate_prior(sample_size=1000, batch_size=250, np=4)
        self.assertFalse(self.samples == None)
        theta = GammaDistribution(1, 0.001)
        tau = GammaDistribution(1.0, 2.0)
        for i in range(1000):
            for j in range(4):
                self.assertEqual(self.samples['a_thetas'][i][j],
                        self.samples['d1_thetas'][i][j])
                self.assertEqual(self.samples['a_thetas'][i][j],
                        self.samples['d2_thetas'][i][j])
        a_theta_ss = SampleSummarizer(self.samples['all_a_thetas'])
        tau_ss = SampleSummarizer(self.samples['unique_taus'])
        self.assertAlmostEqual(tau_ss.mean, tau.mean, places=0)
        self.assertAlmostEqual(tau_ss.variance, tau.variance, places=0)
        self.assertAlmostEqual(a_theta_ss.mean, theta.mean, places=2)
        self.assertAlmostEqual(a_theta_ss.variance, theta.variance, places=2)
        psi_freqs = get_freqs(self.samples['psi'])
        for psi, freq in psi_freqs.iteritems():
            if psi == 2:
                self.assertAlmostEqual(freq, 2/float(5), places=1)
            else:
                self.assertAlmostEqual(freq, 1/float(5), places=1)

    def test_new_uniform_theta001(self):
        if not test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        self.preamble = """
concentrationShape = 0
concentrationScale = 0
thetaShape = 1
thetaScale = 0.001
ancestralThetaShape = 1
ancestralThetaScale = 0.1
thetaParameters = 001
tauShape = 1.0
tauScale = 2.0
bottleProportionShapeA = 0
bottleProportionShapeB = 0
bottleProportionShared = 1
numTauClasses = 0
migrationShape = 0
migrationScale = 0
recombinationShape = 0
recombinationScale = 0
constrain = 0
subParamConstrain = 111111111
"""

        self.generate_prior(sample_size=1000, batch_size=250, np=4)
        self.assertFalse(self.samples == None)
        theta = GammaDistribution(1, 0.001)
        a_theta = GammaDistribution(1, 0.1)
        tau = GammaDistribution(1.0, 2.0)
        not_equal_failures = 0
        for i in range(1000):
            for j in range(4):
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d1_thetas'][i][j]:
                    not_equal_failures += 1
                self.assertEqual(self.samples['d1_thetas'][i][j],
                        self.samples['d2_thetas'][i][j])
        self.assertTrue(not_equal_failures < 10)
        a_theta_ss = SampleSummarizer(self.samples['all_a_thetas'])
        d1_theta_ss = SampleSummarizer(self.samples['all_d1_thetas'])
        tau_ss = SampleSummarizer(self.samples['unique_taus'])
        self.assertAlmostEqual(tau_ss.mean, tau.mean, places=0)
        self.assertAlmostEqual(tau_ss.variance, tau.variance, places=0)
        self.assertAlmostEqual(a_theta_ss.mean, a_theta.mean, places=2)
        self.assertAlmostEqual(a_theta_ss.variance, a_theta.variance, places=2)
        self.assertAlmostEqual(d1_theta_ss.mean, theta.mean, places=2)
        self.assertAlmostEqual(d1_theta_ss.variance, theta.variance, places=2)
        psi_freqs = get_freqs(self.samples['psi'])
        for psi, freq in psi_freqs.iteritems():
            if psi == 2:
                self.assertAlmostEqual(freq, 2/float(5), places=1)
            else:
                self.assertAlmostEqual(freq, 1/float(5), places=1)

    def test_new_uniform_theta011(self):
        if not test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        self.preamble = """
concentrationShape = 0
concentrationScale = 0
thetaShape = 4
thetaScale = 0.005
ancestralThetaShape = 1
ancestralThetaScale = 0.1
thetaParameters = 011
tauShape = 0.5
tauScale = 1.0
bottleProportionShapeA = 0
bottleProportionShapeB = 0
bottleProportionShared = 1
numTauClasses = 0
migrationShape = 0
migrationScale = 0
recombinationShape = 0
recombinationScale = 0
constrain = 0
subParamConstrain = 111111111
"""

        self.generate_prior(sample_size=1000, batch_size=250, np=4)
        self.assertFalse(self.samples == None)
        theta = GammaDistribution(4, 0.005)
        a_theta = theta
        tau = GammaDistribution(0.5, 1.0)
        not_equal_failures = 0
        for i in range(1000):
            for j in range(4):
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d1_thetas'][i][j]:
                    not_equal_failures += 1
                self.assertEqual(self.samples['a_thetas'][i][j],
                        self.samples['d2_thetas'][i][j])
        self.assertTrue(not_equal_failures < 10)
        a_theta_ss = SampleSummarizer(self.samples['all_a_thetas'])
        d1_theta_ss = SampleSummarizer(self.samples['all_d1_thetas'])
        tau_ss = SampleSummarizer(self.samples['unique_taus'])
        self.assertAlmostEqual(tau_ss.mean, tau.mean, places=1)
        self.assertAlmostEqual(tau_ss.variance, tau.variance, places=0)
        self.assertAlmostEqual(a_theta_ss.mean, a_theta.mean, places=2)
        self.assertAlmostEqual(a_theta_ss.variance, a_theta.variance, places=2)
        self.assertAlmostEqual(d1_theta_ss.mean, theta.mean, places=2)
        self.assertAlmostEqual(d1_theta_ss.variance, theta.variance, places=2)
        psi_freqs = get_freqs(self.samples['psi'])
        for psi, freq in psi_freqs.iteritems():
            if psi == 2:
                self.assertAlmostEqual(freq, 2/float(5), places=1)
            else:
                self.assertAlmostEqual(freq, 1/float(5), places=1)

    def test_new_uniform_theta010(self):
        if not test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        self.preamble = """
concentrationShape = 0
concentrationScale = 0
thetaShape = 1
thetaScale = 0.005
ancestralThetaShape = 1
ancestralThetaScale = 0.1
thetaParameters = 010
tauShape = 1.0
tauScale = 2.0
bottleProportionShapeA = 0
bottleProportionShapeB = 0
bottleProportionShared = 1
numTauClasses = 0
migrationShape = 0
migrationScale = 0
recombinationShape = 0
recombinationScale = 0
constrain = 0
subParamConstrain = 111111111
"""

        self.generate_prior(sample_size=1000, batch_size=250, np=4)
        self.assertFalse(self.samples == None)
        theta = GammaDistribution(1, 0.005)
        a_theta = theta
        tau = GammaDistribution(1.0, 2.0)
        not_equal_failures = 0
        for i in range(1000):
            for j in range(4):
                self.assertEqual(self.samples['a_thetas'][i][j],
                        self.samples['d1_thetas'][i][j])
                if self.samples['d1_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
        self.assertTrue(not_equal_failures < 10)
        a_theta_ss = SampleSummarizer(self.samples['all_a_thetas'])
        d2_theta_ss = SampleSummarizer(self.samples['all_d2_thetas'])
        tau_ss = SampleSummarizer(self.samples['unique_taus'])
        self.assertAlmostEqual(tau_ss.mean, tau.mean, places=0)
        self.assertAlmostEqual(tau_ss.variance, tau.variance, places=0)
        self.assertAlmostEqual(a_theta_ss.mean, a_theta.mean, places=2)
        self.assertAlmostEqual(a_theta_ss.variance, a_theta.variance, places=2)
        self.assertAlmostEqual(d2_theta_ss.mean, theta.mean, places=2)
        self.assertAlmostEqual(d2_theta_ss.variance, theta.variance, places=2)
        psi_freqs = get_freqs(self.samples['psi'])
        for psi, freq in psi_freqs.iteritems():
            if psi == 2:
                self.assertAlmostEqual(freq, 2/float(5), places=1)
            else:
                self.assertAlmostEqual(freq, 1/float(5), places=1)

    def test_new_uniform_theta012(self):
        if not test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        self.preamble = """
concentrationShape = 0
concentrationScale = 0
thetaShape = 0.5
thetaScale = 0.01
ancestralThetaShape = 1
ancestralThetaScale = 0.02
thetaParameters = 012
tauShape = 4.0
tauScale = 3.0
bottleProportionShapeA = 0
bottleProportionShapeB = 0
bottleProportionShared = 1
numTauClasses = 0
migrationShape = 0
migrationScale = 0
recombinationShape = 0
recombinationScale = 0
constrain = 0
subParamConstrain = 111111111
"""

        self.generate_prior(sample_size=1000, batch_size=250, np=4)
        self.assertFalse(self.samples == None)
        a_theta = GammaDistribution(1, 0.02)
        theta = GammaDistribution(0.5, 0.01)
        tau = GammaDistribution(4.0, 3.0)
        not_equal_failures = 0
        for i in range(1000):
            for j in range(4):
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d1_thetas'][i][j]:
                    not_equal_failures += 1
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
                if self.samples['d1_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
        self.assertTrue(not_equal_failures < 10)
        d_theta_ss = SampleSummarizer(self.samples['all_d1_thetas'] + \
                self.samples['all_d2_thetas'])
        a_theta_ss = SampleSummarizer(self.samples['all_a_thetas'])
        tau_ss = SampleSummarizer(self.samples['unique_taus'])
        self.assertAlmostEqual(tau_ss.mean, tau.mean, places=0)
        self.assertAlmostEqual(tau_ss.variance, tau.variance, places=0)
        self.assertAlmostEqual(d_theta_ss.mean, theta.mean, places=2)
        self.assertAlmostEqual(d_theta_ss.variance, theta.variance, places=2)
        self.assertAlmostEqual(a_theta_ss.mean, a_theta.mean, places=2)
        self.assertAlmostEqual(a_theta_ss.variance, a_theta.variance, places=2)
        psi_freqs = get_freqs(self.samples['psi'])
        for psi, freq in psi_freqs.iteritems():
            if psi == 2:
                self.assertAlmostEqual(freq, 2/float(5), places=1)
            else:
                self.assertAlmostEqual(freq, 1/float(5), places=1)

    def test_new_uniform_constrained(self):
        if not test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        self.preamble = """
concentrationShape = 0
concentrationScale = 0
thetaShape = 1.0
thetaScale = 0.01
ancestralThetaShape = 0
ancestralThetaScale = 0
thetaParameters = 012
tauShape = 2.0
tauScale = 5.0
bottleProportionShapeA = 0
bottleProportionShapeB = 0
bottleProportionShared = 1
numTauClasses = 4
migrationShape = 0
migrationScale = 0
recombinationShape = 0
recombinationScale = 0
constrain = 0
subParamConstrain = 111111111
"""

        self.generate_prior(sample_size=1000, batch_size=250, np=4)
        self.assertFalse(self.samples == None)
        theta = GammaDistribution(1.0, 0.01)
        tau = GammaDistribution(2.0, 5.0)
        not_equal_failures = 0
        for i in range(1000):
            for j in range(4):
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d1_thetas'][i][j]:
                    not_equal_failures += 1
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
                if self.samples['d1_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
        self.assertTrue(not_equal_failures < 10)
        d_theta_ss = SampleSummarizer(self.samples['all_d1_thetas'] + \
                self.samples['all_d2_thetas'] + \
                self.samples['all_a_thetas'])
        tau_ss = SampleSummarizer(self.samples['unique_taus'])
        self.assertAlmostEqual(tau_ss.mean, tau.mean, places=0)
        self.assertAlmostEqual(tau_ss.variance, tau.variance, places=0)
        self.assertAlmostEqual(d_theta_ss.mean, theta.mean, places=2)
        self.assertAlmostEqual(d_theta_ss.variance, theta.variance, places=2)
        self.assertEqual(self.samples['psi'], [4]*1000)
        for t in self.samples['taus']:
            self.assertEqual(len(set(t)), 4)

    def test_new_uniform_with_uniform_tau(self):
        if not test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        self.preamble = """
concentrationShape = 0
concentrationScale = 0
thetaShape = 1.0
thetaScale = 0.01
ancestralThetaShape = 0
ancestralThetaScale = 0
thetaParameters = 012
tauShape = -1.0
tauScale = -4.0
bottleProportionShapeA = 10
bottleProportionShapeB = 1
bottleProportionShared = 0
numTauClasses = 0
migrationShape = 0
migrationScale = 0
recombinationShape = 0
recombinationScale = 0
constrain = 0
subParamConstrain = 111111111
"""

        self.generate_prior(sample_size=1000, batch_size=250, np=4)
        self.assertFalse(self.samples == None)
        theta = GammaDistribution(1.0, 0.01)
        tau = ContinuousUniformDistribution(1.0, 4.0)
        not_equal_failures = 0
        for i in range(1000):
            for j in range(4):
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d1_thetas'][i][j]:
                    not_equal_failures += 1
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
                if self.samples['d1_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
        self.assertTrue(not_equal_failures < 10)
        d_theta_ss = SampleSummarizer(self.samples['all_d1_thetas'] + \
                self.samples['all_d2_thetas'] + \
                self.samples['all_a_thetas'])
        tau_ss = SampleSummarizer(self.samples['unique_taus'])
        self.assertAlmostEqual(tau_ss.mean, tau.mean, places=0)
        self.assertAlmostEqual(tau_ss.variance, tau.variance, places=0)
        self.assertAlmostEqual(d_theta_ss.mean, theta.mean, places=2)
        self.assertAlmostEqual(d_theta_ss.variance, theta.variance, places=2)
        psi_freqs = get_freqs(self.samples['psi'])
        for psi, freq in psi_freqs.iteritems():
            if psi == 2:
                self.assertAlmostEqual(freq, 2/float(5), places=1)
            else:
                self.assertAlmostEqual(freq, 1/float(5), places=1)

    def test_new_uniform_constrained_with_uniform_tau(self):
        if not test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        self.preamble = """
concentrationShape = 0
concentrationScale = 0
thetaShape = 1.0
thetaScale = 0.01
ancestralThetaShape = 0
ancestralThetaScale = 0
thetaParameters = 012
tauShape = -1.0
tauScale = -4.0
bottleProportionShapeA = 10
bottleProportionShapeB = 1
bottleProportionShared = 0
numTauClasses = 4
migrationShape = 0
migrationScale = 0
recombinationShape = 0
recombinationScale = 0
constrain = 0
subParamConstrain = 111111111
"""

        self.generate_prior(sample_size=1000, batch_size=250, np=4)
        self.assertFalse(self.samples == None)
        theta = GammaDistribution(1.0, 0.01)
        tau = ContinuousUniformDistribution(1.0, 4.0)
        not_equal_failures = 0
        for i in range(1000):
            for j in range(4):
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d1_thetas'][i][j]:
                    not_equal_failures += 1
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
                if self.samples['d1_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
        self.assertTrue(not_equal_failures < 10)
        d_theta_ss = SampleSummarizer(self.samples['all_d1_thetas'] + \
                self.samples['all_d2_thetas'] + \
                self.samples['all_a_thetas'])
        tau_ss = SampleSummarizer(self.samples['unique_taus'])
        self.assertAlmostEqual(tau_ss.mean, tau.mean, places=0)
        self.assertAlmostEqual(tau_ss.variance, tau.variance, places=0)
        self.assertAlmostEqual(d_theta_ss.mean, theta.mean, places=2)
        self.assertAlmostEqual(d_theta_ss.variance, theta.variance, places=2)
        self.assertEqual(self.samples['psi'], [4]*1000)
        for t in self.samples['taus']:
            self.assertEqual(len(set(t)), 4)

    def test_new_uniform_psi(self):
        if not test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        self.preamble = """
concentrationShape = -1
concentrationScale = -1
thetaShape = 0.5
thetaScale = 0.01
ancestralThetaShape = 0
ancestralThetaScale = 0
thetaParameters = 012
tauShape = 1.0
tauScale = 3.0
bottleProportionShapeA = 10
bottleProportionShapeB = 1
bottleProportionShared = 0
numTauClasses = 0
migrationShape = 0
migrationScale = 0
recombinationShape = 0
recombinationScale = 0
constrain = 0
subParamConstrain = 111111111
"""

        self.generate_prior(sample_size=1000, batch_size=250, np=4)
        self.assertFalse(self.samples == None)
        theta = GammaDistribution(0.5, 0.01)
        tau = GammaDistribution(1.0, 3.0)
        not_equal_failures = 0
        for i in range(1000):
            for j in range(4):
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d1_thetas'][i][j]:
                    not_equal_failures += 1
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
                if self.samples['d1_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
        self.assertTrue(not_equal_failures < 10)
        d_theta_ss = SampleSummarizer(self.samples['all_d1_thetas'] + \
                self.samples['all_d2_thetas'] + \
                self.samples['all_a_thetas'])
        tau_ss = SampleSummarizer(self.samples['unique_taus'])
        self.assertAlmostEqual(tau_ss.mean, tau.mean, places=0)
        self.assertAlmostEqual(tau_ss.variance, tau.variance, places=0)
        self.assertAlmostEqual(d_theta_ss.mean, theta.mean, places=2)
        self.assertAlmostEqual(d_theta_ss.variance, theta.variance, places=2)
        psi_freqs = get_freqs(self.samples['psi'])
        for psi, freq in psi_freqs.iteritems():
            self.assertAlmostEqual(freq, 0.25, places=1)

    def test_new_dpp_clustered(self):
        if not test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        self.preamble = """
concentrationShape = 1000
concentrationScale = 0.0000000001
thetaShape = 1.0
thetaScale = 0.01
ancestralThetaShape = 0
ancestralThetaScale = 0
thetaParameters = 012
tauShape = 1.0
tauScale = 5.0
bottleProportionShapeA = 0
bottleProportionShapeB = 0
bottleProportionShared = 1
numTauClasses = 0
migrationShape = 0
migrationScale = 0
recombinationShape = 0
recombinationScale = 0
constrain = 0
subParamConstrain = 111111111
"""

        self.generate_prior(sample_size=1000, batch_size=250, np=4)
        self.assertFalse(self.samples == None)
        theta = GammaDistribution(1.0, 0.01)
        tau = GammaDistribution(1.0, 5.0)
        not_equal_failures = 0
        for i in range(1000):
            for j in range(4):
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d1_thetas'][i][j]:
                    not_equal_failures += 1
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
                if self.samples['d1_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
        self.assertTrue(not_equal_failures < 10)
        d_theta_ss = SampleSummarizer(self.samples['all_d1_thetas'] + \
                self.samples['all_d2_thetas'] + \
                self.samples['all_a_thetas'])
        tau_ss = SampleSummarizer(self.samples['unique_taus'])
        self.assertAlmostEqual(tau_ss.mean, tau.mean, places=0)
        self.assertTrue(math.fabs(tau_ss.variance - tau.variance) < 3.0,
                msg='unique tau variance is {0}; expecting {1}'.format(
                        tau_ss.variance, tau.variance))
        self.assertAlmostEqual(d_theta_ss.mean, theta.mean, places=2)
        self.assertAlmostEqual(d_theta_ss.variance, theta.variance, places=2)
        self.assertEqual(self.samples['psi'], [1]*1000)
        for t in self.samples['taus']:
            self.assertEqual(len(set(t)), 1)

    def test_new_dpp_dispersed(self):
        if not test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        self.preamble = """
concentrationShape = 1000
concentrationScale = 1000
thetaShape = 1.0
thetaScale = 0.01
ancestralThetaShape = 0
ancestralThetaScale = 0
thetaParameters = 012
tauShape = 1.0
tauScale = 5.0
bottleProportionShapeA = 0
bottleProportionShapeB = 0
bottleProportionShared = 1
numTauClasses = 0
migrationShape = 0
migrationScale = 0
recombinationShape = 0
recombinationScale = 0
constrain = 0
subParamConstrain = 111111111
"""

        self.generate_prior(sample_size=1000, batch_size=250, np=4)
        self.assertFalse(self.samples == None)
        theta = GammaDistribution(1.0, 0.01)
        tau = GammaDistribution(1.0, 5.0)
        not_equal_failures = 0
        for i in range(1000):
            for j in range(4):
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d1_thetas'][i][j]:
                    not_equal_failures += 1
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
                if self.samples['d1_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
        self.assertTrue(not_equal_failures < 10)
        d_theta_ss = SampleSummarizer(self.samples['all_d1_thetas'] + \
                self.samples['all_d2_thetas'] + \
                self.samples['all_a_thetas'])
        tau_ss = SampleSummarizer(self.samples['unique_taus'])
        self.assertAlmostEqual(tau_ss.mean, tau.mean, places=0)
        self.assertTrue(math.fabs(tau_ss.variance - tau.variance) < 3.0,
                msg='unique tau variance is {0}; expecting {1}'.format(
                        tau_ss.variance, tau.variance))
        self.assertAlmostEqual(d_theta_ss.mean, theta.mean, places=2)
        self.assertAlmostEqual(d_theta_ss.variance, theta.variance, places=2)
        self.assertEqual(self.samples['psi'], [4]*1000)
        for t in self.samples['taus']:
            self.assertEqual(len(set(t)), 4)

    def test_new_dpp_2(self):
        if not test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        self.preamble = """
concentrationShape = 1000
concentrationScale = 0.00088
thetaShape = 1.0
thetaScale = 0.01
ancestralThetaShape = 0
ancestralThetaScale = 0
thetaParameters = 012
tauShape = 1.0
tauScale = 5.0
bottleProportionShapeA = 0
bottleProportionShapeB = 0
bottleProportionShared = 1
numTauClasses = 0
migrationShape = 0
migrationScale = 0
recombinationShape = 0
recombinationScale = 0
constrain = 0
subParamConstrain = 111111111
"""

        self.generate_prior(sample_size=1000, batch_size=250, np=4)
        self.assertFalse(self.samples == None)
        theta = GammaDistribution(1.0, 0.01)
        tau = GammaDistribution(1.0, 5.0)
        not_equal_failures = 0
        for i in range(1000):
            for j in range(4):
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d1_thetas'][i][j]:
                    not_equal_failures += 1
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
                if self.samples['d1_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
        self.assertTrue(not_equal_failures < 10)
        d_theta_ss = SampleSummarizer(self.samples['all_d1_thetas'] + \
                self.samples['all_d2_thetas'] + \
                self.samples['all_a_thetas'])
        tau_ss = SampleSummarizer(self.samples['unique_taus'])
        self.assertAlmostEqual(tau_ss.mean, tau.mean, places=0)
        self.assertTrue(math.fabs(tau_ss.variance - tau.variance) < 3.0,
                msg='unique tau variance is {0}; expecting {1}'.format(
                        tau_ss.variance, tau.variance))
        self.assertAlmostEqual(d_theta_ss.mean, theta.mean, places=2)
        self.assertAlmostEqual(d_theta_ss.variance, theta.variance, places=2)
        psi_freqs = get_freqs(self.samples['psi'])
        self.assertAlmostEqual(psi_freqs[2], 0.46078, places=1)
        self.assertEqual(mode_list(self.samples['psi']), [2])

    def test_new_dpp_2_with_uniform_tau(self):
        if not test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        self.preamble = """
concentrationShape = 1000
concentrationScale = 0.00088
thetaShape = 1.0
thetaScale = 0.01
ancestralThetaShape = 0
ancestralThetaScale = 0
thetaParameters = 012
tauShape = 2.0
tauScale = -10.0
bottleProportionShapeA = 0
bottleProportionShapeB = 0
bottleProportionShared = 1
numTauClasses = 0
migrationShape = 0
migrationScale = 0
recombinationShape = 0
recombinationScale = 0
constrain = 0
subParamConstrain = 111111111
"""

        self.generate_prior(sample_size=1000, batch_size=250, np=4)
        self.assertFalse(self.samples == None)
        theta = GammaDistribution(1.0, 0.01)
        tau = ContinuousUniformDistribution(2.0, 10.0)
        not_equal_failures = 0
        for i in range(1000):
            for j in range(4):
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d1_thetas'][i][j]:
                    not_equal_failures += 1
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
                if self.samples['d1_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
        self.assertTrue(not_equal_failures < 10)
        d_theta_ss = SampleSummarizer(self.samples['all_d1_thetas'] + \
                self.samples['all_d2_thetas'] + \
                self.samples['all_a_thetas'])
        tau_ss = SampleSummarizer(self.samples['unique_taus'])
        self.assertAlmostEqual(tau_ss.mean, tau.mean, places=0)
        self.assertTrue(math.fabs(tau_ss.variance - tau.variance) < 3.0,
                msg='unique tau variance is {0}; expecting {1}'.format(
                        tau_ss.variance, tau.variance))
        self.assertAlmostEqual(d_theta_ss.mean, theta.mean, places=2)
        self.assertAlmostEqual(d_theta_ss.variance, theta.variance, places=2)
        psi_freqs = get_freqs(self.samples['psi'])
        self.assertAlmostEqual(psi_freqs[2], 0.46078, places=1)
        self.assertEqual(mode_list(self.samples['psi']), [2])

    def test_new_dpp_constrained_with_uniform_tau(self):
        if not test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        self.preamble = """
concentrationShape = 1000
concentrationScale = 0.00088
thetaShape = 1.0
thetaScale = 0.01
ancestralThetaShape = 0
ancestralThetaScale = 0
thetaParameters = 012
tauShape = 2.0
tauScale = -10.0
bottleProportionShapeA = 0
bottleProportionShapeB = 0
bottleProportionShared = 1
numTauClasses = 4
migrationShape = 0
migrationScale = 0
recombinationShape = 0
recombinationScale = 0
constrain = 0
subParamConstrain = 111111111
"""

        self.generate_prior(sample_size=1000, batch_size=250, np=4)
        self.assertFalse(self.samples == None)
        theta = GammaDistribution(1.0, 0.01)
        tau = ContinuousUniformDistribution(2.0, 10.0)
        not_equal_failures = 0
        for i in range(1000):
            for j in range(4):
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d1_thetas'][i][j]:
                    not_equal_failures += 1
                if self.samples['a_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
                if self.samples['d1_thetas'][i][j] == \
                        self.samples['d2_thetas'][i][j]:
                    not_equal_failures += 1
        self.assertTrue(not_equal_failures < 10)
        d_theta_ss = SampleSummarizer(self.samples['all_d1_thetas'] + \
                self.samples['all_d2_thetas'] + \
                self.samples['all_a_thetas'])
        tau_ss = SampleSummarizer(self.samples['unique_taus'])
        self.assertAlmostEqual(tau_ss.mean, tau.mean, places=0)
        self.assertTrue(math.fabs(tau_ss.variance - tau.variance) < 3.0,
                msg='unique tau variance is {0}; expecting {1}'.format(
                        tau_ss.variance, tau.variance))
        self.assertAlmostEqual(d_theta_ss.mean, theta.mean, places=2)
        self.assertAlmostEqual(d_theta_ss.variance, theta.variance, places=2)
        self.assertEqual(self.samples['psi'], [4]*1000)
        for t in self.samples['taus']:
            self.assertEqual(len(set(t)), 4)

if __name__ == '__main__':
    unittest.main()

