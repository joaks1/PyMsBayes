#! /usr/bin/env python

import unittest
import os
from cStringIO import StringIO

from pymsbayes.config import MsBayesConfig
from pymsbayes.utils.probability import *
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class MsBayesConfigTestCase(PyMsBayesTestCase):
    
    def setUp(self):
        self.cfg = StringIO()

    def _update_config(self, cfg, params, multi_locus=False, new_impl=False):
        if new_impl:
            cfg.write("""concentrationShape = {c_shape}
concentrationScale = {c_scale}
thetaShape = {theta_shape}
thetaScale = {theta_scale}
ancestralThetaShape = {a_theta_shape}
ancestralThetaScale = {a_theta_scale}
thetaParameters = {theta_params}
tauShape = {tau_shape}
tauScale = {tau_scale}
migrationShape = {migration_shape}
migrationScale = {migration_scale}
recombinationShape = {recombination_shape}
recombinationScale = {recombination_scale}
numTauClasses = {psi}
constrain = 0
subParamConstrain = 111111111
""".format(**params))
        else:
            cfg.write("""upperTheta = 0.1
lowerTheta = {ltheta}
upperTau = {utau}
numTauClasses = {psi}
upperMig = {umig}
upperRec = {urec}
upperAncPopSize = {atheta}
constrain = 0
subParamConstrain = 111111111
""".format(**params))
        if not multi_locus:
            cfg.write("""
BEGIN SAMPLE_TBL
pair0	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair0.locus0.fasta
pair1	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair2.locus0.fasta
pair2	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair3.locus0.fasta
pair3	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair3.locus0.fasta
END SAMPLE_TBL
""")
        else:
            cfg.write("""
BEGIN SAMPLE_TBL
pair0	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair0.locus0.fasta
pair0	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair0.locus0.fasta
pair1	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair2.locus0.fasta
pair1	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair2.locus0.fasta
pair1	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair2.locus0.fasta
pair1	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair2.locus0.fasta
pair2	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair3.locus0.fasta
pair3	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair3.locus0.fasta
pair3	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair3.locus0.fasta
END SAMPLE_TBL
""")
        cfg.seek(0)

    def test_single_locus(self):
        p = {'ltheta': 0.0001,
             'utheta': 0.1,
             'utau': 10.0,
             'psi': 0,
             'umig': 0.0,
             'urec': 0.0,
             'atheta': 1.0,}
        self._update_config(self.cfg, p)
        _LOG.debug('testing config:\n\n{0}\n'.format(self.cfg.getvalue()))
        c = MsBayesConfig(self.cfg)
        self.assertEqual(c.npairs, 4)
        self.assertIsInstance(c.psi, DiscreteUniformDistribution)
        self.assertIsInstance(c.tau, ContinuousUniformDistribution)
        self.assertIsInstance(c.recombination, ContinuousUniformDistribution)
        self.assertIsInstance(c.migration, ContinuousUniformDistribution)
        self.assertIsInstance(c.a_theta, ContinuousUniformDistribution)
        self.assertIsInstance(c.theta, ContinuousUniformDistribution)
        self.assertIsInstance(c.d_theta, BetaDistribution)
        self.assertEqual(c.div_model_prior, 'psi')
        self.assertEqual(c.dpp_concentration, None)
        self.assertEqual(c.theta_parameters, None)
        self.assertEqual(c.implementation, 'old')

        self.assertSameDistributions(c.psi, DiscreteUniformDistribution(1,4))
        self.assertSameDistributions(c.tau, ContinuousUniformDistribution(0,10))
        self.assertSameDistributions(c.recombination, ContinuousUniformDistribution(0,0))
        self.assertSameDistributions(c.migration, ContinuousUniformDistribution(0,0))
        self.assertSameDistributions(c.a_theta, ContinuousUniformDistribution(0.0001,0.1))
        self.assertSameDistributions(c.theta, ContinuousUniformDistribution(0.0001,0.1))
        self.assertSameDistributions(c.d_theta, BetaDistribution(1,1,2))

    def test_multi_locus(self):
        p = {'ltheta': 0.0001,
             'utheta': 0.1,
             'utau': 10.0,
             'psi': 0,
             'umig': 0.0,
             'urec': 0.0,
             'atheta': 1.0,}
        self._update_config(self.cfg, p, multi_locus=True)
        _LOG.debug('testing config:\n\n{0}\n'.format(self.cfg.getvalue()))
        c = MsBayesConfig(self.cfg)
        self.assertEqual(c.npairs, 4)
        self.assertIsInstance(c.psi, DiscreteUniformDistribution)
        self.assertIsInstance(c.tau, ContinuousUniformDistribution)
        self.assertIsInstance(c.recombination, ContinuousUniformDistribution)
        self.assertIsInstance(c.migration, ContinuousUniformDistribution)
        self.assertIsInstance(c.a_theta, ContinuousUniformDistribution)
        self.assertIsInstance(c.theta, ContinuousUniformDistribution)
        self.assertIsInstance(c.d_theta, BetaDistribution)
        self.assertEqual(c.div_model_prior, 'psi')
        self.assertEqual(c.dpp_concentration, None)
        self.assertEqual(c.theta_parameters, None)
        self.assertEqual(c.implementation, 'old')

        self.assertSameDistributions(c.psi, DiscreteUniformDistribution(1,4))
        self.assertSameDistributions(c.tau, ContinuousUniformDistribution(0,10))
        self.assertSameDistributions(c.recombination, ContinuousUniformDistribution(0,0))
        self.assertSameDistributions(c.migration, ContinuousUniformDistribution(0,0))
        self.assertSameDistributions(c.a_theta, ContinuousUniformDistribution(0.0001,0.1))
        self.assertSameDistributions(c.theta, ContinuousUniformDistribution(0.0001,0.1))
        self.assertSameDistributions(c.d_theta, BetaDistribution(1,1,2))

    def test_new_implementation(self):
        p = {'c_shape': 2,
             'c_scale': 5,
             'theta_shape': 2.0,
             'theta_scale': 0.002,
             'tau_shape': 1.0,
             'tau_scale': 5.0,
             'a_theta_shape': 1,
             'a_theta_scale': 0.01,
             'theta_params': '001',
             'migration_shape': 1.5,
             'migration_scale': 0.1,
             'recombination_shape': 0.0,
             'recombination_scale': 0.0,
             'psi': 0,}
        self._update_config(self.cfg, p, new_impl=True)
        _LOG.debug('testing config:\n\n{0}\n'.format(self.cfg.getvalue()))
        c = MsBayesConfig(self.cfg)
        self.assertEqual(c.npairs, 4)
        self.assertEqual(c.psi, None)
        self.assertIsInstance(c.tau, GammaDistribution)
        self.assertIsInstance(c.recombination, ContinuousUniformDistribution)
        self.assertIsInstance(c.migration, GammaDistribution)
        self.assertIsInstance(c.a_theta, GammaDistribution)
        self.assertEqual(c.theta, None)
        self.assertIsInstance(c.d_theta, GammaDistribution)
        self.assertEqual(c.div_model_prior, 'dpp')
        self.assertIsInstance(c.dpp_concentration, GammaDistribution)
        self.assertEqual(c.theta_parameters, [0, 0, 1])
        self.assertEqual(c.implementation, 'new')

        self.assertSameDistributions(c.dpp_concentration, GammaDistribution(2, 5))
        self.assertSameDistributions(c.tau, GammaDistribution(1, 5))
        self.assertSameDistributions(c.recombination, ContinuousUniformDistribution(0,0))
        self.assertSameDistributions(c.migration, GammaDistribution(1.5, 0.1))
        self.assertSameDistributions(c.a_theta, GammaDistribution(1, 0.01))
        self.assertSameDistributions(c.d_theta, GammaDistribution(2, 0.002))

if __name__ == '__main__':
    unittest.main()

