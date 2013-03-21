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

    def _update_config(self, cfg, params, multi_locus=False):
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

        self.assertSameDistributions(c.psi, DiscreteUniformDistribution(1,4))
        self.assertSameDistributions(c.tau, ContinuousUniformDistribution(0,10))
        self.assertSameDistributions(c.recombination, ContinuousUniformDistribution(0,0))
        self.assertSameDistributions(c.migration, ContinuousUniformDistribution(0,0))
        self.assertSameDistributions(c.a_theta, ContinuousUniformDistribution(0.0001,0.1))
        self.assertSameDistributions(c.theta, ContinuousUniformDistribution(0.0001,0.1))
        self.assertSameDistributions(c.d_theta, BetaDistribution(1,1,2))

if __name__ == '__main__':
    unittest.main()

