#! /usr/bin/env python

import os
import sys
import re
import math
from cStringIO import StringIO
import ConfigParser

from pymsbayes.utils.messaging import get_logger
from pymsbayes.fileio import process_file_arg
from pymsbayes.utils.probability import (ContinuousUniformDistribution,
        BetaDistribution, DiscreteUniformDistribution, GammaDistribution)

_LOG = get_logger(__name__)

class MsBayesConfig(object):
    _begin_pattern = re.compile(r'^begin\s*sample_tbl$', re.IGNORECASE)
    _end_pattern = re.compile(r'^end\s*sample_tbl$', re.IGNORECASE)

    def __init__(self, cfg_file):
        self.theta = None
        self.d_theta = None
        self.a_theta = None
        self.tau = None
        self.migration = None
        self.recombination = None
        self.psi = None
        self.dpp_concentration = None
        self.npairs = None
        self.implementation = 'old'
        self.div_model_prior = None
        self.theta_parameters = None
        self._parse_config(cfg_file)

    @classmethod
    def is_config(cls, cfg_file):
        cfg_stream, close = process_file_arg(cfg_file)
        for i, line in enumerate(cfg_stream):
            if cls._begin_pattern.match(line.strip()):
                if close:
                    cfg_stream.close()
                return True
        if close:
            cfg_stream.close()
        return False

    def _parse_config(self, cfg):
        preamble, table = self._split_config(cfg)
        preamble.seek(0)
        table.seek(0)
        self._parse_table(table)
        self._set_priors(preamble)

    def _split_config(self, cfg):
        cfg_stream, close = process_file_arg(cfg)
        preamble = StringIO()
        table = StringIO()
        preamble.write('[preamble]\n')
        table_started = False
        for i, line in enumerate(cfg_stream):
            if self._end_pattern.match(line.strip()):
                break
            if self._begin_pattern.match(line.strip()):
                table_started = True
                continue
            if table_started:
                table.write(line)
            else:
                preamble.write(line)
        if close:
            cfg_stream.close()
        return preamble, table

    def _parse_table(self, table):
        taxa = set()
        for i, row in enumerate(table):
            r = row.strip()
            if r:
                taxa.add(r.split()[0])
        self.npairs = len(taxa)

    def _set_priors(self, preamble):
        kwargs = self._parse_preamble(preamble)
        psi = int(kwargs.get('numtauclasses', 0))
        for k in ['concentrationshape', 'concentrationscale', 'taushape',
                'tauscale', 'thetashape', 'thetascale', 'ancestralthetashape',
                'ancestralthetascale', 'thetaparameters']:
            if kwargs.has_key(k):
                self.implementation = 'new'
        
        if self.implementation == 'new':
            dpp_concentration_shape = float(kwargs.get('concentrationshape',
                    0.0))
            dpp_concentration_scale = float(kwargs.get('concentrationscale',
                    0.0))
            if psi != 0:
                self.psi = DiscreteUniformDistribution(psi, psi)
                self.dpp_concentration = None
                self.div_model_prior = 'constrained'
            elif (dpp_concentration_shape > 0.0) and (
                    dpp_concentration_scale > 0.0):
                self.div_model_prior = 'dpp'
                self.dpp_concentration = GammaDistribution(
                        dpp_concentration_shape,
                        dpp_concentration_scale)
            elif (dpp_concentration_shape > -1.0) and (
                    dpp_concentration_scale > -1.0):
                self.div_model_prior = 'uniform'
                self.dpp_concentration = None
            else:
                self.div_model_prior = 'psi'
                self.dpp_concentration = None

            tau_shape = float(kwargs.get('taushape', 1.0))
            tau_scale = float(kwargs.get('tauscale', 2.0))
            if tau_shape > 0.0 and tau_scale > 0.0:
                self.tau = GammaDistribution(tau_shape, tau_scale)
            else:
                a = math.fabs(tau_shape)
                b = math.fabs(tau_scale)
                if a < b:
                    self.tau = ContinuousUniformDistribution(a, b)
                else:
                    self.tau = ContinuousUniformDistribution(b, a)
            theta_shape = float(kwargs.get('thetashape', 1.0))
            theta_scale = float(kwargs.get('thetascale', 0.001))
            anc_theta_shape = float(kwargs.get('ancestralthetashape', 0.0))
            anc_theta_scale = float(kwargs.get('ancestralthetascale', 0.0))
            theta_params = kwargs.get('thetaparameters', '012')
            self.theta_parameters = [int(i) for i in list(theta_params)]
            self.theta = None
            self.d_theta = GammaDistribution(theta_shape, theta_scale)
            if ((anc_theta_shape > 0.0) and (anc_theta_scale > 0.0)) and (
                    (theta_params == '001') or (theta_params == '012')):
                self.a_theta = GammaDistribution(anc_theta_shape,
                        anc_theta_scale)
            else:
                self.a_theta = GammaDistribution(theta_shape, theta_scale)
            mig_shape = float(kwargs.get('migrationshape', 0.0))
            mig_scale = float(kwargs.get('migrationscale', 0.0))
            if (mig_shape > 0.0) and (mig_scale > 0.0):
                self.migration = GammaDistribution(mig_shape, mig_scale)
            else:
                self.migration = ContinuousUniformDistribution(0.0, 0.0)
            rec_shape = float(kwargs.get('recombinationshape', 0.0))
            rec_scale = float(kwargs.get('recombinationscale', 0.0))
            if (rec_shape > 0.0) and (rec_scale > 0.0):
                self.recombination = GammaDistribution(rec_shape, rec_scale)
            else:
                self.recombination = ContinuousUniformDistribution(0.0, 0.0)
            return

        ltau = float(kwargs.get('lowertau', 0.0))
        utau = float(kwargs.get('uppertau', 0.0))
        psi = int(kwargs['numtauclasses'])
        ltheta = float(kwargs.get('lowertheta', 0.0))
        utheta = float(kwargs['uppertheta'])
        anc_scalar = float(kwargs['upperancpopsize'])
        lmig = float(kwargs.get('lowermig', 0.0))
        umig = float(kwargs.get('uppermig', 0.0))
        lrec = float(kwargs.get('lowerrec', 0.0))
        urec = float(kwargs.get('upperrec', 0.0))
        if psi == 0:
            self.psi = DiscreteUniformDistribution(1, self.npairs)
            self.div_model_prior = 'psi'
        else:
            self.psi = DiscreteUniformDistribution(psi, psi)
            self.div_model_prior = 'constrained'
        self.tau = ContinuousUniformDistribution(ltau, utau)
        self.theta = ContinuousUniformDistribution(ltheta, utheta)
        self.a_theta = ContinuousUniformDistribution(ltheta, utheta*anc_scalar)
        self.d_theta = BetaDistribution(1, 1, 2)
        self.migration = ContinuousUniformDistribution(lmig, umig)
        self.recombination = ContinuousUniformDistribution(lrec, urec)

    def _parse_preamble(self, preamble):
        cfg = ConfigParser.RawConfigParser()
        cfg.readfp(preamble)
        d = {}
        for tup in cfg.items('preamble'):
            d[tup[0].lower()] = tup[1]
        return d

