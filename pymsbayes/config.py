#! /usr/bin/env python

import os
import sys
import re
import math
from cStringIO import StringIO
import ConfigParser

from pymsbayes import fileio
from pymsbayes.utils import probability
from pymsbayes.utils.messaging import get_logger

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
        self.bottle_proportion = None
        self.bottle_proportion_shared = None
        self.taxa = []
        self.sample_table = []
        self._parse_config(cfg_file)

    @classmethod
    def is_config(cls, cfg_file):
        cfg_stream, close = fileio.process_file_arg(cfg_file)
        for i in range(100):
            try:
                line = cfg_stream.next()
                if cls._begin_pattern.match(line.strip()):
                    if close:
                        cfg_stream.close()
                    return True
            except:
                if close:
                    cfg_stream.close()
                return False
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
        cfg_stream, close = fileio.process_file_arg(cfg)
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
        self.sample_table = []
        for i, row in enumerate(table):
            r = row.strip().split()
            if r and (r[0] not in self.taxa):
                self.taxa.append(r[0])
            self.sample_table.append(r)
        self.npairs = len(self.taxa)
    
    def _get_gamma_or_uniform_distribution(self, shape, scale):
        if shape > 0.0 and scale > 0.0:
            return probability.GammaDistribution(shape, scale)
        else:
            a = math.fabs(shape)
            b = math.fabs(scale)
            if a < b:
                return probability.ContinuousUniformDistribution(a, b)
            else:
                return probability.ContinuousUniformDistribution(b, a)

    def _get_gamma_or_uniform_settings(self, distribution):
        if isinstance(distribution, probability.GammaDistribution):
            return distribution.shape, distribution.scale
        return -distribution.minimum, -distribution.maximum

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
                self.psi = probability.DiscreteUniformDistribution(psi, psi)
                self.dpp_concentration = None
                self.div_model_prior = 'constrained'
            elif (dpp_concentration_shape > 0.0) and (
                    dpp_concentration_scale > 0.0):
                self.div_model_prior = 'dpp'
                self.dpp_concentration = probability.GammaDistribution(
                        dpp_concentration_shape,
                        dpp_concentration_scale)
            elif (dpp_concentration_shape > -1.0) and (
                    dpp_concentration_scale > -1.0):
                self.div_model_prior = 'uniform'
                self.dpp_concentration = None
            else:
                self.div_model_prior = 'psi'
                self.dpp_concentration = None
                self.psi = probability.DiscreteUniformDistribution(1,
                        self.npairs)

            tau_shape = float(kwargs.get('taushape', 1.0))
            tau_scale = float(kwargs.get('tauscale', 2.0))
            self.tau = self._get_gamma_or_uniform_distribution(tau_shape,
                    tau_scale)
            theta_shape = float(kwargs.get('thetashape', 1.0))
            theta_scale = float(kwargs.get('thetascale', 0.001))
            anc_theta_shape = float(kwargs.get('ancestralthetashape', 0.0))
            anc_theta_scale = float(kwargs.get('ancestralthetascale', 0.0))
            theta_params = kwargs.get('thetaparameters', '012')
            self.theta_parameters = [int(i) for i in list(theta_params)]
            self.theta = None
            self.d_theta = self._get_gamma_or_uniform_distribution(theta_shape,
                    theta_scale)
            if ((anc_theta_shape > 0.0) and (anc_theta_scale > 0.0)) and (
                    (theta_params == '001') or (theta_params == '012')):
                self.a_theta = self._get_gamma_or_uniform_distribution(
                        anc_theta_shape,
                        anc_theta_scale)
            else:
                self.a_theta = self._get_gamma_or_uniform_distribution(
                        theta_shape,
                        theta_scale)
            mig_shape = float(kwargs.get('migrationshape', 0.0))
            mig_scale = float(kwargs.get('migrationscale', 0.0))
            if (mig_shape > 0.0) and (mig_scale > 0.0):
                self.migration = probability.GammaDistribution(mig_shape,
                        mig_scale)
            else:
                self.migration = probability.ContinuousUniformDistribution(0.0,
                        0.0)
            rec_shape = float(kwargs.get('recombinationshape', 0.0))
            rec_scale = float(kwargs.get('recombinationscale', 0.0))
            if (rec_shape > 0.0) and (rec_scale > 0.0):
                self.recombination = probability.GammaDistribution(rec_shape,
                        rec_scale)
            else:
                self.recombination = probability.ContinuousUniformDistribution(
                        0.0, 0.0)
            bottle_a = float(kwargs.get('bottleProportionShapeA', 0.0))
            bottle_b = float(kwargs.get('bottleProportionShapeB', 0.0))
            if (bottle_a > 0.0) and (bottle_b > 0.0):
                self.bottle_proportion = probability.BetaDistribution(bottle_a,
                        bottle_b)
            else:
                self.bottle_proportion = None
            bottle_shared = int(kwargs.get('bottleProportionShared', 0))
            self.bottle_proportion_shared = bool(bottle_shared)
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
            self.psi = probability.DiscreteUniformDistribution(1, self.npairs)
            self.div_model_prior = 'psi'
        else:
            self.psi = probability.DiscreteUniformDistribution(psi, psi)
            self.div_model_prior = 'constrained'
        self.tau = probability.ContinuousUniformDistribution(ltau, utau)
        self.theta = probability.ContinuousUniformDistribution(ltheta, utheta)
        self.a_theta = probability.ContinuousUniformDistribution(ltheta,
                utheta*anc_scalar)
        self.d_theta = probability.BetaDistribution(1, 1, 2)
        self.migration = probability.ContinuousUniformDistribution(lmig, umig)
        self.recombination = probability.ContinuousUniformDistribution(lrec,
                urec)

    def _parse_preamble(self, preamble):
        cfg = ConfigParser.RawConfigParser()
        cfg.readfp(preamble)
        d = {}
        for tup in cfg.items('preamble'):
            d[tup[0].lower()] = tup[1]
        return d

    def _get_old_settings(self):
        d = {}
        d['upperTheta'] = self.theta.maximum
        d['lowerTheta'] = self.theta.minimum
        d['upperTau'] = self.tau.maximum
        d['lowerTau'] = self.tau.minimum
        d['upperAncPopSize'] = self.a_theta.maximum /self.theta.maximum
        if self.div_model_prior.lower() == 'constrained':
            assert self.psi.maximum == self.psi.minimum
            d['numTauClasses'] = self.psi.maximum
        else:
            d['numTauClasses'] = 0
        d['upperMig'] = self.migration.maximum
        d['upperRec'] = self.recombination.maximum
        d['constrain'] = 0
        d['subParamConstrain'] = '111111111'
        return d
    
    def _get_new_settings(self):
        d = {}
        if self.dpp_concentration:
            d['concentrationShape'] = self.dpp_concentration.shape
            d['concentrationScale'] = self.dpp_concentration.scale
        elif self.div_model_prior.lower() == 'psi':
            d['concentrationShape'] = -1.0
            d['concentrationScale'] = -1.0
        else:
            d['concentrationShape'] = 0.0
            d['concentrationScale'] = 0.0
        d['thetaShape'], d['thetaScale'] = \
                self._get_gamma_or_uniform_settings(self.d_theta)
        d['ancestralThetaShape'], d['ancestralThetaScale'] = \
                self._get_gamma_or_uniform_settings(self.a_theta)
        d['thetaParameters'] = ''.join([str(x) for x in self.theta_parameters])
        d['tauShape'], d['tauScale'] = \
                self._get_gamma_or_uniform_settings(self.tau)
        if self.bottle_proportion:
            d['bottleProportionShapeA'] = self.bottle_proportion.alpha
            d['bottleProportionShapeB'] = self.bottle_proportion.beta
        else:
            d['bottleProportionShapeA'] = 0.0
            d['bottleProportionShapeB'] = 0.0
        if self.bottle_proportion_shared:
            d['bottleProportionShared'] = 1
        else:
            d['bottleProportionShared'] = 0
        if self.div_model_prior.lower() == 'constrained':
            assert self.psi.maximum == self.psi.minimum
            d['numTauClasses'] = self.psi.maximum
        else:
            d['numTauClasses'] = 0
        d['migrationShape'], d['migrationScale'] = \
                self._get_gamma_or_uniform_settings(self.migration)
        d['recombinationShape'], d['recombinationScale'] = \
                self._get_gamma_or_uniform_settings(self.recombination)
        d['constrain'] = 0
        d['subParamConstrain'] = '111111111'
        return d

    def get_settings(self):
        if self.implementation.lower() == 'old':
            return self._get_old_settings()
        return self._get_new_settings()

    def get_preamble(self):
        s = StringIO()
        settings = self.get_settings()
        for k in sorted(settings.iterkeys()):
            s.write('{0} = {1}\n'.format(k, settings[k]))
        return s.getvalue()

    def get_sample_table(self):
        s = StringIO()
        s.write('BEGIN SAMPLE_TBL\n')
        for row in self.sample_table:
            s.write('{0}\n'.format('\t'.join(row)))
        s.write('END SAMPLE_TBL\n')
        return s.getvalue()
    
    def __str__(self):
        return '\n'.join([self.get_preamble() + self.get_sample_table()])

    def write(self, file_obj):
        out, close = fileio.process_file_arg(file_obj, 'w')
        out.write(self.__str__())
        if close:
            out.close()

