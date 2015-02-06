#! /usr/bin/env python

import os
import sys
import re
import math
from cStringIO import StringIO
import ConfigParser

from pymsbayes import fileio
from pymsbayes.utils import probability, errors
from pymsbayes.utils.containers import OrderedDict
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class MsBayesConfig(object):
    _table_begin_pattern = re.compile(r'^begin\s*sample_tbl$', re.IGNORECASE)
    _table_end_pattern = re.compile(r'^end\s*sample_tbl$', re.IGNORECASE)

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
        self.sample_table = None
        self._parse_config(cfg_file)

    @classmethod
    def is_config(cls, cfg_file):
        cfg_stream, close = fileio.process_file_arg(cfg_file)
        for i in range(100):
            try:
                line = cfg_stream.next()
                if cls._table_begin_pattern.match(line.strip()):
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
            if self._table_end_pattern.match(line.strip()):
                table.write(line)
                break
            if self._table_begin_pattern.match(line.strip()):
                table_started = True
                table.write(line)
                continue
            if table_started:
                table.write(line)
            else:
                preamble.write(line)
        if close:
            cfg_stream.close()
        return preamble, table

    def _parse_table(self, table):
        self.sample_table = SampleTable(table)
        self.taxa = list(self.sample_table.taxa)
        self.npairs = len(self.taxa) 

    def _get_taxa(self):
        return self.sample_table.taxa

    taxa = property(_get_taxa)

    def _get_npairs(self):
        return self.sample_table.npairs

    npairs = property(_get_npairs)
    
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

    def __str__(self):
        return '\n'.join([self.get_preamble() + str(self.sample_table)])

    def write(self, file_obj):
        out, close = fileio.process_file_arg(file_obj, 'w')
        out.write(self.__str__())
        if close:
            out.close()


class SampleTable(object):
    _begin_pattern = re.compile(r'^begin\s*sample_tbl$', re.IGNORECASE)
    _end_pattern = re.compile(r'^end\s*sample_tbl$', re.IGNORECASE)

    def __init__(self, config_file):
        self.alignments = None
        self._ordering = []
        self._parse_table(config_file)

    def _parse_table(self, config_file):
        self.alignments = OrderedDict()
        cfg_stream, close = fileio.process_file_arg(config_file)
        try:
            table_started = False
            row_num = 0
            for i, l in enumerate(cfg_stream):
                line = l.strip()
                if self._end_pattern.match(line):
                    if not table_started:
                        raise errors.SampleTableError(
                                'hit end of sample table before beginning')
                    if len(self.alignments) < 1:
                        raise errors.SampleTableError(
                                'no rows found in sample table')
                    break
                if self._begin_pattern.match(line):
                    table_started = True
                    continue
                if not table_started:
                    continue
                if (line == '') or (line.startswith('#')):
                    continue
                row_num += 1
                try:
                    al = AlignmentConfig(line)
                except errors.SampleTableRowError as e:
                    _LOG.error('sample table row {0} is invalid'.format(
                            row_num))
                    raise e
                if not al.taxon_name in self.alignments:
                    self.alignments[al.taxon_name] = OrderedDict()
                    self.alignments[al.taxon_name][al.locus_name] = al
                    self._ordering.append((al.taxon_name, al.locus_name))
                    continue
                if al.locus_name in self.alignments[al.taxon_name]:
                    raise errors.SampleTableError('locus {0} found twice '
                            'for taxon {1} at row {2} of sample '
                            'table'.format(al.locus_name, al.taxon_name,
                                    row_num))
                self.alignments[al.taxon_name][al.locus_name] = al
                self._ordering.append((al.taxon_name, al.locus_name))
        finally:
            if close:
                cfg_stream.close()

    def _get_taxa(self):
        return self.alignments.keys()

    taxa = property(_get_taxa)

    def _get_loci(self):
        l = []
        for t, d in self.alignments.iteritems():
            for locus in d.iterkeys():
                if not locus in l:
                    l.append(locus)
        return l

    loci = property(_get_loci)

    def _get_number_of_taxa(self):
        return len(self.taxa)

    npairs = property(_get_number_of_taxa)

    def get_sample_table_string(self):    
        return '\n'.join(('BEGIN SAMPLE_TBL',
                '\n'.join((str(self.alignments[t][l]) for t, l in self._ordering)),
                'END SAMPLE_TBL'))

    def __str__(self):
        return self.get_sample_table_string()

    def __eq__(self, other):
        if not isinstance(other, SampleTable):
            return False
        if self.alignments != other.alignments:
            return False
        if self._ordering != self._ordering:
            return False
        return True
                                   

class AlignmentConfig(object):     

    def __init__(self, sample_table_row = None):
        self.taxon_name = None
        self.locus_name = None
        self.ploidy_multiplier = None
        self.mutation_rate_multiplier = None
        self.number_of_gene_copies = None
        self.kappa = None
        self.length = None
        self.base_frequencies = None
        self.path = None
        if sample_table_row:
            self._parse_sample_table_row(sample_table_row)

    def _parse_sample_table_row(self, sample_table_row):
        row = sample_table_row
        if isinstance(sample_table_row, str):
            row = [e.strip() for e in sample_table_row.split()]
        if len(row) != 12:
            raise errors.SampleTableRowError(
                    'sample table row has {0} columns'.format(len(row)))
        self.taxon_name, self.locus_name = row[0:2]
        self.ploidy_multiplier, self.mutation_rate_multiplier = (
                float(x) for x in row[2:4])
        self.number_of_gene_copies = tuple(int(x) for x in row[4:6])
        self.kappa, self.length = float(row[6]), int(row[7])
        self.base_frequencies = tuple(float(x) for x in row[8:11])
        self.path = row[11]

    def get_sample_table_row_elements(self):
        return(self.taxon_name,
                self.locus_name,
                self.ploidy_multiplier,
                self.ploidy_multiplier,
                self.mutation_rate_multiplier,
                self.number_of_gene_copies[0],
                self.number_of_gene_copies[1],
                self.kappa,
                self.base_frequencies[0],
                self.base_frequencies[1],
                self.base_frequencies[2],
                self.path)

    def get_sample_table_row_string(self):
        return '\t'.join([str(e) for e in self.get_sample_table_row_elements])
        
    def __str__(self):
        return self.get_sample_table_row_string()

