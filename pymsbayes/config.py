#! /usr/bin/env python

import os
import sys
import re
from cStringIO import StringIO
import ConfigParser

from pymsbayes.utils.messaging import get_logger
from pymsbayes.utils.functions import process_file_arg
from pymsbayes.utils.probability import (ContinuousUniformDistribution,
        BetaDistribution, DiscreteUniformDistribution)

_LOG = get_logger(__name__)

class MsBayesConfig(object):
    _begin_pattern = re.compile(r'^begin\s*sample_tbl$', re.IGNORECASE)
    _end_pattern = re.compile(r'^end\s*sample_tbl$', re.IGNORECASE)

    def __init__(self, cfg_file):
        self.theta = None
        self.d_theta = None,
        self.a_theta = None,
        self.tau = None,
        self.migration = None,
        self.recombination = None,
        self.psi = None,
        self.npairs = None
        self._parse_config(cfg_file)

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
        return preamble, table

    def _parse_table(self, table):
        taxa = set()
        for i, row in enumerate(table):
            taxa.add(row.strip().split()[0])
        self.npairs = len(taxa)

    def _set_priors(self, preamble):
        kwargs = self._parse_preamble(preamble)
        psi = int(kwargs['numtauclasses'])
        ltheta = float(kwargs.get('lowertheta', 0.0))
        utheta = float(kwargs['uppertheta'])
        anc_scalar = float(kwargs['upperancpopsize'])
        ltau = float(kwargs.get('lowertau', 0.0))
        utau = float(kwargs['uppertau'])
        lmig = float(kwargs.get('lowermig', 0.0))
        umig = float(kwargs.get('uppermig', 0.0))
        lrec = float(kwargs.get('lowerrec', 0.0))
        urec = float(kwargs.get('upperrec', 0.0))

        if psi == 0:
            self.psi = DiscreteUniformDistribution(1, self.npairs)
        else:
            self.psi = DiscreteUniformDistribution(psi, psi)
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
            d[tup[0]] = tup[1]
        return d

