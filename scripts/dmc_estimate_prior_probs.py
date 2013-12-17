#! /usr/bin/env python

"""
Main CLI for PyMsBayes package.
"""

import os
import sys
import re
import math
import glob
import multiprocessing
import random
import argparse
import datetime
import logging

from pymsbayes.utils.argparse_utils import arg_is_dir, arg_is_config

_program_info = {
    'name': os.path.basename(__file__),
    'author': 'Jamie Oaks',
    'version': 'Version 0.1.0',
    'description': __doc__,
    'copyright': 'Copyright (C) 2013 Jamie Oaks',
    'license': 'GNU GPL version 3 or later',}

def main_cli():
    description = '{name} {version}'.format(**_program_info)
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('configs', metavar = 'CONFIG-PATH',
            type = arg_is_config,
            nargs = '+',
            help = ('msBayes config file paths for which to estimate prior '
                    'probabilities.'))
    parser.add_argument('-n', '--num-prior-samples',
            action = 'store',
            type = int,
            default = 1000,
            help = ('The number of prior samples to simulate for '
                    'proabability estimates.'))
    parser.add_argument('--np',
            action = 'store',
            type = int,
            default = multiprocessing.cpu_count(),
            help = ('The maximum number of processes to run in parallel. The '
                    'default is the number of CPUs available on the machine.'))
    parser.add_argument('-e', '--num-divergence-events',
            action = 'store',
            type = int,
            default = 1,
            help = ('The number of divergence events. The prior probability '
                    'for this number of events will be reported for each '
                    'config.'))
    parser.add_argument('-d', '--dispersion-threshold',
            action = 'store',
            type = float,
            default = 0.01,
            help = ('The threshold for the dispersion index of divegence '
                    'times. The estimated prior probability that the '
                    'dispersion index is less than this threshold will '
                    'be reported for each config.'))
    parser.add_argument('--seed',
            action = 'store',
            type = int,
            help = 'Random number seed to use for the analysis.')
    parser.add_argument('--version',
            action = 'version',
            version = '%(prog)s ' + _program_info['version'],
            help = 'Report version and exit.')
    parser.add_argument('--quiet',
            action = 'store_true',
            help = 'Run with verbose messaging.')
    parser.add_argument('--debug',
            action = 'store_true',
            help = 'Run in debugging mode.')

    args = parser.parse_args()

    ##########################################################################
    ## handle args

    from pymsbayes.utils.messaging import (get_logger, LOGGING_LEVEL_ENV_VAR,
            InfoLogger)

    os.environ[LOGGING_LEVEL_ENV_VAR] = "INFO"
    if args.quiet:
        os.environ[LOGGING_LEVEL_ENV_VAR] = "WARNING"
    if args.debug:
        os.environ[LOGGING_LEVEL_ENV_VAR] = "DEBUG"
    log = get_logger()

    from pymsbayes.workers import ModelProbabilityEstimator
    from pymsbayes.manager import Manager
    from pymsbayes.utils import stats, GLOBAL_RNG
    from pymsbayes.utils.functions import long_division
    from pymsbayes.config import MsBayesConfig

    if not args.seed:
        args.seed = random.randint(1, 999999999)
    GLOBAL_RNG.seed(args.seed)

    configs = dict(zip(args.configs,
            [MsBayesConfig(c) for c in args.configs]))
    psi_summaries = dict(zip(args.configs,
            [stats.SampleSummary() for c in args.configs]))
    omega_summaries = dict(zip(args.configs,
            [stats.SampleSummary() for c in args.configs]))

    ##########################################################################
    ## begin analysis --- generate samples

    start_time = datetime.datetime.now()

    if args.np > args.num_prior_samples:
        args.np = args.num_prior_samples
    batch_size, remainder = long_division(args.num_prior_samples, args.np)
    workers = []
    for path, cfg in configs.iteritems():
        for i in range(args.np):
            sample_size = batch_size
            if i == (args.np - 1):
                sample_size += remainder
            w = ModelProbabilityEstimator(
                    config = cfg,
                    num_samples = sample_size,
                    psi_of_interest = args.num_divergence_events,
                    omega_threshold = args.dispersion_threshold,
                    tag = path)
            workers.append(w)

    log.info('Using seed {0}'.format(args.seed))
    log.info('Generating samples...')
    workers = Manager.run_workers(
            workers = workers,
            num_processors = args.np)
    log.info('Done!')
    log.info('Summarizing results...')
    for w in workers:
        psi_summaries[w.tag].update(w.psi_summary)
        omega_summaries[w.tag].update(w.omega_summary)
    for path in args.configs:
        sys.stdout.write('Prior probabilities for model {0}:\n'.format(
                path))
        sys.stdout.write('\tnum of divergence events = {0}: {1}\n'.format(
                args.num_divergence_events,
                psi_summaries[path].mean))
        sys.stdout.write('\tdispersion of div times < {0}: {1}\n'.format(
                args.dispersion_threshold,
                omega_summaries[path].mean))

    stop_time = datetime.datetime.now()
    log.info('[run_stats]')
    log.info('\tstart_time = {0}'.format(str(start_time)))
    log.info('\tstop_time = {0}'.format(str(stop_time)))
    log.info('\ttotal_duration = {0}'.format(str(stop_time - start_time)))

if __name__ == '__main__':
    main_cli()

