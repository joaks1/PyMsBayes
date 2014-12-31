#! /usr/bin/env python

"""
CLI for estimating the prior probabilities of divergence models.
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
    'version': 'Version 0.1.1',
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
    parser.add_argument('-d', '--dispersion-threshold',
            action = 'store',
            type = float,
            default = 0.01,
            help = ('The threshold for the dispersion index of divegence '
                    'times. The estimated prior probability that the '
                    'dispersion index is less than this threshold will '
                    'be reported for each config.'))
    parser.add_argument('-c', '--cv-threshold',
            action = 'store',
            type = float,
            default = 0.01,
            help = ('The threshold for the coefficient of variation (CV) of '
                    'divegence times. The estimated prior probability that the '
                    'CV is less than this threshold will '
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
            help = 'Run without verbose messaging.')
    parser.add_argument('--debug',
            action = 'store_true',
            help = 'Run in debugging mode.')

    args = parser.parse_args()

    ##########################################################################
    ## handle args

    from pymsbayes.utils.messaging import (LoggingControl,
            InfoLogger)

    LoggingControl.set_logging_level("INFO")
    if args.quiet:
        LoggingControl.set_logging_level("WARNING")
    if args.debug:
        LoggingControl.set_logging_level("DEBUG")
    log = LoggingControl.get_logger(__name__)

    from pymsbayes.teams import ModelProbabilityEstimatorTeam
    from pymsbayes.utils import GLOBAL_RNG

    if not args.seed:
        args.seed = random.randint(1, 999999999)
    log.info('Using seed {0}'.format(args.seed))
    GLOBAL_RNG.seed(args.seed)

    ##########################################################################
    ## begin analysis --- generate samples

    start_time = datetime.datetime.now()

    prob_esimator_team = ModelProbabilityEstimatorTeam(
            config_paths = args.configs,
            num_samples = args.num_prior_samples,
            omega_threshold = args.dispersion_threshold,
            cv_threshold = args.cv_threshold,
            num_processors = args.np)
    prob_esimator_team.start()

    for path in args.configs:
        sys.stdout.write('Prior probabilities for model {0}:\n'.format(
                path))
        for k, p in prob_esimator_team.psi_probs[path].iteritems():
            sys.stdout.write('\tnum of divergence events = {0}: {1}\n'.format(
                    k,
                    p))
        sys.stdout.write('\tdispersion of div times < {0}: {1}\n'.format(
                args.dispersion_threshold,
                prob_esimator_team.omega_probs[path]))
        sys.stdout.write('\tCV of div times < {0}: {1}\n'.format(
                args.cv_threshold,
                prob_esimator_team.cv_probs[path]))

    stop_time = datetime.datetime.now()
    log.info('[run_stats]')
    log.info('\tstart_time = {0}'.format(str(start_time)))
    log.info('\tstop_time = {0}'.format(str(stop_time)))
    log.info('\ttotal_duration = {0}'.format(str(stop_time - start_time)))

if __name__ == '__main__':
    main_cli()

