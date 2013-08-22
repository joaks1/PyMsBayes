#! /usr/bin/env python

"""
Main CLI for PyMsBayes package.
"""

import os
import sys
import re
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
    parser.add_argument('-c', '--config',
            type = arg_is_config,
            required = True,
            help = ('msBayes config file to be used to generate saturation '
                    'plot.'))
    parser.add_argument('-n', '--num-prior-samples',
            action = 'store',
            type = int,
            default = 1000,
            help = ('The number of prior samples to simulate for the '
                    'saturation plot.'))
    parser.add_argument('--np',
            action = 'store',
            type = int,
            default = multiprocessing.cpu_count(),
            help = ('The maximum number of processes to run in parallel. The '
                    'default is the number of CPUs available on the machine.'))
    parser.add_argument('--output-dir',
            action = 'store',
            type = arg_is_dir,
            help = ('The directory in which all output files will be written. '
                    'The default is to use the directory of the first observed '
                    'config file.'))
    parser.add_argument('--temp-dir',
            action = 'store',
            type = arg_is_dir,
            help = ('A directory to temporarily stage files. The default is to '
                    'use the output directory.'))
    parser.add_argument('-s', '--stat-prefixes',
            nargs = '*',
            type = str,
            help = ('Prefixes of summary statistics to use in the analyses. '
                    'The prefixes should be separated by spaces. '
                    'Default: `-s pi. wattTheta. pi.net. tajD.denom.`.'))
    parser.add_argument('--compress',
            action = 'store_true',
            help = 'Compress large results files.')
    parser.add_argument('--keep-temps',
            action = 'store_true',
            help = 'Keep all temporary files.')
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

    from pymsbayes.workers import MsBayesWorker
    from pymsbayes.utils.parsing import (get_patterns_from_prefixes,
            DEFAULT_STAT_PATTERNS, get_stats_by_time)
    from pymsbayes.manager import Manager
    from pymsbayes.utils.tempfs import TempFileSystem
    from pymsbayes.utils import probability
    from pymsbayes.utils.functions import long_division
    from pymsbayes.config import MsBayesConfig
    from pymsbayes.utils import GLOBAL_RNG, MSBAYES_SORT_INDEX

    MSBAYES_SORT_INDEX.set_index(0)

    if not args.output_dir:
        args.output_dir = os.path.dirname(args.config)
    info = InfoLogger(os.path.join(args.output_dir, 'pymsbayes-info.txt'))

    if not args.temp_dir:
        args.temp_dir = args.output_dir
    temp_fs = TempFileSystem(parent=args.temp_dir, prefix='temp-files-')
    stat_patterns = DEFAULT_STAT_PATTERNS
    if args.stat_prefixes:
        for i in range(len(args.stat_prefixes)):
            if not args.stat_prefixes[i].endswith('.'):
                args.stat_prefixes[i] += '.'
        stat_patterns = get_patterns_from_prefixes(
                args.stat_prefixes,
                ignore_case=True)
    if not args.seed:
        args.seed = random.randint(1, 999999999)
    GLOBAL_RNG.seed(args.seed)

    cfg = MsBayesConfig(args.config)
    num_taxon_pairs = cfg.npairs
    cfg.div_model_prior = 'constrained'
    cfg.psi = probability.DiscreteUniformDistribution(num_taxon_pairs,
            num_taxon_pairs)
    config_path = temp_fs.get_file_path(prefix='cfg-')
    cfg.write(config_path)

    info.write('[pymsbayes]', log.info)
    info.write('\tprogram_name = {name}'.format(**_program_info), log.info)
    info.write('\tversion = {version}'.format(**_program_info), log.info)
    info.write('\toutput_directory = {0}'.format(args.output_dir), log.info)
    info.write('\ttemp_directory = {0}'.format(temp_fs.base_dir), log.info)
    info.write('\tsort_index = {0}'.format(
            MSBAYES_SORT_INDEX.current_value()), log.info)
    info.write('\tstat_patterns = {0}'.format(
            ', '.join([p.pattern for p in stat_patterns])), log.info)
    info.write('\tseed = {0}'.format(args.seed), log.info)
    info.write('\tnum_prior_samples = {0}'.format(args.num_prior_samples),
            log.info)

    info.write('\t[[config]]', log.debug)
    info.write('{0}'.format(str(cfg)), log.debug)

    ##########################################################################
    ## begin analysis --- generate samples

    start_time = datetime.datetime.now()

    if args.np > args.num_prior_samples:
        args.np = args.num_prior_samples
    batch_size, remainder = long_division(args.num_prior_samples, args.np)
    schema = 'abctoolbox'
    workers = []
    for i in range(args.np):
        sample_size = batch_size
        if i == (args.np - 1):
            sample_size += remainder
        w = MsBayesWorker(
                temp_fs = temp_fs,
                sample_size = sample_size,
                config_path = config_path,
                report_parameters = True,
                schema = schema,
                include_header = True,
                stat_patterns = stat_patterns,
                write_stats_file = False)
        workers.append(w)

    log.info('Generating samples...')
    workers = Manager.run_workers(
            workers = workers,
            num_processors = args.np)
    log.info('Parsing samples...')
    stats_by_time = get_stats_by_time([w.prior_path for w in workers])
    sys.stderr.write('{0}\n'.format(stats_by_time))

    stop_time = datetime.datetime.now()
    log.info('Done!')
    info.write('\t[[run_stats]]', log.info)
    info.write('\t\tstart_time = {0}'.format(str(start_time)), log.info)
    info.write('\t\tstop_time = {0}'.format(str(stop_time)), log.info)
    info.write('\t\ttotal_duration = {0}'.format(str(stop_time - start_time)),
            log.info)

    if not args.keep_temps:
        log.debug('purging temps...')
        temp_fs.purge()

if __name__ == '__main__':
    main_cli()

