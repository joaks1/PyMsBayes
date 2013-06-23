#! /usr/bin/env python

"""
Main CLI for PyMsBayes package.
"""

import os
import sys
import re
import multiprocessing
import random
import argparse
import datetime
import logging

from pymsbayes.fileio import expand_path, process_file_arg, open
from pymsbayes.utils import WORK_FORCE, GLOBAL_RNG
from pymsbayes.utils.messaging import get_logger, LOGGING_LEVEL_ENV_VAR
from pymsbayes.utils.functions import (is_file, is_dir, long_division,
        mk_new_dir)

_LOG = get_logger(__name__)

_program_info = {
    'name': os.path.basename(__file__),
    'author': 'Jamie Oaks',
    'version': 'Version 0.1.0',
    'description': __doc__,
    'copyright': 'Copyright (C) 2013 Jamie Oaks',
    'license': 'GNU GPL version 3 or later',}

class InfoLogger(object):
    def __init__(self, path):
        self.path = path

    def write(self, msg):
        out = open(self.path, 'a')
        out.write(msg)
        out.close()

def arg_is_file(path):
    try:
        if not is_file(path):
            raise
    except:
        msg = '{0!r} is not a file'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return expand_path(path)

def arg_is_dir(path):
    try:
        if not is_dir(path):
            raise
    except:
        msg = '{0!r} is not a directory'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return expand_path(path)

def main_cli():
    description = '{name} {version}'.format(**_program_info)
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('-o', '--observed-config',
            action = 'store',
            type = arg_is_file,
            required = True,
            help = ('The msBayes config file that specifies the observed '
                    'data. If used in combination with `-r` this config will '
                    'be used to simulate pseudo-observed data. If analyzing '
                    'real data, do not use the `-r` option, and the data '
                    'files specified in the config must exist and contain the '
                    'sequence data.'))
    parser.add_argument('-p', '--prior-configs',
            nargs = '+',
            type = arg_is_file,
            required = True,
            help = ('One or more config files to be used to generate prior '
                    'samples. If more than one config is specified, they '
                    'should be separated by spaces.'))
    parser.add_argument('-r', '--reps',
            action = 'store',
            type = int,
            default = 0,
            help = ('This option has two effects. First, it signifies that '
                    'the analysis will be simulation based (i.e., no real '
                    'data will be used). Second, it specifies how many '
                    'simulation replicates to perform (i.e., how many data '
                    'sets to simulate and analyze).'))
    parser.add_argument('-n', '--num-prior-samples',
            action = 'store',
            type = int,
            default = 1000000,
            help = ('The number of prior samples to simulate for each prior '
                    'config specified with `-p`.'))
    parser.add_argument('--prior-batch-size',
            action = 'store',
            type = int,
            default = 10000,
            help = ('The number of prior samples to simulate for each batch.'))
    parser.add_argument('--num-posterior-samples',
            action = 'store',
            type = int,
            default = 1000,
            help = ('The number of posterior samples desired for each '
                    'analysis. Default: 1000.'))
    parser.add_argument('--num-standardizing-samples',
            action = 'store',
            type = int,
            default = 10000,
            help = ('The number of prior samples desired to use for '
                    'standardizing statistics. Default: 10000.'))
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
                    'The default is to use the directory of the observed '
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
                    'Default: `-s pi wattTheta pi.net tajD.denom`.'))
    parser.add_argument('-b', '--bandwidth',
            action = 'store',
            type = float,
            help = ('Smoothing parameter for the posterior kernal density '
                    'estimation. This option is used for the `glm` '
                    'regression method. The default is 2 / '
                    '`num-posterior-samples`.'))
    parser.add_argument('-q', '--num-posterior-quantiles',
            action = 'store',
            type = int,
            default = 1000,
            help = ('The number of equally spaced quantiles at which to '
                    'evaluate the GLM-estimated posterior density. '
                    'Default: 1000.'))
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

    _LOG.setLevel(logging.INFO)
    os.environ[LOGGING_LEVEL_ENV_VAR] = "INFO"
    if args.quiet:
        _LOG.setLevel(logging.WARNING)
        os.environ[LOGGING_LEVEL_ENV_VAR] = "WARNING"
    if args.debug:
        _LOG.setLevel(logging.DEBUG)
        os.environ[LOGGING_LEVEL_ENV_VAR] = "DEBUG"

    from pymsbayes.workers import (MsBayesWorker, merge_prior_files,
            ObsSumStatsWorker)
    from pymsbayes.teams import ABCTeam
    from pymsbayes.config import MsBayesConfig
    from pymsbayes.utils.parsing import (get_patterns_from_prefixes,
        DEFAULT_STAT_PATTERNS, DIV_MODEL_PATTERNS, MODEL_PATTERNS, PSI_PATTERNS,
        MEAN_TAU_PATTERNS, OMEGA_PATTERNS)
    from pymsbayes.manager import Manager
    from pymsbayes.utils.tempfs import TempFileSystem

    if not args.output_dir:
        args.output_dir = os.path.dirname(args.observed_config)
    base_dir = mk_new_dir(os.path.join(args.output_dir, 'pymsbayes-results'))
    if not args.temp_dir:
        args.temp_dir = base_dir
    info = InfoLogger(os.path.join(base_dir, 'pymsbayes-info.txt'))
    info.write('[pymsbayes]\n'.format(base_dir))
    info.write('\tversion = {version}\n'.format(**_program_info))
    info.write('\toutput_directory = {0}\n'.format(base_dir))
    temp_fs = TempFileSystem(parent=args.temp_dir, prefix='temp-files-')
    base_temp_dir = temp_fs.base_dir
    info.write('\ttemp_directory = {0}\n'.format(base_temp_dir))
    if args.reps < 1:
        info.write('\tsimulate_data = False\n')
    else:
        info.write('\tsimulate_data = True\n')
    if len(args.prior_configs) != len(set(args.prior_configs)):
        raise ValueError('All paths to prior config files must be unique') 
    stat_patterns = DEFAULT_STAT_PATTERNS
    if args.stat_prefixes:
        stat_patterns = get_patterns_from_prefixes(
                args.stat_prefixes,
                ignore_case=True)
    parameter_patterns = DIV_MODEL_PATTERNS + MODEL_PATTERNS + \
            PSI_PATTERNS + MEAN_TAU_PATTERNS + OMEGA_PATTERNS
    if not args.bandwidth:
        args.bandwidth = 2 / float(args.num_posterior_samples)
    if not args.seed:
        args.seed = random.randint(1, 999999999)
    GLOBAL_RNG.seed(args.seed)
    observed_path = os.path.join(base_dir, 'observed.txt')
    info.write('\tseed = {0}\n'.format(args.seed))
    info.write('\tnum_processors = {0}\n'.format(args.np))
    info.write('\tbandwidth = {0}\n'.format(args.bandwidth))
    info.write('\tposterior_quantiles = {0}\n'.format(
            args.num_posterior_quantiles))
    info.write('\tposterior_sample_size = {0}\n'.format(
            args.num_posterior_samples))
    info.write('\tnum_standardizing_samples = {0}\n'.format(
            args.num_posterior_samples))
    info.write('\tobserved_path = {0}\n'.format(observed_path))
    info.write('\t[[column_prefixes]]\n')
    info.write('\t\tstat_patterns = {0}\n'.format(
            ', '.join([p.pattern for p in stat_patterns])))
    info.write('\t\tparameter_patterns = {0}\n'.format(
            ', '.join([p.pattern for p in parameter_patterns])))

    models_to_configs = {}
    configs_to_models = {}
    num_taxon_pairs = None
    for i in range(len(args.prior_configs)):
        model_idx = i + 1
        models_to_configs[model_idx] = args.prior_configs[i]
        configs_to_models[args.prior_configs[i]] = model_idx
        cfg = MsBayesConfig(args.prior_configs[i])
        if not num_taxon_pairs:
            num_taxon_pairs = cfg.npairs
        else:
            assert num_taxon_pairs == cfg.npairs
    model_indices = models_to_configs.keys()
    info.write('\t[[prior_configs]]\n')
    for model_idx, cfg in models_to_configs.iteritems():
        info.write('\t\t{0} = {1}\n'.format(model_idx, cfg))

    ##########################################################################
    ## begin analysis --- get observed summary stats

    start_time = datetime.datetime.now()

    observed_temp_fs = TempFileSystem(parent = base_temp_dir,
            prefix = 'observed-temps-')

    if args.reps < 1:
        ss_worker = ObsSumStatsWorker(
                temp_fs = observed_temp_fs,
                config_path = args.observed_config,
                output_path = observed_path,
                schema = 'abctoolbox',
                stat_patterns = stat_patterns)
        ss_worker.start()
    else:
        num_observed_workers = min([args.reps, args.np])
        if args.reps <= args.np:
            observed_batch_size = 1
        else:
            observed_batch_size, remainder = long_division(args.reps, args.np)
        observed_model_idx = configs_to_models.get(args.observed_config, 0)
        info.write('\t[[observed_configs]]\n')
        info.write('\t\t{0} = {1}\n'.format(observed_model_idx,
            args.observed_config))
        schema = 'abctoolbox'
        msbayes_workers = []
        for i in range(num_observed_workers):
            worker = MsBayesWorker(
                    temp_fs = observed_temp_fs,
                    sample_size = observed_batch_size,
                    config_path = args.observed_config,
                    model_index = observed_model_idx,
                    report_parameters = True,
                    schema = schema,
                    include_header = True,
                    stat_patterns = stat_patterns,
                    write_stats_file = True,
                    staging_dir = None)
            WORK_FORCE.put(worker)
            msbayes_workers.append(worker)
        if remainder > 0:
            worker = MsBayesWorker(
                    temp_fs = observed_temp_fs,
                    sample_size = remainder,
                    config_path = args.observed_config,
                    model_index = observed_model_idx,
                    report_parameters = True,
                    schema = schema,
                    include_header = True,
                    stat_patterns = stat_patterns,
                    write_stats_file = True,
                    staging_dir = None)
            WORK_FORCE.put(worker)
            msbayes_workers.append(worker)

        # run parallel msbayes processes
        msbayes_result_queue = multiprocessing.Queue()
        msbayes_managers = []
        for i in range(args.np):
            m = Manager(work_queue = WORK_FORCE,
                    result_queue = msbayes_result_queue)
            m.start()
            msbayes_managers.append(m)
        for i in range(len(msbayes_workers)):
            msbayes_workers[i] = msbayes_result_queue.get()
        for m in msbayes_managers:
            m.join()
        assert WORK_FORCE.empty()
        assert msbayes_result_queue.empty()

        # merged simulated observed data into one file
        merge_prior_files([w.prior_stats_path for w in msbayes_workers],
                observed_path)
        lc = line_count(observed_path, ignore_headers=True)
        if lc != args.reps:
            raise Exception('The number of observed simulations ({0}) does '
                    'not match the number of reps ({1})'.format(lc,
                            args.reps))
    if not args.keep_temps:
        _LOG.debug('purging observed temps...')
        observed_temp_fs.purge()

    ##########################################################################
    ## Begin ABC analyses

    abc_team = ABCTeam(
            temp_fs = temp_fs,
            observed_stats_file = observed_path,
            num_taxon_pairs = num_taxon_pairs,
            model_indices_to_config_paths = models_to_configs,
            num_prior_samples = args.num_prior_samples,
            num_processors = args.np,
            num_standardizing_samples = args.num_standardizing_samples,
            num_posterior_samples = args.num_posterior_samples,
            num_posterior_density_quantiles = args.num_posterior_quantiles,
            batch_size = args.prior_batch_size,
            output_dir = base_dir,
            output_prefix = '',
            rng = GLOBAL_RNG,
            sort_index = None,
            report_parameters = True,
            stat_patterns = stat_patterns,
            abctoolbox_bandwidth = args.bandwidth,
            omega_threshold = 0.01,
            compress = True,
            keep_temps = args.keep_temps,
            global_estimate_only = False,
            work_queue = WORK_FORCE)
    abc_team.run()

    stop_time = datetime.datetime.now()
    info.write('\t[[run_stats]]\n')
    info.write('\t\tstart_time = {0}\n'.format(str(start_time)))
    info.write('\t\tstop_time = {0}\n'.format(str(stop_time)))
    info.write('\t\ttotal_duration = {0}\n'.format(str(stop_time - start_time)))

    if not args.keep_temps:
        _LOG.debug('purging temps...')
        temp_fs.purge()

if __name__ == '__main__':
    main_cli()

