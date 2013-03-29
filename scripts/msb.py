#! /usr/bin/env python

"""
Main CLI for PyMsBayes package.
"""

import os
import sys
import re
import multiprocessing
import argparse

from pymsbayes.workers import (MsBayesWorker, MsRejectWorker, RegressionWorker,
        merge_priors, assemble_msreject_workers, get_parameter_indices,
        get_stat_indices, parse_header, get_patterns_from_prefixes)
from pymsbayes.workers import (DEFAULT_STAT_PATTERNS, PSI_PATTERNS,
        MODEL_PATTERNS, MEAN_TAU_PATTERNS, OMEGA_PATTERNS)
from pymsbayes.manager import Manager
from pymsbayes.utils import WORK_FORCE, GLOBAL_RNG
from pymsbayes.utils.functions import (is_file, is_dir, expand_path,
        long_division, mk_new_dir, line_count)
from pymsbayes.utils.tempfs import TempFileSystem
from pymsbayes.utils.messaging import get_logger

_program_info = {
    'name': os.path.basename(__file__),
    'author': 'Jamie Oaks',
    'version': 'Version 0.1.0',
    'description': __doc__,
    'copyright': 'Copyright (C) 2013 Jamie Oaks',
    'license': 'GNU GPL version 3 or later',}

_LOG = get_logger(__name__, 'INFO')


def arg_is_file(path):
    if not is_file(path):
        msg = '{0!r} is not a file'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return expand_path(path)

def arg_is_dir(path):
    if not is_dir(path):
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
    parser.add_argument('-n', '--nsamples',
            action = 'store',
            type = int,
            default = 1000000,
            help = ('The number of prior samples to simulate for each prior '
                    'config specified with `-p`.'))
    parser.add_argument('--np',
            action = 'store',
            type = int,
            default = multiprocessing.cpu_count(),
            help = ('The maximum number of processes to run in parallel. The '
                    'default is the number of CPUs available on the machine.')
    parser.add_argument('--output-dir',
            action = 'store',
            type = arg_is_dir,
            help = ('The directory in which all output files will be written. '
                    'The default is to use the directory of the observed '
                    'config file.'))
    parser.add_argument('-s', '--stat-prefixes',
            nargs = '*',
            type = str,
            help = ('Prefixes of summary statistics to use in the analyses. '
                    'The prefixes should be separated by spaces. '
                    'Default: `-s pi wattTheta pi.net tajD.denom`.'))
    parser.add_argument('-c', '--continuous-prefixes',
            nargs = '*',
            type = str,
            help = ('Prefixes of continuous parameters on which regression '
                    'adjustment should be performed. '
                    'The prefixes should be separated by spaces. '
                    'Default: `-c PRI.E.t PRI.omega`.'))
    parser.add_argument('-d', '--discrete-prefixes',
            nargs = '*',
            type = str,
            help = ('Prefixes of discrete parameters on which regression '
                    'adjustment should be performed. '
                    'The prefixes should be separated by spaces. '
                    'Default: `-c PRI.model PRI.Psi`.'))
    parser.add_argument('-m', '--merge-priors',
            action = 'store_true',
            default = 'False',
            help = ('Merge the priors generated from the specified configs '
                    'together before analysis. Otherwise, the observed data '
                    'sets will be analyzed under each prior separately.'))
    parser.add_argument('--repress-parameters',
            action = 'store_false',
            dest = 'report_parameters',
            help = ('Do not report divergence time and population size '
                    'parameters (i.e. only report the summaries of '
                    'parameters [E(t), Var(t), omega, psi] as in the standard '
                    'msBayes.'))
    parser.add_argument('-k', '--keep-priors',
            action = 'store_true',
            default = 'False',
            help = 'Keep prior-sample files.')
    parser.add_argument('--keep-temps',
            action = 'store_true',
            default = 'False',
            help = 'Keep all temporary files.')
    parser.add_argument('--seed',
            action = 'store',
            type = int,
            help = 'Random number seed to use for the analysis.')
    parser.add_argument('--dry-run',
            action = 'store_true',
            default = 'False',
            help = 'Report configuration of the analysis and exit.')
    parser.add_argument('--version',
            action = 'version',
            version = '%(prog)s ' + _program_info['version'],
            help = 'Report version and exit.')
    parser.add_argument('-v', '--verbose',
            action = 'store_true',
            default = 'False',
            help = 'Run with verbose messaging.')
    parser.add_argument('--debug',
            action = 'store_true',
            default = 'False',
            help = 'Run in debugging mode.')

    args = parser.parse_args()

    if not args.output_dir:
        args.output_dir = os.path.dirname(args.observed_config)
    base_dir = mk_new_dir(os.path.join(args.output_dir, 'pymsbayes-output'))
    base_temp_dir = mk_new_dir(os.path.join(base_dir, 'temp-files'))
    if not args.reps < 1:
        raise NotImplementedError('Sorry, analyzing real data is not yet '
                'implemented')
    if len(args.prior_configs) != len(set(args.prior_configs)):
        raise ValueError('All paths to prior config files must be unique') 
    stat_patterns = DEFAULT_STAT_PATTERNS
    if args.stat_prefixes:
        stat_patterns = get_patterns_from_prefixes(
                args.stat_prefixes,
                ignore_case=True)
    continuous_patterns = OMEGA_PATTERNS + MEAN_TAU_PATTERNS
    if args.continuous_prefixes:
        continuous_patterns = get_patterns_from_prefixes(
                args.continuous_prefixes,
                ignore_case=True)
    discrete_patterns = PSI_PATTERNS + MODEL_PATTERNS
    if args.discrete_prefixes:
        discrete_patterns = get_patterns_from_prefixes(
                args.discrete_prefixes,
                ignore_case=True)
    if args.seed:
        GLOBAL_RNG.seed(args.seed)

    prior_subsample_size = 100000
    total_nsamples = args.nsamples * len(args.prior_configs)
    if args.reps > 0:
        total_nsamples += args.reps
    max_proc_samples = total_nsamples / args.np
    if max_proc_samples < prior_subsample_size:
        prior_subsample_size = max_proc_samples
    num_prior_workers, remainder_samples = long_division(args.nsamples,
            prior_subsample_size)
    num_observed_workers, remainder_observed_reps = long_division(args.reps,
            prior_subsample_size)

    working_prior_temp_fs = TempFileSystem(parent = base_temp_dir,
            prefix = 'pymsbayes-working-prior-files-')
    working_observed_temp_fs = TempFileSystem(parent = base_temp_dir,
            prefix = 'pymsbayes-working-observed-files-')
    models_to_configs = {}
    configs_to_models = {}
    msbayes_workers = []
    # put prior-generating msbayes workers in the queue
    for i in range(len(args.prior_configs)):
        model_idx = i + 1
        models_to_configs[model_idx] = args.prior_configs[i]
        configs_to_models[args.prior_configs[i]] = model_idx
        for j in range(num_prior_workers):
            worker = MsBayeWorker(
                    temp_fs = working_prior_temp_fs,
                    sample_size = prior_subsample_size,
                    config_path = args.prior_configs[i],
                    model_index = model_idx,
                    report_parameters = args.report_parameters,
                    include_header = False,
                    stat_patterns = stat_patterns)
            WORK_FORCE.put(worker)
            msbayes_workers.append(worker)
        if remainder_samples > 0:
            worker = MsBayeWorker(
                    temp_fs = working_prior_temp_fs,
                    sample_size = remainder_samples,
                    config_path = args.prior_configs[i],
                    model_index = model_idx,
                    report_parameters = args.report_parameters,
                    include_header = False,
                    stat_patterns = stat_patterns)
            WORK_FORCE.put(worker)
            msbayes_workers.append(worker)

    # put observed-data-generating msbayes workers in the queue
    observed_model_idx = configs_to_models.get(args.observed_config, 0)
    for i in range(num_observed_workers):
        worker = MsBayeWorker(
                temp_fs = working_observed_temp_fs,
                sample_size = prior_subsample_size,
                config_path = args.observed_config,
                model_index = observed_model_idx,
                report_parameters = args.report_parameters,
                include_header = True,
                stat_patterns = stat_patterns)
        WORK_FORCE.put(worker)
        msbayes_workers.append(worker)
    if remainder_observed_reps > 0:
        worker = MsBayeWorker(
                temp_fs = working_observed_temp_fs,
                sample_size = remainder_observed_reps,
                config_path = args.observed_config,
                model_index = observed_model_idx,
                report_parameters = args.report_parameters,
                include_header = False,
                stat_patterns = stat_patterns)
        WORK_FORCE.put(worker)
        msbayes_workers.append(worker)

    # run parallel processes
    msbayes_result_queue = multiprocessing.Queue()
    msbayes_managers = []
    for i in range(args.np):
        m = Manager(work_queue = WORK_FORCE,
                result_queue = msbayes_result_queue)
        m.start()
        msbayes_managers.append(m)
    for i in len(msbayes_workers):
        msbayes_workers[i] = msbayes_result_queue.get()
    for m in msbayes_managers:
        m.join()
    assert WORK_FORCE.empty()
    assert msbayes_result_queue.empty()

    msbayes_observed_workers = []
    msbayes_prior_workers = {}
    for w in msbayes_workers:
        if w.include_header:
            msbayes_observed_workers.append(w)
        else:
            if w.model_index in msbayes_prior_workers.keys():
                msbayes_prior_workers[w.model_index].append(w)
            else:
                msbayes_prior_workers[w.model_index] = [w]
    assert len(msbayes_prior_workers) == len(args.prior_configs)
    
    observed_path = os.path.join(base_dir, 'observed.txt')
    if args.reps > 0:
        obs_header_path = working_observed_temp_fs.get_file_path(
                prefix = 'observed-header-')
        merge_priors(workers = msbayes_observed_workers,
                prior_path = observed_path,
                header_path = obs_header_path,
                include_header = True)
        lc = line_count(observed_path)
        if lc != args.reps + 1:
            raise Exception('The number of observed simulations ({0}) does '
                    'not match the number of reps ({1})'.format(lc-1,
                            args.reps))
        if not args.keep_temps:
            working_observed_temp_fs.purge()
    prior_temp_fs = TempFileSystem(parent = base_temp_dir,
            prefix = 'pymsbayes-prior-files-')
    prior_paths = {}
    prior_paths['header'] = prior_temp_fs.get_file_path(
            prefix = 'prior-header-')
    for mod_idx in msbayes_prior_workers.iterkeys():
        prior_path = prior_temp_fs.get_file_path(
                prefix = 'prior-{0}-{1}-'.format(mod_idx, args.nsamples))
        merge_priors(workers = msbayes_prior_workers[mod_idx],
                prior_path = prior_path,
                header_path = prior_paths['header'],
                include_header = False)
        lc = line_count(prior_path)
        if lc != args.nsamples:
            raise Exception('The number of prior samples ({0}) for model '
                    '{1} does not match nsamples ({2})'.format(
                            lc, model_idx, args.nsamples))
    if not args.keep_temps:
        working_prior_temp_fs.purge()
    if args.merge_priors:

                            


def merge_priors(workers, prior_path, header_path=None, include_header=False):

if __name__ == '__main__':
    main_cli()

