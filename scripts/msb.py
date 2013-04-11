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

from pymsbayes.fileio import expand_path, process_file_arg, open
from pymsbayes.workers import (MsBayesWorker, MsRejectWorker, RegressionWorker,
        merge_priors, assemble_rejection_workers, get_parameter_indices,
        get_stat_indices, parse_header, get_patterns_from_prefixes,
        merge_prior_files, line_count)
from pymsbayes.workers import (DEFAULT_STAT_PATTERNS, PSI_PATTERNS,
        MODEL_PATTERNS, MEAN_TAU_PATTERNS, OMEGA_PATTERNS)
from pymsbayes.manager import Manager
from pymsbayes.utils import WORK_FORCE, GLOBAL_RNG
from pymsbayes.utils.functions import (is_file, is_dir, long_division,
        mk_new_dir, get_tolerance)
from pymsbayes.utils.tempfs import TempFileSystem
from pymsbayes.utils.messaging import get_logger

_program_info = {
    'name': os.path.basename(__file__),
    'author': 'Jamie Oaks',
    'version': 'Version 0.1.0',
    'description': __doc__,
    'copyright': 'Copyright (C) 2013 Jamie Oaks',
    'license': 'GNU GPL version 3 or later',}

_LOG = get_logger(__name__)

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
    parser.add_argument('--num-posterior-samples',
            action = 'store',
            type = int,
            default = 1000,
            help = ('The number of posterior samples desired for each '
                    'analysis. Default: 1000.'))
    parser.add_argument('--np',
            action = 'store',
            type = int,
            default = multiprocessing.cpu_count(),
            help = ('The maximum number of processes to run in parallel. The '
                    'default is the number of CPUs available on the machine.'))
    parser.add_argument('--rejection-tool',
            action = 'store',
            type = str,
            choices = ['msreject', 'abctoolbox'],
            default = 'abctoolbox',
            help = ('The tool to use for accepting posterior samples. '
                    'The default is `abctoolbox`.'))
    parser.add_argument('--regression-method',
            action = 'store',
            type = str,
            choices = ['llr', 'glm'],
            default = 'glm',
            help = ('The regression method to use for post-rejection '
                    'adjustment of the posterior sample: '
                    'local-linear/multinomial-logistic regression (`llr`) '
                    'or general-linear model regression (`glm`; Default).'))
    parser.add_argument('--output-dir',
            action = 'store',
            type = arg_is_dir,
            help = ('The directory in which all output files will be written. '
                    'The default is to use the directory of the observed '
                    'config file.'))
    parser.add_argument('--staging-dir',
            action = 'store',
            type = arg_is_dir,
            help = ('A directory to temporarily stage prior files. This option '
                    'can be useful on clusters to speed up I/O while '
                    'generating prior samples. You can designate a local temp '
                    'directory on a compute node to avoid constant writing to '
                    'a shared drive.'))
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
    parser.add_argument('--no-regression',
            action = 'store_false',
            dest = 'regression',
            help = 'Do not perform regression on posterior.')
    parser.add_argument('-b', '--bandwidth',
            action = 'store',
            type = float,
            help = ('Smoothing parameter for the posterior kernal density '
                    'estimation. This option is only used if `glm` is the '
                    'regression method. The default is 2 / '
                    '`num-posterior-samples`.'))
    parser.add_argument('-q', '--num-posterior-quantiles',
            action = 'store',
            type = int,
            default = 1000,
            help = ('The number of equally spaced quantiles at which to '
                    'evaluate the GLM-estimated posterior density. '
                    'This option is only used if `glm` is the regression '
                    'method. Default: 1000.'))
    parser.add_argument('-k', '--keep-priors',
            action = 'store_true',
            help = 'Keep prior-sample files.')
    parser.add_argument('--keep-temps',
            action = 'store_true',
            help = 'Keep all temporary files.')
    parser.add_argument('--seed',
            action = 'store',
            type = int,
            help = 'Random number seed to use for the analysis.')
    # parser.add_argument('--dry-run',
    #         action = 'store_true',
    #         help = 'Report configuration of the analysis and exit.')
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

    if args.quiet:
        _LOG.setlevel(logging.WARNING)
    if args.debug:
        _LOG.setlevel(logging.DEBUG)
    if not args.output_dir:
        args.output_dir = os.path.dirname(args.observed_config)
    staging_dir = None
    if args.staging_dir:
        staging_dir = args.staging_dir
    base_dir = mk_new_dir(os.path.join(args.output_dir, 'pymsbayes-output'))
    info = InfoLogger(os.path.join(base_dir, 'pymsbayes-info.txt'))
    info.write('[pymsbayes]\n'.format(base_dir))
    info.write('\toutput_directory = {0}\n'.format(base_dir))
    base_temp_dir = mk_new_dir(os.path.join(base_dir, 'temp-files'))
    info.write('\ttemp_directory = {0}\n'.format(base_temp_dir))
    if args.reps < 1:
        info.write('\tsimulate_data = False\n')
        raise NotImplementedError('Sorry, analyzing real data is not yet '
                'implemented')
    else:
        info.write('\tsimulate_data = True\n')
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
    if not args.bandwidth:
        args.bandwidth = 2 / float(args.num_posterior_samples)
    if not args.seed:
        args.seed = random.randint(1, 999999999)
    GLOBAL_RNG.seed(args.seed)
    info.write('\tseed = {0}\n'.format(args.seed))
    info.write('\tnum_processors = {0}\n'.format(args.np))
    info.write('\t[[column_prefixes]]\n')
    info.write('\t\tstat_patterns = {0}\n'.format(
            ', '.join([p.pattern for p in stat_patterns])))
    info.write('\t\tcontinuous_patterns = {0}\n'.format(
            ', '.join([p.pattern for p in continuous_patterns])))
    info.write('\t\tdiscrete_patterns = {0}\n'.format(
            ', '.join([p.pattern for p in discrete_patterns])))

    # calculate decent prior chunk size from user settings
    prior_subsample_size = 20000
    total_num_prior_samples = args.num_prior_samples * len(args.prior_configs)
    if args.reps > 0:
        total_num_prior_samples += args.reps
    max_proc_samples = total_num_prior_samples / args.np
    if max_proc_samples < prior_subsample_size:
        prior_subsample_size = max_proc_samples
    num_prior_workers, remainder_samples = long_division(args.num_prior_samples,
            prior_subsample_size)
    num_observed_workers, remainder_observed_reps = long_division(args.reps,
            prior_subsample_size)

    ##########################################################################
    ## begin analysis --- generate priors and observed simulations

    start_time = datetime.datetime.now()
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
            worker = MsBayesWorker(
                    temp_fs = working_prior_temp_fs,
                    sample_size = prior_subsample_size,
                    config_path = args.prior_configs[i],
                    model_index = model_idx,
                    report_parameters = args.report_parameters,
                    schema = args.rejection_tool,
                    include_header = False,
                    stat_patterns = stat_patterns,
                    staging_dir = staging_dir)
            WORK_FORCE.put(worker)
            msbayes_workers.append(worker)
        if remainder_samples > 0:
            worker = MsBayesWorker(
                    temp_fs = working_prior_temp_fs,
                    sample_size = remainder_samples,
                    config_path = args.prior_configs[i],
                    model_index = model_idx,
                    report_parameters = args.report_parameters,
                    schema = args.rejection_tool,
                    include_header = False,
                    stat_patterns = stat_patterns,
                    staging_dir = staging_dir)
            WORK_FORCE.put(worker)
            msbayes_workers.append(worker)
    model_indices = models_to_configs.keys()
    info.write('\t[[prior_configs]]\n')
    for model_idx, cfg in models_to_configs.iteritems():
        info.write('\t\t{0} = {1}\n'.format(model_idx, cfg))

    # put observed-data-generating msbayes workers in the queue
    observed_model_idx = configs_to_models.get(args.observed_config, 0)
    schema = args.rejection_tool
    for i in range(num_observed_workers):
        worker = MsBayesWorker(
                temp_fs = working_observed_temp_fs,
                sample_size = prior_subsample_size,
                config_path = args.observed_config,
                model_index = observed_model_idx,
                report_parameters = args.report_parameters,
                schema = schema,
                include_header = True,
                stat_patterns = stat_patterns,
                staging_dir = staging_dir)
        WORK_FORCE.put(worker)
        msbayes_workers.append(worker)
    if remainder_observed_reps > 0:
        worker = MsBayesWorker(
                temp_fs = working_observed_temp_fs,
                sample_size = remainder_observed_reps,
                config_path = args.observed_config,
                model_index = observed_model_idx,
                report_parameters = args.report_parameters,
                schema = schema,
                include_header = True,
                stat_patterns = stat_patterns,
                staging_dir = staging_dir)
        WORK_FORCE.put(worker)
        msbayes_workers.append(worker)
    info.write('\t[[observed_configs]]\n')
    info.write('\t\t{0} = {1}\n'.format(observed_model_idx,
        args.observed_config))

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

    # sort out the finished workers
    msbayes_observed_workers = []
    msbayes_prior_workers = {}
    observed_header = None
    prior_header = None
    for w in msbayes_workers:
        if w.include_header:
            msbayes_observed_workers.append(w)
            if not observed_header:
                observed_header = w.header
        else:
            if w.model_index in msbayes_prior_workers.keys():
                msbayes_prior_workers[w.model_index].append(w)
            else:
                msbayes_prior_workers[w.model_index] = [w]
            if not prior_header:
                prior_header = w.header
    assert len(msbayes_prior_workers) == len(args.prior_configs)
    assert prior_header and prior_header == observed_header

    # merged simulated observed data into one file
    observed_path = os.path.join(base_dir, 'observed.txt')
    info.write('\t[[observed_paths]]\n')
    info.write('\t\t{0} = {1}\n'.format(observed_model_idx, observed_path))
    if args.reps > 0:
        obs_header_path = working_observed_temp_fs.get_file_path(
                prefix = 'observed-header-',
                create = False)
        merge_priors(workers = msbayes_observed_workers,
                prior_path = observed_path,
                header_path = obs_header_path,
                include_header = True)
        lc = line_count(observed_path, ignore_headers=True)
        if lc != args.reps:
            raise Exception('The number of observed simulations ({0}) does '
                    'not match the number of reps ({1})'.format(lc,
                            args.reps))
    if not args.keep_temps:
        _LOG.debug('purging working observed...')
        working_observed_temp_fs.purge()

    # merge model priors into one file per model
    prior_temp_fs = TempFileSystem(parent = base_temp_dir,
            prefix = 'pymsbayes-prior-files-')
    prior_paths = {}
    prior_paths['header'] = prior_temp_fs.get_file_path(
            prefix = 'prior-header-',
            create = False)
    include_header = False
    if args.rejection_tool == 'abctoolbox':
        include_header = True
    for mod_idx in model_indices:
        prior_path = prior_temp_fs.get_file_path(
                prefix = 'prior-{0}-{1}-'.format(mod_idx,
                        args.num_prior_samples),
                create = False)
        merge_priors(workers = msbayes_prior_workers[mod_idx],
                prior_path = prior_path,
                header_path = prior_paths['header'],
                include_header = include_header)
        lc = line_count(prior_path, ignore_headers=True)
        if lc != args.num_prior_samples:
            raise Exception('The number of prior samples ({0}) for model '
                    '{1} does not match num_prior_samples ({2})'.format(
                            lc, model_idx, args.num_prior_samples))
        prior_paths[mod_idx] = prior_path

    if not args.keep_temps:
        _LOG.debug('purging working observed...')
        working_prior_temp_fs.purge()

    # merge all model priors into one file if requested
    if args.merge_priors:
        ntotal = args.num_prior_samples * len(args.prior_configs)
        merged_path = prior_temp_fs.get_file_path(
                prefix = 'prior-merged-{0}-'.format(ntotal),
                create = False)
        merge_prior_files(
                paths = [prior_paths[i] for i in model_indices],
                dest_path = merged_path)
        lc = line_count(merged_path, ignore_headers=True)
        if lc != ntotal:
            raise Exception('The number of prior samples ({0}) in the '
                    'merged prior file does not match the expected '
                    'number of samples ({1})'.format(lc, ntotal))
        prior_paths['merged'] = merged_path

    info.write('\t[[prior_paths]]\n')
    for idx, path in prior_paths.iteritems():
        info.write('\t\t{0} = {1}\n'.format(idx, path))

    ##########################################################################
    ## begin rejection and regression

    if staging_dir:
        rejection_temp_fs = TempFileSystem(parent = staging_dir,
                prefix = 'pymsbayes-rejection-temp-')
    else:
        rejection_temp_fs = TempFileSystem(parent = base_temp_dir,
                prefix = 'pymsbayes-rejection-temp-')
    stat_indices = get_stat_indices(prior_header,
            stat_patterns)
    continuous_parameter_indices = get_parameter_indices(prior_header,
            continuous_patterns)
    discrete_parameter_indices = get_parameter_indices(prior_header,
            discrete_patterns)
    info.write('\t[[column_indices]]\n')
    info.write('\t\tstat_indices = {0}\n'.format(', '.join(
            [str(i) for i in stat_indices])))
    info.write('\t\tcontinuous_indices = {0}\n'.format(', '.join(
            [str(i) for i in continuous_parameter_indices])))
    info.write('\t\tdiscrete_indices = {0}\n'.format(', '.join(
            [str(i) for i in discrete_parameter_indices])))
    reject_workers = []
    if args.merge_priors:
        reject_workers.extend(assemble_rejection_workers(
                temp_fs = rejection_temp_fs,
                observed_sims_file = observed_path,
                prior_path = prior_paths['merged'],
                num_posterior_samples = args.num_posterior_samples,
                num_prior_samples = ntotal,
                results_dir = base_dir,
                posterior_prefix = 'unadjusted-posterior',
                stat_indices = stat_indices,
                continuous_parameter_indices = continuous_parameter_indices,
                discrete_parameter_indices = discrete_parameter_indices,
                regress = args.regression,
                rejection_tool = args.rejection_tool,
                regression_method = args.regression_method,
                keep_temps = args.keep_temps,
                bandwidth = args.bandwidth,
                num_posterior_quantiles = args.num_posterior_quantiles))
    else:
        for i in model_indices:
            reject_workers.extend(assemble_rejection_workers(
                    temp_fs = rejection_temp_fs,
                    observed_sims_file = observed_path,
                    prior_path = prior_paths[i],
                    num_posterior_samples = args.num_posterior_samples,
                    num_prior_samples = args.num_prior_samples,
                    results_dir = base_dir,
                    posterior_prefix = 'unadjusted-posterior',
                    stat_indices = stat_indices,
                    continuous_parameter_indices = continuous_parameter_indices,
                    discrete_parameter_indices = discrete_parameter_indices,
                    regress = args.regression,
                    rejection_tool = args.rejection_tool,
                    regression_method = args.regression_method,
                    keep_temps = args.keep_temps,
                    bandwidth = args.bandwidth,
                    num_posterior_quantiles = args.num_posterior_quantiles))

    # run parallel rejection processes
    for w in reject_workers:
        WORK_FORCE.put(w)
    reject_result_queue = multiprocessing.Queue()
    reject_managers = []
    for i in range(args.np):
        m = Manager(work_queue = WORK_FORCE,
                result_queue = reject_result_queue)
        m.start()
        reject_managers.append(m)
    for i in range(len(reject_workers)):
        reject_workers[i] = reject_result_queue.get()
    for m in reject_managers:
        m.join()
    assert WORK_FORCE.empty()
    assert reject_result_queue.empty()

    if not args.keep_temps:
        rejection_temp_fs.purge()
    if not args.keep_priors:
        prior_temp_fs.purge()

    stop_time = datetime.datetime.now()
    info.write('\t[[run_stats]]\n')
    info.write('\t\tstart_time = {0}\n'.format(str(start_time)))
    info.write('\t\tstop_time = {0}\n'.format(str(stop_time)))
    info.write('\t\tduration = {0}\n'.format(str(stop_time - start_time)))


if __name__ == '__main__':
    main_cli()

