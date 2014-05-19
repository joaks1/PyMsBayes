#! /usr/bin/env python

"""
CLI for generating saturation plots.
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
    parser.add_argument('-o', '--output-dir',
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
            default = ['pi', 'pi.net', 'wattTheta', 'tajD.denom'],
            help = ('Prefixes of summary statistics to use in the analyses. '
                    'The prefixes should be separated by spaces. '
                    'Default: `-s pi pi.net wattTheta tajD.denom`.'))
    parser.add_argument('--vertical-lines',
            nargs = '*',
            type = float,
            default = [],
            help = ('Positions along x-axis where vertical lines are to be '
                    'drawn. Default is to draw no vertical lines.'))
    parser.add_argument('--compress',
            action = 'store_true',
            help = 'Compress plot data file.')
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

    from pymsbayes.workers import MsBayesWorker
    from pymsbayes.utils.parsing import (get_patterns_from_prefixes,
            DEFAULT_STAT_PATTERNS, get_stats_by_time, dict_line_iter)
    from pymsbayes.manager import Manager
    from pymsbayes.utils.tempfs import TempFileSystem
    from pymsbayes.utils import probability
    from pymsbayes.utils.functions import long_division
    from pymsbayes.config import MsBayesConfig
    from pymsbayes.utils import GLOBAL_RNG, MSBAYES_SORT_INDEX
    from pymsbayes.fileio import process_file_arg
    from pymsbayes.plotting import MATPLOTLIB_AVAILABLE, SaturationPlotGrid

    MSBAYES_SORT_INDEX.set_index(0)

    if not args.output_dir:
        args.output_dir = os.path.dirname(args.config)
    info = InfoLogger(os.path.join(args.output_dir, 'pymsbayes-info.txt'))

    stats_by_time_path = os.path.join(args.output_dir, 'stats-by-time.txt')
    if args.compress:
        stats_by_time_path += '.gz'
    plot_path = os.path.join(args.output_dir, 'saturation-plot.pdf')

    if not args.temp_dir:
        args.temp_dir = args.output_dir
    temp_fs = TempFileSystem(parent=args.temp_dir, prefix='temp-files-')
    args.stat_prefixes = [s.rstrip('.') for s in args.stat_prefixes]
    stat_patterns = get_patterns_from_prefixes(
            [s + '.' for s in args.stat_prefixes],
            ignore_case=True)
    if not args.seed:
        args.seed = random.randint(1, 999999999)
    GLOBAL_RNG.seed(args.seed)
    compress_level = None
    if args.compress:
        compress_level = 9

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
    info.write('\tinvocation = {0!r}'.format(' '.join(sys.argv)), log.info)
    info.write('\toutput_directory = {0!r}'.format(args.output_dir), log.info)
    info.write('\ttemp_directory = {0!r}'.format(temp_fs.base_dir), log.info)
    info.write('\tsort_index = {0}'.format(
            MSBAYES_SORT_INDEX.current_value()), log.info)
    info.write('\tstat_patterns = {0!r}'.format(
            ', '.join([p.pattern for p in stat_patterns])), log.info)
    info.write('\tseed = {0}'.format(args.seed), log.info)
    info.write('\tnum_prior_samples = {0}'.format(args.num_prior_samples),
            log.info)
    info.write('\tstats_by_time_path = {0!r}'.format(stats_by_time_path))

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
    stat_keys = stats_by_time.keys()
    stat_keys.remove('PRI.t')
    for prefix in args.stat_prefixes:
        if not prefix in stat_keys:
            raise Exception('stat prefix {0!r} not found in simulated stats:'
                    '\n\t{1}'.format(prefix, ', '.join(stat_keys)))
    header = ['PRI.t'] + args.stat_prefixes
    log.info('Writing stats-by-time matrix...')
    out, close = process_file_arg(stats_by_time_path, 'w',
            compresslevel = compress_level)
    for row in dict_line_iter(stats_by_time, sep = '\t', header = header):
        out.write(row)
    if close:
        out.close()

    log.info('Creating plots...')

    if not MATPLOTLIB_AVAILABLE:
        log.warning(
                '`matplotlib` could not be imported, so the plot can not be\n'
                'produced. The data to create the plot can be found in:\n\t'
                '{0!r}'.format(stats_by_time_path))
    else:
        y_labels = {'pi': r'$\pi$',
                   'pi.net': r'$\pi_{net}$',
                   'wattTheta': r'$\theta_W$',
                   'tajD.denom': r'$SD(\pi - \theta_W)$'}
        spg = SaturationPlotGrid(stats_by_time,
                x_key = 'PRI.t',
                y_keys = args.stat_prefixes,
                y_labels = y_labels,
                num_columns = 2,
                vertical_line_positions = args.vertical_lines)
        fig = spg.create_grid()
        fig.savefig(plot_path)
        # fig = plt.figure()
        # ncols = 2
        # nrows = get_num_rows(len(stat_keys))
        # for i, k in enumerate(stat_keys):
        #     if fig.axes:
        #         ax = fig.add_subplot(nrows, ncols, i + 1, sharex = fig.axes[0])
        #     else:
        #         ax = fig.add_subplot(nrows, ncols, i + 1)
        #     ax.plot(stats_by_time['PRI.t'], stats_by_time[k])
        #     ax.set_ylabel(y_labels[k])
        # plt.setp([a.lines for a in fig.axes],
        #         marker = 'o',
        #         linestyle='',
        #         markerfacecolor = 'none',
        #         markeredgecolor = '0.4',
        #         markeredgewidth = 0.7)
        # fig.suptitle(r'Divergence time $\tau$ in $4N_C$ generations',
        #         verticalalignment = 'bottom',
        #         y = 0.001)
        # fig.tight_layout(pad = 0.25, # out side margin
        #                  h_pad = None, # height padding between subplots
        #                  w_pad = None, # width padding between subplots
        #                  rect = (0,0.05,1,1))
        # fig.savefig(plot_path)

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

def get_num_rows(n, ncols = 2):
    return int(math.ceil(n / float(ncols)))

if __name__ == '__main__':
    main_cli()

