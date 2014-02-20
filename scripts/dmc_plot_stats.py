#! /usr/bin/env python

"""
CLI for generating histograms of summary statistics sampled from a prior.
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

    from pymsbayes.utils.messaging import (get_logger, LOGGING_LEVEL_ENV_VAR,
            InfoLogger)

    os.environ[LOGGING_LEVEL_ENV_VAR] = "INFO"
    if args.quiet:
        os.environ[LOGGING_LEVEL_ENV_VAR] = "WARNING"
    if args.debug:
        os.environ[LOGGING_LEVEL_ENV_VAR] = "DEBUG"
    log = get_logger(__name__)

    from pymsbayes.workers import MsBayesWorker
    from pymsbayes.utils.parsing import (get_patterns_from_prefixes,
            DEFAULT_STAT_PATTERNS, get_dict_from_spreadsheets, dict_line_iter)
    from pymsbayes.manager import Manager
    from pymsbayes.utils.tempfs import TempFileSystem
    from pymsbayes.utils import probability, stats
    from pymsbayes.utils.functions import long_division
    from pymsbayes.config import MsBayesConfig
    from pymsbayes.utils import GLOBAL_RNG, MSBAYES_SORT_INDEX
    from pymsbayes.fileio import process_file_arg
    from pymsbayes import plotting

    MSBAYES_SORT_INDEX.set_index(0)

    if not args.output_dir:
        args.output_dir = os.path.dirname(args.config)
    info = InfoLogger(os.path.join(args.output_dir, 'pymsbayes-info.txt'))

    sample_path = os.path.join(args.output_dir, 'prior-sample.txt')
    if args.compress:
        sample_path += '.gz'

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
    info.write('\tsample_path = {0!r}'.format(sample_path))

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
                config_path = args.config,
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
    sample = get_dict_from_spreadsheets([w.prior_path for w in workers])

    log.info('Writing prior samples...')
    out, close = process_file_arg(sample_path, 'w',
            compresslevel = compress_level)
    for row in dict_line_iter(sample, sep = '\t'):
        out.write(row)
    if close:
        out.close()

    log.info('Creating plots...')

    if not plotting.MATPLOTLIB_AVAILABLE:
        _LOG.warning(
                '`matplotlib` could not be imported, so the plot can not be\n'
                'produced. The data to create the plot can be found in:\n\t'
                '{0!r}'.format(sample_path))
        sys.exit(1)

    for stat_pattern in stat_patterns:
        found = False
        for stat, values in sample.iteritems():
            if stat_pattern.match(stat):
                values = [float(v) for v in values]
                found = True
                plot_path = os.path.join(args.output_dir,
                        'plot-{0}.pdf'.format(stat))
                summary = stats.get_summary(values)
                s = r'mean = {0:.4f} ({1:.4f}-{2:.4f})'.format(
                        summary['mean'],
                        summary['qi_95'][0],
                        summary['qi_95'][1])
                hd = plotting.HistData(x = values,
                        normed = True,
                        bins = 20,
                        histtype = 'bar',
                        align = 'mid',
                        orientation = 'vertical',
                        zorder = 0)
                hist = plotting.ScatterPlot(hist_data_list = [hd],
                        right_text = s)
                hist.left_text_size = 12.0
                hist.right_text_size = 12.0
                xticks = [i for i in hist.ax.get_xticks()]
                xtick_labels = [i for i in xticks]
                yticks = [i for i in hist.ax.get_yticks()]
                ytick_labels = [i for i in yticks]
                if len(xtick_labels) >= 8:
                    for i in range(1, len(xtick_labels), 2):
                        xtick_labels[i] = ''
                if len(ytick_labels) >= 8:
                    for i in range(1, len(ytick_labels), 2):
                        ytick_labels[i] = ''
                xticks_obj = plotting.Ticks(ticks = xticks,
                        labels = xtick_labels,
                        horizontalalignment = 'center')
                yticks_obj = plotting.Ticks(ticks = yticks,
                        labels = ytick_labels)
                hist.xticks_obj = xticks_obj
                hist.yticks_obj = yticks_obj

                plot_grid = plotting.PlotGrid(subplots = [hist],
                        num_columns = 1,
                        label_schema = None,
                        title = stat,
                        title_size = 14.0,
                        title_top = False,
                        y_title = 'Density',
                        y_title_position = 0.001,
                        y_title_size = 14.0,
                        height = 4.0,
                        width = 6.0,
                        auto_height = False)
                plot_grid.auto_adjust_margins = False
                plot_grid.margin_left = 0.04
                plot_grid.margin_bottom = 0.04 
                plot_grid.margin_right = 1.0 
                plot_grid.margin_top = 0.97
                plot_grid.reset_figure()
                plot_grid.savefig(plot_path)

        if not found:
            raise Exception('stat pattern {0!r} not found in simulated stats:'
                    '\n\t{1}'.format(stat_pattern, ', '.join(sample.keys())))

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

