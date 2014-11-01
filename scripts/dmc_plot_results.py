#! /usr/bin/env python

"""
CLI for generating plots of empirical results from dmc.py.
"""

import os
import sys
import random
import multiprocessing
import argparse

from pymsbayes.utils import argparse_utils

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
    parser.add_argument('info_path',
            metavar='PYMSBAYES-INFO-FILE',
            type=argparse_utils.arg_is_file,
            help=('Path to `pymsbayes-info.txt` file.'))
    parser.add_argument('-n', '--num-prior-samples',
            action = 'store',
            type = argparse_utils.arg_is_positive_int,
            default = 100000,
            help = ('The number of prior samples to simulate for estimating'
                    'prior probabilities.'))
    parser.add_argument('-i', '--iteration-index',
            action = 'store',
            type = argparse_utils.arg_is_positive_int,
            help = ('The batch iteration index of results to be summarized. '
                    'Output files should have a consistent schema. For '
                    'example, a results file for divergence models might look '
                    'something like '
                    '`d1-m1-s1-151-div-model-results.txt`. In this example, '
                    'the batch iteration index is "151". The default is to '
                    'use the largest batch iteration index, which is probably '
                    'what you want.'))
    parser.add_argument('-o', '--output-dir',
            action = 'store',
            type = argparse_utils.arg_is_dir,
            help = ('The directory in which all output plots will be written. '
                    'The default is to use the directory of the pymsbayes info '
                    'file.'))
    parser.add_argument('--np',
            action = 'store',
            type = argparse_utils.arg_is_positive_int,
            default = multiprocessing.cpu_count(),
            help = ('The maximum number of processes to run in parallel. The '
                    'default is the number of CPUs available on the machine.'))
    parser.add_argument('-m', '--mu',
            action = 'store',
            type = argparse_utils.arg_is_positive_float,
            default = None,
            help = ('The mutation rate with which to scale time to units of '
                    'generations. The default is to keep the timescale in '
                    'units of Nc generations.'))
    parser.add_argument('--seed',
            action = 'store',
            type = argparse_utils.arg_is_positive_int,
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

    from pymsbayes import plotting
    from pymsbayes.utils import sumresults
    from pymsbayes.utils import GLOBAL_RNG

    if not plotting.MATPLOTLIB_AVAILABLE:
        log.error(
                '`matplotlib` could not be imported, so plots can not be\n'
                'produced. Please install `matplotlib` and try again.')
        sys.exit(1)

    if not args.seed:
        args.seed = random.randint(1, 999999999)
    GLOBAL_RNG.seed(args.seed)

    if not args.output_dir:
        args.output_dir = os.path.dirname(args.info_path)
    args.output_dir = os.path.join(args.output_dir, 'plots')
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    results = sumresults.DMCSimulationResults(args.info_path)
    if results.num_sim_reps > 1:
        log.error('Results appear to be from simulation-based analysis, '
                'for which this plotting script is not appropriate.')
        sys.exit(1)

    observed_indices = sorted(results.observed_index_to_config.keys())
    prior_indices = sorted(results.prior_index_to_config.keys())
    for obs_idx in observed_indices:
        for prior_idx in prior_indices:
            result_indices = results.get_result_indices(obs_idx, prior_idx, 1)
            result_idx = max(result_indices)
            result_path_prefix = '{0}{1}-'.format(
                    results.get_result_path_prefix(obs_idx, prior_idx, 1),
                    result_idx)
            result_dir = os.path.dirname(result_path_prefix)
            out_prefix = os.path.join(args.output_dir, os.path.basename(
                    result_path_prefix))
            prior_cfg = results.prior_configs[prior_idx]
            posterior_summary_path = get_result_path(result_path_prefix,
                    'posterior-summary')
            div_model_path = get_result_path(result_path_prefix,
                    'div-model-results')
            config_path = results.prior_index_to_config[prior_idx]
            time_multiplier = 1.0
            if args.mu is not None:
                try:
                    mean_theta = prior_cfg.theta.mean
                except:
                    mean_theta = prior_cfg.d_theta.mean
                time_multiplier = mean_theta / args.mu

            if results.sort_index == 0:
                #plot marginal times
                if not posterior_summary_path:
                    log.warning('Could not find {0}{1}.txt(.gz); '
                            'Skipping marginal times plot...'.format(
                                    result_path_prefix,
                                    'posterior-summary'))
                else:
                    label_dimension = (0.34 * (prior_cfg.npairs + 1)) + 0.56
                    marginal_times_plot = plotting.get_marginal_divergence_time_plot(
                            config_path = config_path,
                            posterior_summary_path = posterior_summary_path,
                            labels = None,
                            estimate = 'median',
                            interval = 'HPD_95_interval',
                            time_multiplier = time_multiplier,
                            horizontal = True,
                            label_dimension = label_dimension,
                            measure_dimension = 8.0,
                            label_size = 12.0,
                            measure_tick_label_size = 12.0,
                            measure_axis_label = 'Divergence time',
                            measure_axis_label_size = 14.0,
                            label_axis_label = 'Taxon pair',
                            label_axis_label_size = 14.0,
                            usetex = False)
                    marginal_times_path = '{0}{1}'.format(out_prefix,
                            'marginal-divergence-times.pdf')
                    marginal_times_plot.savefig(marginal_times_path)

                #plot top ordered models
                if not div_model_path:
                    log.warning('Could not find {0}{1}.txt(.gz); '
                            'Skipping ordered div model plot...'.format(
                                    result_path_prefix,
                                    'div-model-results'))
                else:
                    width = (0.38 * prior_cfg.npairs) + 1.5
                    div_model_plot = plotting.OrderedDivergenceModelPlotGrid(
                            div_model_results_path = div_model_path,
                            config_path = config_path,
                            num_top_models = 10,
                            time_multiplier = time_multiplier,
                            height = 12.0,
                            width = width,
                            plot_label_schema = 'uppercase',
                            plot_label_offset = 0,
                            plot_label_size = 12.0,
                            y_title = 'Divergence time',
                            y_title_size = 14.0,
                            y_tick_label_size = 10.0,
                            right_text_size = 10.0,
                            margin_left = 0.03,
                            margin_bottom = 0.0,
                            margin_right = 1,
                            margin_top = 0.99,
                            padding_between_vertical = 0.8)
                    plot = div_model_plot.create_grid()
                    div_model_plot_path = '{0}{1}'.format(out_prefix,
                            'ordered-div-models.pdf')
                    plot.savefig(div_model_plot_path)

            else:
                #plot top unordered models
                if not div_model_path:
                    log.warning('Could not find {0}{1}.txt(.gz); '
                            'Skipping unordered div model plot...'.format(
                                    result_path_prefix,
                                    'div-model-results'))
                else:
                    width = (0.38 * prior_cfg.npairs) + 1.5
                    div_model_plot = plotting.UnorderedDivergenceModelPlotGrid(
                            div_model_results_path = div_model_path,
                            num_top_models = 10,
                            time_multiplier = time_multiplier,
                            height = 10.0,
                            width = width,
                            data_label_size = 10.0,
                            plot_label_schema = 'uppercase',
                            plot_label_offset = 0,
                            plot_label_size = 12.0,
                            y_title = 'Divergence time',
                            y_title_size = 14.0,
                            y_tick_label_size = 10.0,
                            right_text_size = 10.0,
                            margin_left = 0.03,
                            margin_bottom = 0.0,
                            margin_right = 1,
                            margin_top = 0.99,
                            padding_between_vertical = 0.8)
                    plot = div_model_plot.create_grid()
                    div_model_plot_path = '{0}{1}'.format(out_prefix,
                            'ordered-div-models.pdf')
                    plot.savefig(div_model_plot_path)

            #plot ndiv plot
            psi_path = get_result_path(result_path_prefix,
                    'psi-results')
            if not psi_path:
                log.warning('Could not find {0}{1}.txt(.gz); '
                        'Skipping number of divergences plot...'.format(
                                result_path_prefix,
                                'psi-results'))
            else:
                width = (0.25 * prior_cfg.npairs) + 0.55
                if width < 2.8:
                    width = 2.8
                num_div_summary = plotting.NumberOfDivergencesSummary(
                        config_path = results.prior_index_to_config[prior_idx],
                        psi_results_path = psi_path,
                        posterior_summary_path = posterior_summary_path,
                        num_prior_samples = args.num_prior_samples,
                        num_processors = args.np)
                num_div_summary.create_plot(
                        plot_label_size = 10.0,
                        right_text_size = 10.0,
                        x_label_size = 10.0,
                        y_label_size = 10.0,
                        xtick_label_size = 10.0,
                        ytick_label_size = 8.0,
                        height = 6.0,
                        width = width,
                        margin_bottom = 0.0,
                        margin_left = 0.0,
                        margin_top = 0.97,
                        margin_right = 1.0,
                        padding_between_vertical = 1.0)
                num_div_plot_path = '{0}{1}'.format(out_prefix,
                        'number-of-divergences.pdf')
                num_div_summary.save_plot(num_div_plot_path)

                bf_plot_path = '{0}{1}'.format(out_prefix,
                        'number-of-divergences-bayes-factors-only.pdf')
                num_div_summary.save_bf_plot(bf_plot_path)
                
                num_div_bf_path = '{0}{1}'.format(out_prefix,
                        'number-of-divergences-bayes-factors.txt')
                with open(num_div_bf_path, 'w') as out:
                    out.write('num_of_divs\t2ln(bf)\n')
                    for n in sorted(num_div_summary.psi_bayes_factors.keys()):
                        out.write('{0}\t{1}\n'.format(n,
                                num_div_summary.psi_bayes_factors[n]))

def get_result_path(path_prefix, name):
    p = '{0}{1}.txt'.format(path_prefix, name)
    if os.path.exists(p):
        return p
    if os.path.exists(p + '.gz'):
        return p + '.gz'
    return None


if __name__ == '__main__':
    main_cli()

