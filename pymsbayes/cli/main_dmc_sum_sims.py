#! /usr/bin/env python

"""
CLI program for summarizing the results of PyMsBayes simulation-based analyses.
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

from pymsbayes.utils import argparse_utils
from pymsbayes.utils.parsing import spreadsheet_iter
from pymsbayes.utils.sumresults import DMCSimulationResults
from pymsbayes.config import MsBayesConfig
from pymsbayes.plotting import (PowerPlotGrid, ProbabilityPowerPlotGrid,
        AccuracyPowerPlotGrid)

_program_info = {
    'name': os.path.basename(__file__),
    'author': 'Jamie Oaks',
    'version': 'Version 0.1',
    'description': __doc__,
    'copyright': 'Copyright (C) 2015 Jamie Oaks',
    'license': 'GNU GPL version 3 or later',}

class PowerResult(object):
    def __init__(self):
        self.true = []
        self.mode = []
        self.mode_glm = []
        self.median = []
        self.prob = []
        self.prob_glm = []

    @classmethod
    def parse_result_summary(cls, result_summary_path):
        psi = cls()
        cv = cls()
        for d in spreadsheet_iter([result_summary_path]):
            psi.true.append(int(d['psi_true']))
            psi.mode.append(int(d['psi_mode']))
            psi.mode_glm.append(int(round(float(d['psi_mode_glm']))))
            psi.prob.append(float(d['psi_1_prob']))
            psi.prob_glm.append(float(d['psi_1_prob_glm']))
            cv.true.append(float(d['cv_true']))
            cv.mode.append(float(d['cv_mode']))
            cv.median.append(float(d['cv_median']))
            cv.mode_glm.append(float(d['cv_mode_glm']))
            cv.prob.append(float(d['cv_prob_less']))
            cv.prob_glm.append(float(d['cv_prob_less_glm']))
        return psi, cv

def parse_results(dmc_sim):
    psi_results = {}
    cv_results = {}
    observed_index_to_name = {}
    observed_configs = {}
    for k, v in dmc_sim.observed_index_to_config.iteritems():
        observed_index_to_name[k] = os.path.splitext(os.path.basename(v))[0]
        observed_configs[k] = MsBayesConfig(v)

    prior_index_to_name = {}
    prior_index_to_model_type = {}
    for k, v in dmc_sim.prior_index_to_config.iteritems():
        n = os.path.splitext(os.path.basename(v))[0]
        prior_index_to_name[k] = n
        c = MsBayesConfig(v)
        if c.implementation == 'old':
            prior_index_to_model_type[k] = 'psi'
        else:
            prior_index_to_model_type[k] = 'dpp'
        
    for observed_index, cfg in observed_configs.iteritems():
        observed_name = observed_index_to_name[observed_index]
        if not psi_results.has_key(observed_name):
            psi_results[observed_name] = {}
        if not cv_results.has_key(observed_name):
            cv_results[observed_name] = {}
        for prior_index in prior_index_to_name.iterkeys():
            if not psi_results[observed_name].has_key(prior_index):
                psi_results[observed_name][prior_index] = {}
            if not cv_results[observed_name].has_key(prior_index):
                cv_results[observed_name][prior_index] = {}
            result_path = dmc_sim.get_result_summary_path(observed_index,
                    prior_index)
            psi, cv = PowerResult.parse_result_summary(result_path)
            psi_results[observed_name][prior_index][cfg] = psi
            cv_results[observed_name][prior_index][cfg] = cv
    return (psi_results, cv_results, prior_index_to_name,
            prior_index_to_model_type)

def create_plots(info_path):
    dmc_sim = DMCSimulationResults(info_path)
    prior_configs = {}
    for k, v in dmc_sim.prior_index_to_config.iteritems():
        prior_configs[k] = MsBayesConfig(v)
    output_dir = os.path.join(os.path.dirname(info_path), 'plots')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    (psi_res, cv_res, prior_index_to_name,
            prior_index_to_model_type) = parse_results(dmc_sim)
    for prior_index, prior_name in prior_index_to_name.iteritems():
        cfg_to_psi = {}
        cfg_to_psi_prob = {}
        cfg_to_cv = {}
        cfg_to_cv_prob = {}
        cfg_to_cv_true_ests = {}
        for observed_name in psi_res.iterkeys():
            div_model_prior = prior_index_to_model_type[prior_index]
            dpp_concentration_mean = None
            if div_model_prior == 'dpp':
                dpp_concentration_mean = prior_configs[
                        prior_index].dpp_concentration.mean
            psi_results = psi_res[observed_name][prior_index]
            cv_results = cv_res[observed_name][prior_index]
            prefix = '_'.join([observed_name, prior_name])
            for cfg, psi in psi_results.iteritems():
                cfg_to_psi[cfg] = psi.mode
                cfg_to_psi_prob[cfg] = psi.prob
            for cfg, cv in cv_results.iteritems():
                cfg_to_cv[cfg] = cv.median
                cfg_to_cv_prob[cfg] = cv.prob
                cfg_to_cv_true_ests[cfg] = {'x': cv.true, 'y': cv.median}

        psi_plot = PowerPlotGrid(
                observed_config_to_estimates = cfg_to_psi,
                variable = 'psi',
                variable_symbol = r'|\mathbf{\tau}|',
                num_columns = 2,
                margin_top = 0.975)
        fig = psi_plot.create_grid()
        fig.savefig(os.path.join(output_dir,
                prefix + '_power_psi_mode.pdf'))

        psi_prob_plot = ProbabilityPowerPlotGrid(
                observed_config_to_estimates = cfg_to_psi_prob,
                variable = 'psi',
                variable_symbol = r'|\mathbf{\tau}|',
                div_model_prior = div_model_prior,
                dpp_concentration_mean = dpp_concentration_mean,
                bayes_factor = 10,
                num_columns = 2)
        fig = psi_prob_plot.create_grid()
        fig.savefig(os.path.join(output_dir,
                prefix + '_power_psi_prob.pdf'))

        cv_accuracy_plot = AccuracyPowerPlotGrid(
                observed_config_to_estimates = cfg_to_cv_true_ests,
                variable_symbol = r'CV_T',
                num_columns = 2,
                padding_between_vertical = 2.0,
                margin_left = 0.04,
                margin_bottom = 0.03)
        fig = cv_accuracy_plot.create_grid()
        fig.savefig(os.path.join(output_dir,
                prefix + '_power_accuracy_cv_median.pdf'))

def main_cli(argv = sys.argv):
    description = '{name} {version}'.format(**_program_info)
    parser = argparse.ArgumentParser(description = description,
            formatter_class = argparse_utils.SmartHelpFormatter)

    parser.add_argument('info_path',
            metavar='PYMSBAYES-INFO-PATH',
            type=argparse_utils.arg_is_file,
            help=('Path to the "pymsbayes-info.txt" file.'))
    parser.add_argument('--plot',
            action = 'store_true',
            help = 'Create plots from result summaries.')
    parser.add_argument('--quiet',
            action = 'store_true',
            help = 'Run without verbose messaging.')
    parser.add_argument('--debug',
            action = 'store_true',
            help = 'Run in debugging mode.')

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    ##########################################################################
    ## handle args

    from pymsbayes.utils.messaging import LoggingControl

    LoggingControl.set_logging_level("INFO")
    if args.quiet:
        LoggingControl.set_logging_level("WARNING")
    if args.debug:
        LoggingControl.set_logging_level("DEBUG")
    log = LoggingControl.get_logger(__name__)

    from pymsbayes.utils import sumresults

    results = sumresults.DMCSimulationResults(args.info_path)
    prior_indices = results.prior_index_to_config.keys()
    test_path = results.get_result_summary_path(
            results.observed_index_to_path.keys()[0],
            prior_indices[0])
    if os.path.exists(test_path):
        log.warning('summary files already exists; skipping summaries!')
    else:
        results.write_result_summaries(
                prior_indices = prior_indices,
                include_tau_exclusion_info = False)
    if args.plot:
        create_plots(args.info_path)



