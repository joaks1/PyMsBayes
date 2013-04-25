#! /usr/bin/env python

import os
import sys

from pymsbayes.workers import (MsBayesWorker, ABCToolBoxRejectWorker,
        ABCToolBoxRegressWorker, RegressionWorker)
from pymsbayes.utils import GLOBAL_RNG
from pymsbayes.utils.functions import long_division, least_common_multiple
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

def assemble_abc_teams(
        temp_fs,
        observed_sims_file,
        prior_config_dict, # model index by cfg path dict
        num_prior_samples,
        num_processors,
        num_posterior_samples = 1000,
        rng = None,
        sort_index = None,
        report_parameters = True,
        stat_patterns = DEFAULT_STAT_PATTERNS,
        parameter_patterns = PARAMETER_PATTERNS,
        msbayes_exe_path = None,
        abctoolbox_exe_path = None,
        keep_temps = False):
    if not rng:
        rng = GLOBAL_RNG
    num_teams = least_common_multiple([num_processors, len(prior_config_dict)])
    num_teams_per_config, r = long_division(num_teams, len(prior_config_dict))
    assert r == 0
    num_samples_per_team, extra_samples = long_division(num_prior_samples,
            num_teams_per_config)

def assemble_final_rejection_workers_from_abc_teams(
        abc_teams,
        regress = True,
        regression_method = 'glm'):
    # create rejection workers for final reject/regression
    pass

class ABCTeam(object):
    def __init__(self,
            temp_fs,
            observed_paths,
            prior_config_path,
            num_prior_samples,
            seeds,
            num_posterior_samples = 1000,
            model_index = None,
            sort_index = None,
            report_parameters = True,
            stat_patterns = DEFAULT_STAT_PATTERNS,
            parameter_patterns = PARAMETER_PATTERNS,
            msbayes_exe_path = None,
            abctoolbox_exe_path = None,
            keep_temps = False):
        self.temp_fs = temp_fs
        self.observed_paths = observed_paths
        self.config_path = prior_config_path
        self.num_prior_samples = num_prior_samples
        self.num_posterior_samples = num_posterior_samples
        self.model_index = model_index
        self.sort_index = sort_index
        self.report_parameters = report_parameters
        self.msbayes_exe_path = msbayes_exe_path
        self.abctoolbox_exe_path = abctoolbox_exe_path
        batch_size, remainder = long_division(self.num_prior_samples,
                len(seeds))
        batch_sizes = [batch_sizes] * len(seeds)
        if remainder > 0:
            batch_sizes[-1] += remainder
        self.seeds = dict((seeds[i], batch_sizes[i]) for i in range(len(seeds)))
        self.prior_workers = []
        self.rejection_workers = []
        self.keep_temps = keep_temps

    def _pre_process(self):
        for seed, nsamples in self.seeds:
            self.prior_workers.append(MsBayesWorker(
                    temp_fs = self.temp_fs,
                    self.sample_size = nsamples,
                    config_path = self.config_path,
                    exe_path = self.msbayes_exe_path,
                    model_index = self.model_index,
                    sort_index = self.sort_index,
                    report_parameters = self.report_parameters,
                    seed = seed,
                    schema = 'abctoolbox',
                    stat_patterns = self.stat_patterns,
                    parameter_patterns = self.parameter_patterns,
                    write_stats_file = False))
        for i, obs_path in enumerate(self.observed_paths):
            self.rejection_workers.append(ABCToolBoxRejectWorker(
                    temp_fs = temp_fs,
                    observed_path = obs_path,
                    prior_path = self.prior_workers[0].prior_path,
                    num_posterior_samples = self.num_posterior_samples,
                    posterior_path = None,
                    regression_worker = None,
                    exe_path = self.abctoolbox_exe_path,
                    keep_temps = self.keep_temps,
                    max_read_sims = max(self.seeds.values()) + 10,
                    tag = str(i)))

    def start(self):
        self._pre_process()
        for i in range(len(self.prior_workers)):
            pw = self.prior_workers[i]
            pw.start()
            for j in range(len(self.rejection_workers)):
                rw = self.rejection_workers[j]
                if rw.finished:
                    rw.finished = False
                    new_prior_path = self.temp_fs.get_file_path(
                            prefix = 'prior-{0}-{1}-{2}-'.format(rw.name,
                                    rw.tag, i),
                            create = False)
                    new_posterior_path = self.temp_fs.get_file_path(
                            prefix = 'posterior-{0}-{1}-{2}-'.format(rw.name,
                                    rw.tag, i),
                            create = False)
                    merge_prior_files(
                            paths = [pw.prior_path, rw.posterior_path],
                            dest_path = new_prior_path)
                    if not self.keep_temps:
                        os.remove(rw.prior_path)
                        os.remove(rw.posterior_path)
                    rw.prior_path = new_prior_path
                    rw.posterior_path = new_posterior_path
                    rw.start()
                else:
                    rw.prior_path = pw.prior_path
                    rw.posterior_path = self.temp_fs.get_file_path(
                            prefix = 'posterior-{0}-{1}-{2}-'.format(rw.name,
                                    rw.tag, i),
                            create = False)
                    rw.start()

