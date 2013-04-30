#! /usr/bin/env python

import os
import sys

from pymsbayes.workers import (MsBayesWorker, ABCToolBoxRejectWorker,
        ABCToolBoxRegressWorker, RegressionWorker, get_stat_indices,
        parse_header, DEFAULT_STAT_PATTERNS, PARAMETER_PATTERNS,
        merge_prior_files)
from pymsbayes.manager import Manager
from pymsbayes.fileio import process_file_arg
from pymsbayes.utils import GLOBAL_RNG, WORK_FORCE
from pymsbayes.utils.functions import (long_division, least_common_multiple,
        get_random_int, list_splitter)
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class ABCTeam(object):
    count = 0
    def __init__(self,
            temp_fs,
            observed_sims_file,
            model_indices_to_config_paths,
            num_prior_samples,
            num_processors,
            num_posterior_samples = 1000,
            batch_size = 10000,
            rng = None,
            sort_index = None,
            report_parameters = True,
            stat_patterns = DEFAULT_STAT_PATTERNS,
            parameter_patterns = PARAMETER_PATTERNS,
            msbayes_exe_path = None,
            abctoolbox_exe_path = None,
            keep_temps = False,
            work_queue = WORK_FORCE):
        if not rng:
            rng = GLOBAL_RNG
        self.rng = rng
        self.__class__.count += 1
        self.name = 'RejectionTeam-' + str(self.count)
        self.temp_fs = temp_fs
        self.observed_temp_dir = temp_fs.create_subdir(
                prefix = 'observed-files-')
        self.observed_sims_path = observed_sims_path
        self.models = model_indices_to_config_paths
        self.num_prior_samples = num_prior_samples
        self.num_processors = num_processors
        self.num_posterior_samples = num_posterior_samples
        self.batch_size = batch_size
        self.total_samples = self.num_prior_samples * len(self.models)
        max_samples_per_process = self.total_samples / self.num_processors
        if max_samples_per_process < self.batch_size:
            self.batch_size = max_samples_per_process
        if self.batch_size < self.num_posterior_samples:
            self.batch_size = self.num_posterior_samples
        self.num_batches, self.num_extra_samples = long_division(
                self.num_prior_samples, 
                self.batch_size)
        self.sort_index = sort_index
        self.report_parameters = report_parameters
        self.msbayes_exe_path = msbayes_exe_path
        self.abctoolbox_exe_path = abctoolbox_exe_path
        self.rejection_teams = []
        self.prior_workers = []
        self.keep_temps = keep_temps
        self.work_queue = work_queue
        self.result_queue = multiprocessing.Queue()
        self._assemble_rejection_teams()
        self._assemble_prior_workers()

    def _assemble_rejection_teams(self):
        obs_file, close = process_file_arg(self.observed_sims_path)
        header = parse_header(obs_file)
        all_stat_indices = get_stat_indices(header,
                stat_patterns=ALL_STAT_PATTERNS)
        obs_file.next() # header line
        for i, line in enumerate(obs_file):
            obs_path = temp_fs.get_file_path(parent = self.observed_temp_dir,
                    prefix = 'observed-{0}-'.format(i+1),
                    create = False)
            l = line.strip().split()
            with open(obs_path, 'w') as out:
                out.write('{0}\n{1}\n'.format(
                        '\t'.join([header[idx] for idx in all_stat_indices]),
                        '\t'.join([str(l[idx]) for idx in all_stat_indices])))
            self.rejection_teams.append(RejectionTeam(
                    temp_fs = self.temp_fs,
                    prior_workers = [],
                    observed_path = obs_path,
                    num_posterior_samples = self.num_posterior_samples,
                    abctoolbox_exe_path = self.abctoolbox_exe_path,
                    keep_temps = self.keep_temps))
        if close:
            obs_file.close()
    
    def _assemble_prior_workers(self):
        prior_workers = []
        for model_index, config_path in self.models:
            for i in range(self.num_batches):
                seed = get_random_int(self.rng)
                worker = MsBayesWorker(
                        temp_fs = self.temp_fs,
                        sample_size = self.batch_size,
                        config_path = config_path,
                        model_index = model_index,
                        report_parameters = self.report_parameters,
                        seed = seed,
                        schema = 'abctoolbox',
                        stat_patterns = self.stat_patterns)
                prior_workers.append(worker)
            if self.extra_samples > 0:
                worker = MsBayesWorker(
                        temp_fs = self.temp_fs,
                        sample_size = self.extra_samples,
                        config_path = config_path,
                        model_index = model_index,
                        report_parameters = self.report_parameters,
                        seed = seed,
                        schema = 'abctoolbox',
                        stat_patterns = self.stat_patterns)
                prior_workers.append(worker)
        self.prior_workers = list(list_splitter(prior_workers,
                self.num_processors))

    def run(self):
        for prior_worker_batch in self.prior_workers:
            for pw in prior_worker_batch:
                self.work_queue.put(pw)
            managers = []
            for i in range(self.num_processors):
                m = Manager(work_queue = self.work_queue,
                        result_queue = self.result_queue)
                m.start()
                managers.append(m)
            for i in range(len(prior_worker_batch)):
                prior_worker_batch[i] = self.result_queue.get()
            for m in managers:
                m.join()
            assert self.work_queue.empty()
            assert self.result_queue.empty()
            for rt in self.rejection_teams:
                rt.prior_workers = prior_worker_batch
                self.work_queue.put(rt)
            managers = []
            for i in range(self.num_processors):
                m = Manager(work_queue = self.work_queue,
                        result_queue = self.result_queue)
                m.start()
                managers.append(m)
            for i in range(len(self.rejection_teams)):
                self.rejection_teams[i] = self.result_queue.get()
            for m in managers:
                m.join()
            assert self.work_queue.empty()
            assert self.result_queue.empty()


class RejectionTeam(object):
    count = 0
    def __init__(self,
            temp_fs,
            prior_workers,
            observed_path,
            num_posterior_samples = 1000,
            abctoolbox_exe_path = None,
            keep_temps = False):
        self.__class__.count += 1
        self.name = 'RejectionTeam-' + str(self.count)
        self.temp_fs = temp_fs
        self.observed_path = observed_path
        self.num_posterior_samples = num_posterior_samples
        self.abctoolbox_exe_path = abctoolbox_exe_path
        self.prior_workers = prior_workers
        self.keep_temps = keep_temps
        self.posterior_path = None
        self.num_calls = 0

    def start(self):
        self.num_calls += 1
        for i in range(len(self.prior_workers)):
            # pw = self.prior_workers.pop(0)
            pw = self.prior_workers[i]
            merge_paths = [pw.prior_path]
            if self.posterior_path:
                merge_paths.append(self.posterior_path)
            new_prior_path = self.temp_fs.get_file_path(
                    prefix = 'prior-{0}-{1}-{2}-'.format(self.name, 
                            self.num_calls, i+1),
                    create = False)
            merge_prior_files(
                    paths = merge_paths,
                    dest_path = new_prior_path)
            new_posterior_path = self.temp_fs.get_file_path(
                    prefix = 'posterior-{0}-{1}-{2}-'.format(self.name, 
                            self.num_calls, i+1),
                    create = False)
            if self.posterior_path and not self.keep_temps:
                os.remove(self.posterior_path)
            self.posterior_path = new_posterior_path
            rw = ABCToolBoxRejectWorker(
                    temp_fs = self.temp_fs,
                    observed_path = self.observed_path,
                    prior_path = new_prior_path,
                    num_posterior_samples = self.num_posterior_samples,
                    posterior_path = new_posterior_path,
                    regression_worker = None,
                    exe_path = self.abctoolbox_exe_path,
                    keep_temps = self.keep_temps,
                    max_read_sims = pw.sample_size + 10)
            rw.start()
            if not self.keep_temps:
                os.remove(rw.prior_path)

