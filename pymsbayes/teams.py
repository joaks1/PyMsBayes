#! /usr/bin/env python

import os
import sys
import copy
import multiprocessing

from pymsbayes.workers import (MsBayesWorker, ABCToolBoxRegressWorker,
        EuRejectSummaryMerger, EuRejectWorker, PosteriorWorker)
from pymsbayes.utils.parsing import (get_stat_indices, parse_header,
        DEFAULT_STAT_PATTERNS, ALL_STAT_PATTERNS)
from pymsbayes.manager import Manager
from pymsbayes.fileio import process_file_arg, expand_path
from pymsbayes.utils import GLOBAL_RNG, WORK_FORCE
from pymsbayes.utils.functions import (long_division, least_common_multiple,
        get_random_int, list_splitter, mk_new_dir)
from pymsbayes.utils.stats import SampleSummaryCollection
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class ABCTeam(object):
    count = 0
    def __init__(self,
            temp_fs,
            observed_stats_file,
            num_taxon_pairs,
            model_indices_to_config_paths,
            num_prior_samples,
            num_processors,
            num_standardizing_samples = 10000,
            num_posterior_samples = 1000,
            num_posterior_density_quantiles = 1000,
            batch_size = 10000,
            output_dir = None,
            output_prefix = '',
            rng = None,
            sort_index = None,
            report_parameters = True,
            stat_patterns = DEFAULT_STAT_PATTERNS,
            eureject_exe_path = None,
            abctoolbox_exe_path = None,
            msbayes_exe_path = None,
            abctoolbox_bandwidth = None,
            omega_threshold = 0.01,
            compress = True,
            keep_temps = False,
            global_estimate_only = False,
            work_queue = WORK_FORCE):
        if not rng:
            rng = GLOBAL_RNG
        self.rng = rng
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.temp_fs = temp_fs
        self.temp_output_dir = temp_fs.create_subdir(
                prefix = self.name + '-temp-output-')
        self.observed_temp_dir = temp_fs.create_subdir(
                parent = self.temp_output_dir,
                prefix = 'observed-files-')
        self.summary_temp_dir = temp_fs.create_subdir(
                parent = self.temp_output_dir,
                prefix = 'summary-files-')
        self.observed_stats_path = expand_path(observed_stats_file)
        if not output_dir:
            output_dir = os.path.dirname(self.observed_stats_path)
        self.output_dir = mk_new_dir(os.path.join(expand_path(output_dir),
                'pymsbayes-output'))
        self.output_prefix = str(output_prefix)
        self.global_estimate_only = global_estimate_only
        self.num_taxon_pairs = num_taxon_pairs
        self.models = model_indices_to_config_paths
        self.num_prior_samples = num_prior_samples
        self.num_processors = num_processors
        self.num_standardizing_samples = num_standardizing_samples
        self.num_posterior_samples = num_posterior_samples
        self.num_posterior_density_quantiles = num_posterior_density_quantiles
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
        self.stat_patterns = stat_patterns
        self.eureject_exe_path = eureject_exe_path
        self.abctoolbox_exe_path = abctoolbox_exe_path
        self.msbayes_exe_path = msbayes_exe_path
        self.abctoolbox_bandwidth = abctoolbox_bandwidth
        self.omega_threshold = omega_threshold
        self.compress = compress
        self.model_dirs = dict(zip(self.models.keys(),
                [mk_new_dir(os.path.join(self.output_dir,
                        'm' + str(k))) for k in self.models.keys()]))
        if len(self.models) > 1:
            self.model_dirs.update({'combined': mk_new_dir(os.path.join(
                    self.output_dir,
                    'm' + ''.join([str(i) for i in sorted(self.models.keys())])))})
        else:
            self.global_estimate_only = False
        self.summary_paths = {}
        for k, d in self.model_dirs.iteritems():
            self.summary_paths[k] = os.path.join(d,
                    self.output_prefix + os.path.basename(d) + \
                            '-stat-means-and-std-devs.txt')
        self.rejection_teams = {}
        if not self.global_estimate_only:
            self.rejection_teams = dict(zip(self.models.keys(),
                    [[] for i in range(len(self.models.keys()))]))
        if len(self.models) > 1:
            self.rejection_teams.update({'combined': []})
        self.prior_summary_workers = []
        self.prior_workers = []
        self.num_observed = 0
        self.keep_temps = keep_temps
        self.work_queue = work_queue
        self.result_queue = multiprocessing.Queue()
        self._assemble_rejection_teams()
        self._assemble_prior_workers()
        self.finished = False

    def _assemble_prior_workers(self):
        prior_workers = []
        for model_index, config_path in self.models.iteritems():
            to_summarize = self.num_standardizing_samples
            for i in range(self.num_batches):
                sum_worker = None
                seed = get_random_int(self.rng)
                worker = MsBayesWorker(
                        temp_fs = self.temp_fs,
                        sample_size = self.batch_size,
                        config_path = config_path,
                        exe_path = self.msbayes_exe_path,
                        model_index = model_index,
                        report_parameters = self.report_parameters,
                        seed = seed,
                        schema = 'abctoolbox',
                        stat_patterns = self.stat_patterns,
                        summary_worker = sum_worker,
                        tag = model_index)
                if to_summarize > 0:
                    if (to_summarize - self.batch_size) >= 0:
                        nsum = self.batch_size
                    else:
                        nsum = to_summarize
                    to_summarize -= nsum
                    sum_path = self.temp_fs.get_file_path(
                            parent = self.summary_temp_dir,
                            prefix = 'summary-{0}-{1}-'.format(model_index,
                                    i+1),
                            create = False)
                    post_path = self.temp_fs.get_file_path(
                            parent = self.summary_temp_dir,
                            prefix = 'posterior-{0}-{1}-'.format(model_index,
                                    i+1),
                            create = False)
                    sum_worker = EuRejectWorker(
                            temp_fs = self.temp_fs,
                            observed_path = self.observed_stats_path,
                            prior_paths = [worker.prior_path],
                            posterior_path = post_path,
                            num_posterior_samples = 0,
                            num_standardizing_samples = nsum,
                            summary_out_path = sum_path,
                            exe_path = self.eureject_exe_path,
                            keep_temps = self.keep_temps,
                            tag = model_index)
                    worker.summary_worker = sum_worker
                if sum_worker:
                    self.prior_summary_workers.append(worker)
                else:
                    prior_workers.append(worker)
            if self.num_extra_samples > 0:
                sum_worker = None
                seed = get_random_int(self.rng)
                if to_summarize > 0:
                    if (to_summarize - self.num_extra_samples) >= 0:
                        nsum = self.num_extra_samples
                    else:
                        nsum = to_summarize
                    to_summarize -= nsum
                    sum_path = temp_fs.get_file_path(
                            parent = self.summary_temp_dir,
                            prefix = 'summary-{0}-{1}-'.format(model_index,
                                    i+1),
                            create = False)
                    post_path = self.temp_fs.get_file_path(
                            parent = self.summary_temp_dir,
                            prefix = 'posterior-{0}-{1}-'.format(model_index,
                                    i+1),
                            create = False)
                    sum_worker = EuRejectWorker(
                            temp_fs = self.temp_fs,
                            observed_path = self.observed_stats_path,
                            prior_paths = [worker.prior_path],
                            posterior_path = post_path,
                            num_posterior_samples = 0,
                            num_standardizing_samples = nsum,
                            summary_out_path = sum_path,
                            exe_path = self.eureject_exe_path,
                            keep_temps = self.keep_temps,
                            tag = model_index)
                    self.summary_workers.append(sum_worker)
                worker = MsBayesWorker(
                        temp_fs = self.temp_fs,
                        sample_size = self.num_extra_samples,
                        config_path = config_path,
                        exe_path = self.msbayes_exe_path,
                        model_index = model_index,
                        report_parameters = self.report_parameters,
                        seed = seed,
                        schema = 'abctoolbox',
                        stat_patterns = self.stat_patterns,
                        summary_worker = sum_worker,
                        tag = model_index)
                if sum_worker:
                    self.prior_summary_workers.append(worker)
                else:
                    prior_workers.append(worker)
        self.prior_workers = list(list_splitter(prior_workers,
                self.num_processors, by_size=True))

    def _assemble_rejection_teams(self):
        obs_file, close = process_file_arg(self.observed_stats_path)
        header = parse_header(obs_file)
        all_stat_indices = get_stat_indices(header,
                stat_patterns=ALL_STAT_PATTERNS)
        obs_file.next() # header line
        for i, line in enumerate(obs_file):
            self.num_observed += 1
            index = i + 1
            tmp_obs_path = self.temp_fs.get_file_path(
                    parent = self.observed_temp_dir,
                    prefix = 'observed-{0}-'.format(index),
                    create = False)
            l = line.strip().split()
            with open(tmp_obs_path, 'w') as out:
                out.write('{0}\n{1}\n'.format(
                        '\t'.join([header[idx] for idx in all_stat_indices]),
                        '\t'.join([str(l[idx]) for idx in all_stat_indices])))
            for model_idx, rejection_team_list in self.rejection_teams.iteritems():
                if model_idx == 'combined':
                    m_indices = self.models.keys()
                else:
                    m_indices = [model_idx]
                rejection_team_list.append(RejectionTeam(
                        temp_fs = self.temp_fs,
                        observed_path = tmp_obs_path,
                        output_dir = self.model_dirs[model_idx],
                        output_prefix = self.output_prefix,
                        num_taxon_pairs = self.num_taxon_pairs,
                        prior_paths = [],
                        model_indices = m_indices,
                        summary_in_path = self.summary_paths[model_idx],
                        num_posterior_samples = self.num_posterior_samples,
                        num_posterior_density_quantiles = \
                                self.num_posterior_density_quantiles,
                        run_regression = False,
                        eureject_exe_path = self.eureject_exe_path,
                        abctoolbox_exe_path = self.abctoolbox_exe_path,
                        abctoolbox_bandwidth = self.abctoolbox_bandwidth,
                        omega_threshold = self.omega_threshold,
                        compress = self.compress,
                        keep_temps = self.keep_temps,
                        index = index,
                        tag = model_idx))
        if close:
            obs_file.close()
        if self.num_observed == 1:
            for model_idx, rt_list in self.rejection_teams.iteritems():
                assert len(rt_list) == 1
                rt = rt_list[0]
                for i in range(self.num_processors - 1):
                    rt_list.append(RejectionTeam(
                            temp_fs = self.temp_fs,
                            observed_path = rt.observed_path,
                            output_dir = rt.output_dir,
                            output_prefix = self.output_prefix,
                            num_taxon_pairs = self.num_taxon_pairs,
                            prior_paths = [],
                            model_indices = rt.model_indices,
                            summary_in_path = rt.summary_in_path,
                            num_posterior_samples = self.num_posterior_samples,
                            num_posterior_density_quantiles = \
                                    self.num_posterior_density_quantiles,
                            run_regression = False,
                            eureject_exe_path = self.eureject_exe_path,
                            abctoolbox_exe_path = self.abctoolbox_exe_path,
                            abctoolbox_bandwidth = self.abctoolbox_bandwidth,
                            omega_threshold = self.omega_threshold,
                            compress = self.compress,
                            keep_temps = self.keep_temps,
                            index = rt.index,
                            tag = model_idx))
    
    def _run_workers(self, workers):
        assert self.work_queue.empty()
        assert self.result_queue.empty()
        for w in workers:
            self.work_queue.put(w)
        managers = []
        for i in range(self.num_processors):
            m = Manager(work_queue = self.work_queue,
                    result_queue = self.result_queue)
            m.start()
            managers.append(m)
        for i in range(len(workers)):
            workers[i] = self.result_queue.get()
        for m in managers:
            m.join()
        assert self.work_queue.empty()
        assert self.result_queue.empty()
        return workers

    def _run_prior_workers(self, prior_worker_batch):
        return self._run_workers(prior_worker_batch)

    def _merge_summaries(self):
        summary_workers = dict(zip(self.models.keys(),
                [[] for i in range(len(self.models.keys()))]))
        sum_mergers = []
        summaries = {}
        for pw in self.prior_summary_workers:
            assert pw.tag == pw.summary_worker.tag
            summary_workers[pw.tag].append(pw.summary_worker)
        for model_idx, sum_workers in summary_workers.iteritems():
            summary_merger = EuRejectSummaryMerger(sum_workers,
                    tag = model_idx)
            sum_mergers.append(summary_merger)
        sum_mergers = self._run_workers(sum_mergers)
        for sm in sum_mergers:
            summaries[sm.tag] = sm.sample_sum_collection
        if len(self.models) > 1:
            summaries['combined'] = SampleSummaryCollection.merge(
                    summaries.itervalues())
        for k, summary in summaries.iteritems():
            summary.write(self.summary_paths[k])

    def _load_duplicate_rejection_teams(self, prior_worker_batch, model_idx):
        assert self.num_observed == 1
        paths = [pw.prior_path for pw in prior_worker_batch]
        if not paths:
            return
        split_paths = list(list_splitter(
                paths,
                len(self.rejection_teams[model_idx])))
        for i in range(len(split_paths)):
            self.rejection_teams[model_idx][i].prior_paths.extend(
                    split_paths[i])

    def _load_rejection_teams(self, prior_worker_batch):
        if self.global_estimate_only:
            if self.num_observed == 1:
                self._load_duplicate_rejection_teams(prior_worker_batch, 'combined')
                return
            for rt in self.rejection_teams['combined']:
                rt.prior_paths.extend(
                        [pw.prior_path for pw in prior_worker_batch])
            return
        model_indices = set(pw.tag for pw in prior_worker_batch)
        prior_workers = dict(zip(model_indices,
                [[] for i in range(len(model_indices))]))
        for pw in prior_worker_batch:
            prior_workers[pw.tag].append(pw)
        for model_idx, pw_list in prior_workers.iteritems():
            if self.num_observed == 1:
                self._load_duplicate_rejection_teams(pw_list, model_idx)
            else:
                for rt in self.rejection_teams[model_idx]:
                    rt.prior_paths.extend([pw.prior_path for pw in pw_list])
        
    def _run_rejection_teams(self):
        reject_teams = []
        for model_idx, rt_list in self.rejection_teams.iteritems():
            reject_teams.extend(rt_list)
        reject_teams = self._run_workers(reject_teams)
        rteams = dict(zip(self.rejection_teams.keys(),
                [[] for i in range(len(self.rejection_teams))]))
        for rt in reject_teams:
            rteams[rt.tag].append(rt)
        self.rejection_teams = rteams

    def _purge_priors(self, prior_workers):
        for pw in prior_workers:
            os.remove(pw.prior_path)
            
    def _merge_rejection_teams(self):
        if self.num_observed > 1:
            return
        paths_to_purge = []
        for model_idx, rt_list in self.rejection_teams.iteritems():
            for i in range(len(rt_list)-1):
                rt = rt_list.pop(-1)
                if rt.posterior_path:
                    rt_list[0].prior_paths.append(rt.posterior_path)
                    paths_to_purge.append(rt.posterior_path)
            assert len(rt_list) == 1
        self._run_rejection_teams()
        if paths_to_purge:
            for p in paths_to_purge:
                os.remove(p)

    def _sort_rejection_teams(self):
        for k in self.rejection_teams.iterkeys():
            self.rejection_teams[k] = sorted(self.rejection_teams[k],
                    key = lambda x : x.observed_path)

    def _run_final_rejections(self):
        if self.global_estimate_only:
            return
        if len(self.models) < 2:
            return
        self._sort_rejection_teams()
        rej_teams = {}
        rt_length = len(self.rejection_teams['combined'])
        for k in self.models.iterkeys():
            assert len(self.rejection_teams[k]) == rt_length
            for i in range(len(self.rejection_teams[k])):
                assert self.rejection_teams[k][i].observed_path == \
                       self.rejection_teams['combined'][i].observed_path
                self.rejection_teams['combined'][i].prior_paths.append(
                        self.rejection_teams[k][i].posterior_path)
        self.rejection_teams['combined'] = self._run_workers(
                self.rejection_teams['combined'])

    def _run_regressions(self):
        for model_idx, rt_list in self.rejection_teams.iteritems():
            for rt in rt_list:
                rt.run_regression = True
        self._run_rejection_teams()

    def run(self):
        _LOG.debug('running prior summary workers...')
        self.prior_summary_workers = self._run_prior_workers(
                self.prior_summary_workers)
        _LOG.debug('merging summaries...')
        self._merge_summaries()
        _LOG.debug('loading summary workers for rejection...')
        self._load_rejection_teams(self.prior_summary_workers)
        _LOG.debug('running rejection on summary worker priors...')
        self._run_rejection_teams()
        _LOG.debug('purging summary worker priors...')
        self._purge_priors(self.prior_summary_workers)
        for i in range(len(self.prior_workers)):
            _LOG.debug('running prior worker batch...')
            self.prior_workers[i] = self._run_prior_workers(
                    self.prior_workers[i])
            _LOG.debug('loading priors for rejection ...')
            self._load_rejection_teams(self.prior_workers[i])
            _LOG.debug('running rejection teams ...')
            self._run_rejection_teams()
            _LOG.debug('purging priors...')
            self._purge_priors(self.prior_workers[i])
        _LOG.debug('merging rejection team posteriors...')
        self._merge_rejection_teams()
        _LOG.debug('running final rejections...')
        self._run_final_rejections()
        _LOG.debug('running regressions...')
        self._run_regressions()
        self.finished = True


class RejectionTeam(object):
    count = 0
    def __init__(self,
            temp_fs,
            num_taxon_pairs,
            observed_path,
            output_dir,
            output_prefix = '',
            prior_paths = [],
            model_indices = None,
            summary_in_path = None,
            num_posterior_samples = 1000,
            num_posterior_density_quantiles = 1000,
            run_regression = False,
            eureject_exe_path = None,
            abctoolbox_exe_path = None,
            abctoolbox_bandwidth = None,
            omega_threshold = 0.01,
            compress = True,
            keep_temps = False,
            index = 1,
            tag = None):
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.temp_fs = temp_fs
        self.observed_path = expand_path(observed_path)
        self.output_dir = expand_path(output_dir)
        self.num_taxon_pairs = int(num_taxon_pairs)
        if model_indices:
            model_indices = list(set(model_indices))
        self.model_indices = model_indices
        self.num_posterior_samples = int(num_posterior_samples)
        self.num_posterior_density_quantiles = int(
                num_posterior_density_quantiles)
        if eureject_exe_path:
            eureject_exe_path = expand_path(eureject_exe_path)
        self.eureject_exe_path = eureject_exe_path
        self.abctoolbox_exe_path = abctoolbox_exe_path
        self.abctoolbox_bandwidth = abctoolbox_bandwidth
        self.omega_threshold = omega_threshold
        self.compress = compress
        self.prior_paths = [p for p in prior_paths]
        if summary_in_path:
            summary_in_path = expand_path(summary_in_path)
        self.summary_in_path = summary_in_path
        self.keep_temps = bool(keep_temps)
        self.keep_priors = True
        self.posterior_path = None
        self.run_regression = run_regression
        self.index = int(index)
        self.output_prefix = os.path.join(self.output_dir,
                str(output_prefix) + str(self.index))
        self.tag = tag
        self.div_model_results_path = None
        self.model_results_path = None
        self.psi_results_path = None
        self.omega_results_path = None
        self.posterior_summary_path = None
        self.regress_summary_path = None
        self.regress_posterior_path = None
        self.num_calls = 0

    def rejection_ready(self):
        if not self.summary_in_path or not os.path.exists(self.summary_in_path):
            return False
        if len(self.prior_paths) < 1:
            return False
        return True

    def regression_ready(self):
        if self.prior_paths:
            return False
        if not self.posterior_path:
            return False
        return True

    def start(self):
        if self.run_regression:
            if not self.regression_ready():
                return
            self._run_regression_workers()
        else:
            if not self.rejection_ready():
                return
            self._run_rejection_workers()

    def _run_rejection_workers(self):
        self.num_calls += 1
        p_paths = []
        old_posterior_path = None
        prior_paths = copy.deepcopy(self.prior_paths)
        for i in range(len(self.prior_paths)):
            path = self.prior_paths.pop(0)
            p_paths.append(path)
        if self.posterior_path:
            old_posterior_path = str(self.posterior_path)
            p_paths.append(old_posterior_path)
        self.posterior_path = self.temp_fs.get_file_path(
                prefix = 'posterior-{0}-{1}-{2}-'.format(self.name, 
                        self.num_calls, i+1),
                create = False)
        rw = EuRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = self.observed_path,
                prior_paths = p_paths,
                num_posterior_samples = self.num_posterior_samples,
                num_standardizing_samples = 0,
                summary_in_path = self.summary_in_path,
                summary_out_path = None,
                posterior_path = self.posterior_path,
                regression_worker = None,
                exe_path = self.eureject_exe_path,
                keep_temps = self.keep_temps,
                tag = self.tag)
        rw.start()
        if not self.keep_temps:
            if old_posterior_path:
                os.remove(old_posterior_path)
        if not self.keep_priors:
            for pp in prior_paths:
                os.remove(pp)
        assert len(self.prior_paths) == 0

    def _run_regression_workers(self):
        post_path = self.output_prefix + '-posterior-sample.txt'
        pw = PosteriorWorker(
                temp_fs = self.temp_fs,
                observed_path = self.observed_path,
                posterior_path = self.posterior_path,
                num_taxon_pairs = self.num_taxon_pairs,
                posterior_out_path = post_path,
                output_prefix = self.output_prefix,
                model_indices = self.model_indices,
                abctoolbox_exe_path = self.abctoolbox_exe_path,
                abctoolbox_bandwidth = self.abctoolbox_bandwidth,
                abctoolbox_num_posterior_quantiles = \
                        self.num_posterior_density_quantiles,
                omega_threshold = self.omega_threshold,
                compress = self.compress,
                keep_temps = self.keep_temps,
                tag = self.tag)
        pw.start()
        os.remove(self.posterior_path)
        self.posterior_path = post_path
        self.div_model_results_path = pw.div_model_results_path
        self.model_results_path = pw.model_results_path
        self.psi_results_path = pw.psi_results_path
        self.omega_results_path = pw.omega_results_path
        self.posterior_summary_path = pw.posterior_summary_path
        self.regress_summary_path = pw.regress_summary_path
        self.regress_posterior_path = pw.regress_posterior_path

