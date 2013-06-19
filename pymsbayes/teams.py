#! /usr/bin/env python

import os
import sys

from pymsbayes.workers import (MsBayesWorker, ABCToolBoxRejectWorker,
        ABCToolBoxRegressWorker, RegressionWorker, get_stat_indices,
        parse_header, DEFAULT_STAT_PATTERNS, PARAMETER_PATTERNS,
        merge_prior_files, MsRejectWorker, EuRejectSummaryMerger,
        EuRejectWorker)
from pymsbayes.manager import Manager
from pymsbayes.fileio import process_file_arg, expand_path
from pymsbayes.utils import GLOBAL_RNG, WORK_FORCE
from pymsbayes.utils.functions import (long_division, least_common_multiple,
        get_random_int, list_splitter)
from pymsbayes.utils.stats import (merge_sample_summary_mappings
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class ABCTeam(object):
    count = 0
                        # num_taxon_pairs = self.num_taxon_pairs,
                        # num_posterior_density_quantiles = \
                        #         self.num_posterior_density_quantiles,
                        # abctoolbox_exe_path = self.abctoolbox_exe_path,
                        # abctoolbox_bandwidth = self.abctoolbox_bandwidth,
                        # omega_threshold = self.omega_threshold,
                        # compress = self.compress,
    def __init__(self,
            temp_fs,
            observed_sims_file,
            num_taxon_pairs,
            model_indices_to_config_paths,
            num_prior_samples,
            num_processors,
            num_standardizing_samples = 10000,
            num_posterior_samples = 1000,
            batch_size = 10000,
            output_dir = None,
            rng = None,
            sort_index = None,
            report_parameters = True,
            stat_patterns = DEFAULT_STAT_PATTERNS,
            parameter_patterns = PARAMETER_PATTERNS,
            msbayes_exe_path = None,
            abctoolbox_exe_path = None,
            eureject_exe_path = None,
            keep_temps = False,
            global_estimate_only = False,
            work_queue = WORK_FORCE):
        if not rng:
            rng = GLOBAL_RNG
        self.rng = rng
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.temp_fs = temp_fs
        self.observed_temp_dir = temp_fs.create_subdir(
                prefix = 'observed-files-')
        self.summary_temp_dir = temp_fs.create_subdir(
                prefix = 'summary-files-')
        self.observed_sims_path = expand_path(observed_sims_path)
        if not output_dir:
            output_dir = os.path.dirname(self.observed_sims_file)
        self.output_dir = mk_new_dir(os.path.join(expand_path(output_dir),
                'pymsbayes-output'))
        self.global_estimate_only = global_estimate_only
        self.models = model_indices_to_config_paths
        self.num_prior_samples = num_prior_samples
        self.num_processors = num_processors
        self.num_standardizing_samples = num_standardizing_samples
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
        self.eureject_exe_path = eureject_exe_path
        self.model_dirs = dict(zip(self.models.keys(),
                [mk_new_dir(os.path.join(self.output_dir,
                        'm' + str(k))) for k in self.models.keys()]))
        self.model_dirs.update({'combined': mk_new_dir(os.path.join(
                self.output_dir,
                'm' + ''.join([i for i in sorted(self.models.keys())])))})
        self.summary_paths = {}
        for k, d in self.model_dirs.iteritems():
            self.summary_paths[k] = os.path.join(d,
                    os.path.basename(d) + '-stat-means-and-std-devs.txt')
        self.rejection_teams = {}
        if not self.global_estimate_only:
            self.rejection_teams = dict(zip(self.models.keys(),
                    [[] for i in range(len(self.models.keys()))]))
        self.rejection_teams.update({'combined': []})
        self.prior_summary_workers = []
        self.prior_workers = []
        self.num_observed = 0
        self.keep_temps = keep_temps
        self.work_queue = work_queue
        self.result_queue = multiprocessing.Queue()
        self._assemble_rejection_teams()
        self._assemble_prior_workers()

    def _assemble_prior_workers(self):
        prior_workers = []
        for model_index, config_path in self.models:
            to_summarize = self.num_standardizing_samples
            for i in range(self.num_batches):
                sum_worker = None
                if to_summarize > 0:
                    if (to_summarize - self.batch) >= 0:
                        nsum = self.batch_size
                    else:
                        nsum = to_summarize
                    to_summarize -= nsum
                    sum_path = temp_fs.get_file_path(
                            parent = self.summary_temp_dir,
                            prefix = 'summary-{0}-{1}-'.format(model_index,
                                    i+1),
                            create = False)
                    sum_worker = EuRejectWorker(
                            temp_fs = self.temp_fs,
                            observed_path = self.observed_sims_path,
                            prior_paths = [worker.prior_path],
                            num_posterior_samples = 0,
                            num_standardizing_samples = nsum,
                            summary_out_path = sum_path,
                            exe_path = self.eureject_exe_path,
                            keep_temps = self.keep_temps,
                            tag = str(model_index))
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
                        summary_worker = sum_worker
                        tag = str(model_index))
                if sum_worker:
                    self.prior_summary_workers.append(worker)
                else:
                    prior_workers.append(worker)
            if self.extra_samples > 0:
                sum_worker = None
                if to_summarize > 0:
                    if (to_summarize - self.extra_samples) >= 0:
                        nsum = self.extra_samples
                    else:
                        nsum = to_summarize
                    to_summarize -= nsum
                    sum_path = temp_fs.get_file_path(
                            parent = self.summary_temp_dir,
                            prefix = 'summary-{0}-{1}-'.format(model_index,
                                    i+1),
                            create = False)
                    sum_worker = EuRejectWorker(
                            temp_fs = self.temp_fs,
                            observed_path = self.observed_sims_path,
                            prior_paths = [worker.prior_path],
                            num_posterior_samples = 0,
                            num_standardizing_samples = nsum,
                            summary_out_path = sum_path,
                            exe_path = self.eureject_exe_path,
                            keep_temps = self.keep_temps,
                            tag = str(model_index))
                    self.summary_workers.append(sum_worker)
                worker = MsBayesWorker(
                        temp_fs = self.temp_fs,
                        sample_size = self.extra_samples,
                        config_path = config_path,
                        exe_path = self.msbayes_exe_path,
                        model_index = model_index,
                        report_parameters = self.report_parameters,
                        seed = seed,
                        schema = 'abctoolbox',
                        stat_patterns = self.stat_patterns,
                        summary_worker = sum_worker
                        tag = str(model_index))
                if sum_worker:
                    self.prior_summary_workers.append(worker)
                else:
                    prior_workers.append(worker)
        self.prior_workers = list(list_splitter(prior_workers,
                self.num_processors))

    def _assemble_rejection_teams(self):
        obs_file, close = process_file_arg(self.observed_sims_path)
        header = parse_header(obs_file)
        all_stat_indices = get_stat_indices(header,
                stat_patterns=ALL_STAT_PATTERNS)
        obs_file.next() # header line
        for i, line in enumerate(obs_file):
            self.num_observed += 1
            index = i + 1
            tmp_obs_path = temp_fs.get_file_path(
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
                    m_indices = models.keys()
                else:
                    m_indices = [model_idx]
                rejection_team_list.append(EuRejectTeam(
                        temp_fs = self.temp_fs,
                        observed_path = tmp_obs_path,
                        output_dir = self.model_dirs[model_idx],
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
                        tag = str(model_idx)))
        if close:
            obs_file.close()
    
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

    def _run_prior_workers(self, prior_worker_batch):
        self._run_workers(prior_worker_batch)

    def _merge_summaries(self):
        summary_workers = dict(zip(self.models.keys(),
                [[] for i in range(len(self.models.keys()))]))
        sum_mergers = []
        for pw in self.prior_summary_workers:
            assert pw.tag == pw.summary_worker.tag
            summary_workers[pw.tag] = pw.summary_worker
        for model_idx, sum_workers in summary_workers.iteritems():
            summary_merger = EuRejectSummaryMerger(sum_workers,
                    tag = str(model_idx))
            sum_mergers.append(summary_merger)
        self._run_workers(sum_mergers)
        for sm in sum_mergers:
            self.summaries[sm.tag] = sm.sample_sum_collection
        self.summaries['combined'] = SampleSummaryCollection.merge(
                self.summaries.itervalues())
        for k, summary in self.summaries.iteritems():
            summary.write(self.summary_paths[k])

    def _load_rejection_teams(self, prior_worker_batch):
        teams_to_run = []
        if self.global_estimate_only:
            for rt in self.rejection_teams['combined']:
                rt.prior_paths.extend(
                        [pw.prior_path for pw in prior_worker_batch])
                teams_to_run.append(rt)
            return teams_to_run
        model_indices = set(pw.tag for pw in prior_worker_batch)
        prior_workers = dict(zip(model_indices,
                [[] for i in range(len(model_indices))]))
        for pw in prior_worker_batch:
            prior_workers[pw.tag].append(pw)
        for model_idx, pw_list in prior_workers.iteritems():
            for rt in self.rejection_teams[model_idx]:
                rt.prior_paths.extend([pw.prior_path for pw in pw_list])
                teams_to_run.append(rt)
        return teams_to_run
        
    def _run_rejection_teams(self, prior_worker_batch):
        teams_to_run = self._load_rejection_teams(prior_worker_batch)
        self._run_workers(teams_to_run)

    def _run_final_rejections(self):
        if self.global_estimate_only:
            return
        rej_teams = {}
        rt_length = length(self.rejection_teams['combined'])
        for k in self.models_keys():
            assert len(self.rejection_teams[k]) == rt_length
            for i in range(len(self.rejection_teams[k])):
                assert self.rejection_teams[k][i].observed_path == \
                       self.rejection_teams['combined'][i].observed_path
                self.rejection_teams['combined'][i].prior_paths.append(
                        self.rejection_teams[k][i].posterior_path)
        self._run_workers(self.rejection_teams['combined'])

    def _run_regressions(self):
        to_run = []
        for model_idx, rt_list in self.rejection_teams.iteritems():
            for rt in rt_list:
                rt.run_regression = True
                to_run.append(rt)
        self._run_workers(to_run)

    def run(self):
        self._run_prior_workers(self.prior_summary_workers)
        self._merge_summaries()
        self._run_rejection_teams(self.prior_summary_workers)
        for prior_worker_batch in self.prior_workers:
            self._run_prior_workers(prior_worker_batch)
            self._run_rejection_teams(prior_worker_batch)
        self._run_final_rejections()
        self._run_regressions()


class EuRejectTeam(object):
    count = 0
    def __init__(self,
            temp_fs,
            observed_path,
            output_dir,
            num_taxon_pairs,
            prior_paths = [],
            model_indices = None
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
            tag = ''):
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
        self.posterior_path = None
        self.run_regression = run_regression
        self.index = int(index)
        self.output_prefix = os.path.join(self.output_dir, str(self.index))
        self.num_calls = 0

    def _check_rejection_ready(self):
        if not summary_in_path or not os.path.exists(summary_in_path):
            raise Exception('{0} started without a summary path'.format(
                    self.name))
        if len(prior_paths) < 1:
            raise Exception('{0} started without priors'.format(
                    self.name))

    def _check_regression_ready(self):
        if not self.posterior_path:
            raise Exception('regression called for {0}, but there is no '
                    'posterior yet'.format(self.name))

    def start(self):
        if self.run_regression:
            self._run_regression_workers()
        else:
            self._run_rejection_workers()

    def _run_rejection_workers(self):
        self._check_rejection_ready():
        self.num_calls += 1
        p_paths = []
        for i in range(len(self.prior_paths)):
            path = self.prior_paths.pop(0)
            p_paths.append(path)
        if self.posterior_path:
            p_paths.append(str(self.posterior_path))
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
                posterior_path = new_posterior_path,
                regression_worker = None,
                exe_path = self.eureject_exe_path,
                keep_temps = self.keep_temps,
                tag = self.tag)
        rw.start()
        if not self.keep_temps:
            for pp in p_paths:
                os.remove(pp)
        assert len(self.prior_paths) == 0

    def _run_regression_workers(self):
        self._check_regression_ready(self)
        post_path = self.output_prefix + '-posterior-sample.txt.gz'
        pw = PosteriorWorker(
                temp_fs = self.temp_fs,
                observed_path = self.observed_path,
                posterior_path = self.posterior_path,
                num_taxon_pairs = self.num_taxon_pairs,
                posterior_out_path = post_path,
                output_prefix = out_prefix,
                model_indices = self.model_indices,
                abctoolbox_exe_path = self.abctoolbox_exe_path,
                abctoolbox_bandwidth = self.abctoolbox_bandwidth,
                abctoolbox_num_posterior_quantiles = \
                        self.num_posterior_density_quantiles
                omega_threshold = self.omega_threshold,
                compress = self.compress,
                keep_temps = self.keep_temps,
                tag = self.tag)
        pw.start()
        shutil.remove(self.posterior_path)
        self.posterior_path = post_path


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
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.temp_fs = temp_fs
        self.observed_path = observed_path
        self.num_posterior_samples = num_posterior_samples
        self.abctoolbox_exe_path = abctoolbox_exe_path
        self.prior_workers = prior_workers
        self.keep_temps = keep_temps
        self.posterior_path = None
        self.num_calls = 0

    # 2/10
    # def start(self):
    #     self.num_calls += 1
    #     for i in range(len(self.prior_workers)):
    #         # pw = self.prior_workers.pop(0)
    #         pw = self.prior_workers[i]
    #         temp_posterior_path = self.temp_fs.get_file_path(
    #                 prefix = 'temp-posterior-{0}-{1}-{2}-'.format(self.name, 
    #                         self.num_calls, i+1),
    #                 create = False)
    #         rw = ABCToolBoxRejectWorker(
    #                 temp_fs = self.temp_fs,
    #                 observed_path = self.observed_path,
    #                 prior_path = pw.prior_path,
    #                 num_posterior_samples = self.num_posterior_samples,
    #                 posterior_path = temp_posterior_path,
    #                 regression_worker = None,
    #                 exe_path = self.abctoolbox_exe_path,
    #                 keep_temps = self.keep_temps,
    #                 max_read_sims = pw.sample_size + 10)
    #         rw.start()
    #         if self.posterior_path:
    #             new_prior_path = self.temp_fs.get_file_path(
    #                     prefix = 'prior-{0}-{1}-{2}-'.format(self.name, 
    #                             self.num_calls, i+1),
    #                     create = False)
    #             merge_prior_files(
    #                     paths = [self.posterior_path, temp_posterior_path],
    #                     dest_path = new_prior_path)
    #             new_posterior_path = self.temp_fs.get_file_path(
    #                     prefix = 'posterior-{0}-{1}-{2}-'.format(self.name, 
    #                             self.num_calls, i+1),
    #                     create = False)
    #             if not self.keep_temps:
    #                 os.remove(temp_posterior_path)
    #                 os.remove(self.posterior_path)
    #                 os.remove(rw.prior_path)
    #             rw = ABCToolBoxRejectWorker(
    #                     temp_fs = self.temp_fs,
    #                     observed_path = self.observed_path,
    #                     prior_path = new_prior_path,
    #                     num_posterior_samples = self.num_posterior_samples,
    #                     posterior_path = new_posterior_path,
    #                     regression_worker = None,
    #                     exe_path = self.abctoolbox_exe_path,
    #                     keep_temps = self.keep_temps,
    #                     max_read_sims = (2 * self.num_posterior_samples) + 10)
    #             rw.start()
    #             self.posterior_path = new_posterior_path
    #         else:
    #             self.posterior_path = temp_posterior_path
    #0/10
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
                    max_read_sims = pw.sample_size + \
                            self.num_posterior_samples + 10)
            rw.start()
            if not self.keep_temps:
                os.remove(rw.prior_path)

class MsRejectionTeam(object):
    count = 0
    def __init__(self,
            temp_fs,
            prior_workers,
            observed_path,
            num_posterior_samples = 1000,
            msreject_exe_path = None,
            keep_temps = False):
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.temp_fs = temp_fs
        self.observed_path = observed_path
        self.num_posterior_samples = num_posterior_samples
        self.msreject_exe_path = msreject_exe_path
        self.prior_workers = prior_workers
        self.keep_temps = keep_temps
        self.posterior_path = None
        self.num_calls = 0

    #6/10
    def start(self):
        self.num_calls += 1
        for i in range(len(self.prior_workers)):
            # pw = self.prior_workers.pop(0)
            pw = self.prior_workers[i]
            nsamples = pw.sample_size
            merge_paths = [pw.prior_path]
            if self.posterior_path:
                nsamples += self.num_posterior_samples
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
            tolerance = self.num_posterior_samples / float(nsamples)
            rw = MsRejectWorker(
                    header = pw.header,
                    observed_path = self.observed_path,
                    prior_path = new_prior_path,
                    tolerance = tolerance,
                    posterior_path = new_posterior_path,
                    regression_worker = None,
                    exe_path = self.msreject_exe_path)
            rw.start()
            if not self.keep_temps:
                os.remove(rw.prior_path)
    #1/10
    # def start(self):
    #     self.num_calls += 1
    #     for i in range(len(self.prior_workers)):
    #         # pw = self.prior_workers.pop(0)
    #         pw = self.prior_workers[i]
    #         temp_posterior_path = self.temp_fs.get_file_path(
    #                 prefix = 'temp-posterior-{0}-{1}-{2}-'.format(self.name, 
    #                         self.num_calls, i+1),
    #                 create = False)
    #         tolerance = self.num_posterior_samples /float(pw.sample_size)
    #         rw = MsRejectWorker(
    #                 header = pw.header,
    #                 observed_path = self.observed_path,
    #                 prior_path = pw.prior_path,
    #                 tolerance = tolerance,
    #                 posterior_path = temp_posterior_path,
    #                 regression_worker = None,
    #                 exe_path = self.msreject_exe_path)
    #         rw.start()
    #         if self.posterior_path:
    #             new_prior_path = self.temp_fs.get_file_path(
    #                     prefix = 'prior-{0}-{1}-{2}-'.format(self.name, 
    #                             self.num_calls, i+1),
    #                     create = False)
    #             merge_prior_files(
    #                     paths = [self.posterior_path, temp_posterior_path],
    #                     dest_path = new_prior_path)
    #             new_posterior_path = self.temp_fs.get_file_path(
    #                     prefix = 'posterior-{0}-{1}-{2}-'.format(self.name, 
    #                             self.num_calls, i+1),
    #                     create = False)
    #             if not self.keep_temps:
    #                 os.remove(temp_posterior_path)
    #                 os.remove(self.posterior_path)
    #                 os.remove(rw.prior_path)
    #             rw = MsRejectWorker(
    #                     header = pw.header,
    #                     observed_path = self.observed_path,
    #                     prior_path = new_prior_path,
    #                     tolerance = 0.5,
    #                     posterior_path = new_posterior_path,
    #                     regression_worker = None,
    #                     exe_path = self.msreject_exe_path)
    #             rw.start()
    #             self.posterior_path = new_posterior_path
    #         else:
    #             self.posterior_path = temp_posterior_path

