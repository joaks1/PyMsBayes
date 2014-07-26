#! /usr/bin/env python

import os
import sys
import re
import copy
import shutil
import math
import multiprocessing
import logging
import glob
import Queue

from pymsbayes import config
from pymsbayes import workers
from pymsbayes import manager
from pymsbayes import fileio
from pymsbayes.utils import (GLOBAL_RNG, WORK_FORCE, DUMP_DEBUG_INFO,
        dump_memory_info)
from pymsbayes.utils import parsing
from pymsbayes.utils import sumresults
from pymsbayes.utils import stats
from pymsbayes.utils import probability
from pymsbayes.utils import functions
from pymsbayes.utils import stats
from pymsbayes.utils import tempfs
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class ABCTeam(object):
    count = 0
    summary_file_pattern = re.compile(r'^.*m(?P<model_index>\d+)-?'
            '(?P<combined>combined)?-stat-means-and-std-devs\.txt$')
    prior_file_pattern = re.compile(r'^.*m(?P<model_index>\d+)-?'
            '(?P<combined>combined)?-prior-sample\.txt$')

    def __init__(self,
            temp_fs,
            num_taxon_pairs,
            observed_stats_files,
            num_processors,
            config_paths = None,
            previous_prior_dir = None,
            num_prior_samples = 1000000,
            num_standardizing_samples = 10000,
            num_posterior_samples = 1000,
            num_posterior_density_quantiles = 1000,
            batch_size = 10000,
            output_dir = None,
            output_prefix = '',
            prior_temp_dir = None,
            rng = None,
            report_parameters = True,
            stat_patterns = parsing.DEFAULT_STAT_PATTERNS,
            eureject_exe_path = None,
            abctoolbox_exe_path = None,
            msbayes_exe_path = None,
            abctoolbox_bandwidth = None,
            omega_threshold = 0.01,
            compress = True,
            keep_temps = False,
            reporting_frequency = None,
            global_estimate = True,
            global_estimate_only = False,
            generate_prior_samples_only = False,
            start_from_observed_index = 0,
            start_from_simulation_index = 0):
        if ((config_paths) and (previous_prior_dir)) or \
                ((not config_paths) and (not previous_prior_dir)):
            raise ValueError('`config_paths` or `previous_prior_dir` must '
                    'be provided, but not both')
        if global_estimate_only and (not global_estimate):
            raise ValueError('`global_estimate_only` cannot be true if '
                    '`global_estimate` is false`')
        if previous_prior_dir and generate_prior_samples_only:
            raise ValueError('`generate_prior_samples_only` must be false '
                    'when `previous_prior_dir` is provided')
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
        self.posterior_temp_dir = temp_fs.create_subdir(
                parent = self.temp_output_dir,
                prefix = 'posterior-files-')
        self.old_posterior_temp_dir = temp_fs.create_subdir(
                parent = self.temp_output_dir,
                prefix = 'old-posterior-files-')
        self.prior_temp_fs = self.temp_fs
        if prior_temp_dir:
            self.prior_temp_fs = tempfs.TempFileSystem(parent = prior_temp_dir,
                    prefix = 'pymsbayes-temp-priors-')
        if not output_dir:
            output_dir = os.path.dirname(observed_stats_files[0])
        self.output_dir = functions.mk_new_dir(os.path.join(
                fileio.expand_path(output_dir),
                'pymsbayes-output'))
        self.output_prefix = str(output_prefix)
        self.global_estimate = global_estimate
        self.global_estimate_only = global_estimate_only
        self.num_taxon_pairs = num_taxon_pairs
        self.generate_prior_samples_only = generate_prior_samples_only
        self.observed_stats_paths = dict(zip(
                [i + 1 for i in range(len(observed_stats_files))],
                [fileio.expand_path(p) for p in observed_stats_files]))
        self.observed_dirs = dict(zip(
                self.observed_stats_paths.keys(),
                [functions.mk_new_dir(os.path.join(self.output_dir, 'd' + str(
                        k))) for k in self.observed_stats_paths.iterkeys()]))
        self.report_parameters = report_parameters
        self.stat_patterns = stat_patterns
        self.eureject_exe_path = eureject_exe_path
        self.abctoolbox_exe_path = abctoolbox_exe_path
        self.msbayes_exe_path = msbayes_exe_path
        self.abctoolbox_bandwidth = abctoolbox_bandwidth
        self.omega_threshold = omega_threshold
        self.compress = compress
        self.num_processors = num_processors
        self.num_posterior_samples = num_posterior_samples
        self.num_posterior_density_quantiles = num_posterior_density_quantiles
        self.start_from_simulation = start_from_simulation_index
        self.start_from_observed = start_from_observed_index

        self.models = None
        self.model_strings = None
        self.model_samples = None
        self.summary_dir = None
        self.summary_paths = None
        self.previous_prior_dir = None
        self.use_previous_priors = False
        if previous_prior_dir:
            self.use_previous_priors = True
            self.previous_prior_dir = fileio.expand_path(previous_prior_dir)
            self._parse_previous_prior_dir()
            self.num_prior_samples = 0
            self.num_standardizing_samples = 0
            self.batch_size = 0
            self.total_samples = 0
            self.num_batches = 0
            self.num_extra_samples = 0
            self.num_batches_remaining = 0
            self.num_prior_batch_iters = 1
            self.reporting_frequency = None
            self.reporting_indices = None
        else:
            self.models = dict(zip(
                    [i + 1 for i in range(len(config_paths))],
                    [fileio.expand_path(p) for p in config_paths]))
            self.model_strings = dict(zip(self.models.keys(),
                    ['m' + str(i) for i in self.models.iterkeys()]))
            if len(self.models) > 1:
                self.model_strings['combined'] = 'm' + ''.join(
                        [str(i) for i in sorted(self.models.iterkeys())]) + \
                                '-combined'
            self.summary_dir = functions.mk_new_dir(os.path.join(
                    self.output_dir,
                    'prior-stats-summaries'))
            self.summary_paths = {}
            for k, v in self.model_strings.iteritems():
                self.summary_paths[k] = os.path.join(self.summary_dir,
                        self.output_prefix + v + '-stat-means-and-std-devs.txt')
            self.num_prior_samples = num_prior_samples
            self.num_standardizing_samples = num_standardizing_samples
            if self.num_standardizing_samples > self.num_prior_samples:
                self.num_standardizing_samples = self.num_prior_samples
            self.batch_size = batch_size
            self.total_samples = self.num_prior_samples * len(self.models)
            max_samples_per_process = self.total_samples / self.num_processors
            if max_samples_per_process < self.batch_size:
                self.batch_size = max_samples_per_process
            if self.batch_size < self.num_posterior_samples:
                self.batch_size = self.num_posterior_samples
            self.num_batches, self.num_extra_samples = functions.long_division(
                    self.num_prior_samples, 
                    self.batch_size)
            self.num_batches_remaining, self.num_extra_samples_remaining = \
                    self.num_batches, self.num_extra_samples
            self.num_prior_batch_iters = int(math.ceil(
                    (self.num_prior_samples - (int(math.ceil(
                        self.num_standardizing_samples / float(self.batch_size))) \
                                * self.batch_size)) / 
                    float(self.batch_size * self.num_processors)
                ))
            if self.num_prior_batch_iters < 1:
                self.num_prior_batch_iters = 1
            self.reporting_frequency = reporting_frequency
            self.reporting_indices = None
            if self.reporting_frequency and self.reporting_frequency > 0:
                rep_indices = range(0, self.num_prior_batch_iters,
                        self.reporting_frequency)
                if rep_indices[-1] >= (self.num_prior_batch_iters - 1):
                    rep_indices.pop(-1)
                self.reporting_indices = rep_indices
        self.model_dirs = {}
        for obs_idx, obs_dir in self.observed_dirs.iteritems():
            self.model_dirs[obs_idx] = {}
            for model_idx, model_str in self.model_strings.iteritems():
                self.model_dirs[obs_idx][model_idx] = functions.mk_new_dir(
                        os.path.join(obs_dir, model_str))
        if len(self.models) < 2:
            self.global_estimate_only = False
            self.global_estimate = False
        self.num_observed = 0
        self.n_observed_map = dict(zip(self.observed_stats_paths.keys(),
                [0 for i in self.observed_stats_paths.iterkeys()]))
        self.num_samples_processed = None
        self.keep_temps = keep_temps
        self.duplicate_rejection_workers = False
        self.num_samples_generated = 0
        self.num_samples_summarized = 0
        self.model_key_path = os.path.join(self.output_dir,
                self.output_prefix + 'model-key.txt')
        self.data_key_path = os.path.join(self.output_dir,
                self.output_prefix + 'data-key.txt')
        self._write_keys()
        self.prior_summary_worker_iter = self._prior_summary_worker_iter()
        self.prior_worker_iter = self._prior_worker_iter()
        self.iter_count = 0
        self.finished = False

    def _parse_previous_prior_dir(self):
        self.generate_prior_samples_only = False
        self.summary_dir = self.previous_prior_dir
        self.summary_paths = {}
        model_key_pattern = os.path.abspath(os.path.join(
                os.path.dirname(self.summary_dir),
                '*model-key.txt'))
        model_key_path = glob.glob(model_key_pattern)
        assert len(model_key_path) == 1
        model_key_path = model_key_path[0]
        self.models = sumresults.parse_model_key_file(model_key_path)
        self.model_samples = {}
        self.model_strings = {}
        for f in os.listdir(self.previous_prior_dir):
            m1 = self.prior_file_pattern.match(f)
            m2 = self.summary_file_pattern.match(f)
            if m1:
                if m1.group('combined'):
                    raise Exception('unexpected file {0} in previous prior '
                            'directory {1}'.format(f, self.previous_prior_dir))
                model_idx = int(m1.group('model_index'))
                if self.model_samples.has_key(model_idx):
                    raise Exception('model index {0} found more than once in '
                            'previous prior directory {1}'.format(model_idx,
                                self.previous_prior_dir))
                self.model_samples[model_idx] = os.path.join(
                        self.previous_prior_dir,
                        f)
                self.model_strings[model_idx] = 'm' + str(model_idx)
            if m2:
                if m2.group('combined'):
                    if self.summary_paths.has_key('combined'):
                        raise Exception('multiple combined summary files '
                                'found in previous prior directory {0}'.format(
                                    self.previous_prior_dir))
                    self.model_strings['combined'] = 'm' + \
                            m2.group('model_index') + '-combined'
                    self.summary_paths['combined'] = os.path.join(
                            self.previous_prior_dir, f)
                else:
                    model_idx = int(m2.group('model_index'))
                    if self.summary_paths.has_key(model_idx):
                        raise Exception('summary index {0} found more than '
                                'once in previous prior directory {1}'.format(
                                    model_idx, self.previous_prior_dir))
                    self.summary_paths[model_idx] = os.path.join(
                            self.previous_prior_dir, f)
        expected_combined_str = 'm' + ''.join(
                [str(i) for i in sorted(self.model_samples.iterkeys())]) + \
                        '-combined'
        if expected_combined_str != \
                self.model_strings.get('combined', expected_combined_str):
            raise Exception('expecting combined string {0}, but found {1} '
                    'in previous prior directory {2}'.format(
                        expected_combined_str,
                        self.model_strings['combined'],
                        self.previous_prior_dir))
        mkeys = self.model_samples.keys()
        if sorted(mkeys) != sorted(self.models.keys()):
            raise Exception('problem parsing info from previous prior '
                    'directory {0} and model key file {1}'.format(
                        self.previous_prior_dir, model_key_path))
        if self.model_strings.has_key('combined'):
            mkeys.append('combined')
        if (sorted(mkeys) != sorted(self.model_strings.keys())) or (sorted(
                    self.model_strings.keys()) != sorted(
                        self.summary_paths.keys())):
            raise Exception('problem parsing info from previous prior '
                    'directory {0}'.format(self.previous_prior_dir))

    def get_observed_path(self, observed_index, simulation_index):
        return os.path.join(self.observed_temp_dir, 
                '{0}-d{1}-s{2}-observed.txt'.format(self.temp_fs.token_id,
                        observed_index, simulation_index + 1))

    def get_posterior_paths(self, observed_index, simulation_index,
            per_processor = False):
        if per_processor:
            d = dict(zip(self.model_strings.keys(),
                [[] for i in range(len(self.model_strings))]))
            for model_idx, paths in d.iteritems():
                for i in range(self.num_processors):
                    paths.append(os.path.join(self.posterior_temp_dir,
                            '{0}-d{1}-{2}-s{3}-{4}-posterior.txt'.format(
                                self.temp_fs.token_id,
                                observed_index,
                                self.model_strings[model_idx],
                                simulation_index + 1,
                                i)))
            return d
        return dict(zip(self.model_strings.keys(),
                [[os.path.join(self.posterior_temp_dir,
                        '{0}-d{1}-{2}-s{3}-posterior.txt'.format(
                                self.temp_fs.token_id,
                                observed_index,
                                self.model_strings[k],
                                simulation_index + 1
                        ))] for k in self.model_strings.iterkeys()]))

    def get_observed_posterior_paths(self, observed_index, simulation_index,
            per_processor = False):
        return (self.get_observed_path(observed_index, simulation_index),
                self.get_posterior_paths(observed_index, simulation_index,
                        per_processor = per_processor))

    def _get_prior_path(self, model_index):
        return os.path.join(os.path.dirname(self.summary_paths[model_index]),
                    ''.join([self.output_prefix, self.model_strings[model_index],
                            '-prior-sample.txt']))

    def _prior_summary_worker_iter(self):
        to_summarize = self.num_standardizing_samples
        i = -1
        for i in range(self.num_batches):
            if to_summarize <= 0:
                break
            if (to_summarize - self.batch_size) >= 0:
                nsum = self.batch_size
            else:
                nsum = to_summarize
            to_summarize -= nsum
            for model_index, config_path in self.models.iteritems():
                prior_path = self.prior_temp_fs.get_file_path(
                        parent = self.prior_temp_fs.base_dir,
                        prefix = 'prior-to-sum-{0}-{1}-'.format(
                                self.model_strings[model_index],
                                i + 1),
                        register = False)
                sum_path = self.temp_fs.get_file_path(
                        parent = self.summary_temp_dir,
                        prefix = 'summary-{0}-{1}-'.format(
                                model_index,
                                i + 1),
                        register = False)
                post_path = self.temp_fs.get_file_path(
                        parent = self.summary_temp_dir,
                        prefix = 'posterior-{0}-{1}-'.format(
                                model_index,
                                i + 1),
                        register = False)
                sum_worker = workers.EuRejectWorker(
                        temp_fs = self.temp_fs,
                        observed_path = self.get_observed_path(
                                observed_index = min(
                                        self.observed_stats_paths.iterkeys()),
                                simulation_index = 0),
                        prior_paths = [prior_path],
                        posterior_path = post_path,
                        num_posterior_samples = 0,
                        num_standardizing_samples = nsum,
                        summary_out_path = sum_path,
                        exe_path = self.eureject_exe_path,
                        keep_temps = self.keep_temps,
                        tag = model_index)
                seed = functions.get_random_int(self.rng)
                worker = workers.MsBayesWorker(
                        temp_fs = self.prior_temp_fs,
                        sample_size = self.batch_size,
                        config_path = config_path,
                        prior_path = prior_path,
                        exe_path = self.msbayes_exe_path,
                        model_index = model_index,
                        report_parameters = self.report_parameters,
                        seed = seed,
                        schema = 'abctoolbox',
                        stat_patterns = self.stat_patterns,
                        summary_worker = sum_worker,
                        write_header_file = False,
                        tag = model_index)
                yield worker
            self.num_batches_remaining -= 1
        if (to_summarize > 0) and (self.num_extra_samples > 0):
            if (to_summarize - self.num_extra_samples) >= 0:
                nsum = self.num_extra_samples
            else:
                nsum = to_summarize
            to_summarize -= nsum
            self.num_extra_samples_remaining = 0
            for model_index, config_path in self.models.iteritems():
                prior_path = self.prior_temp_fs.get_file_path(
                        parent = self.prior_temp_fs.base_dir,
                        prefix = 'prior-to-sum-{0}-{1}-'.format(
                                self.model_strings[model_index],
                                i + 2),
                        register = False)
                sum_path = self.temp_fs.get_file_path(
                        parent = self.summary_temp_dir,
                        prefix = 'summary-{0}-{1}-'.format(
                                model_index,
                                i + 2),
                        register = False)
                post_path = self.temp_fs.get_file_path(
                        parent = self.summary_temp_dir,
                        prefix = 'posterior-{0}-{1}-'.format(
                                model_index,
                                i + 2),
                        register = False)
                sum_worker = workers.EuRejectWorker(
                        temp_fs = self.temp_fs,
                        observed_path = self.get_observed_path(
                                observed_index = min(
                                        self.observed_stats_paths.iterkeys()),
                                simulation_index = 0),
                        prior_paths = [prior_path],
                        posterior_path = post_path,
                        num_posterior_samples = 0,
                        num_standardizing_samples = nsum,
                        summary_out_path = sum_path,
                        exe_path = self.eureject_exe_path,
                        keep_temps = self.keep_temps,
                        tag = model_index)
                seed = functions.get_random_int(self.rng)
                worker = workers.MsBayesWorker(
                        temp_fs = self.prior_temp_fs,
                        sample_size = self.num_extra_samples,
                        config_path = config_path,
                        prior_path = prior_path,
                        exe_path = self.msbayes_exe_path,
                        model_index = model_index,
                        report_parameters = self.report_parameters,
                        seed = seed,
                        schema = 'abctoolbox',
                        stat_patterns = self.stat_patterns,
                        summary_worker = sum_worker,
                        write_header_file = False,
                        tag = model_index)
                yield worker

    def _prior_worker_iter(self):
        n = 0
        prior_workers = []
        i = -1
        for i in range(self.num_batches_remaining):
            n += 1
            for model_index, config_path in self.models.iteritems():
                prior_path = self.prior_temp_fs.get_file_path(
                        parent = self.prior_temp_fs.base_dir,
                        prefix = 'prior-{0}-{1}-'.format(
                                self.model_strings[model_index],
                                i + 1),
                        register = False)
                seed = functions.get_random_int(self.rng)
                worker = workers.MsBayesWorker(
                        temp_fs = self.prior_temp_fs,
                        sample_size = self.batch_size,
                        config_path = config_path,
                        prior_path = prior_path,
                        exe_path = self.msbayes_exe_path,
                        model_index = model_index,
                        report_parameters = self.report_parameters,
                        seed = seed,
                        schema = 'abctoolbox',
                        stat_patterns = self.stat_patterns,
                        summary_worker = None,
                        write_header_file = False,
                        tag = model_index)
                prior_workers.append(worker)
            if n >= self.num_processors:
                n = 0
                yield [prior_workers.pop(0) for i in range(len(prior_workers))]
        if self.num_extra_samples_remaining > 0:
            for model_index, config_path in self.models.iteritems():
                prior_path = self.prior_temp_fs.get_file_path(
                        parent = self.prior_temp_fs.base_dir,
                        prefix = 'prior-{0}-{1}-'.format(
                                self.model_strings[model_index],
                                i + 2),
                        register = False)
                seed = functions.get_random_int(self.rng)
                worker = workers.MsBayesWorker(
                        temp_fs = self.prior_temp_fs,
                        sample_size = self.num_extra_samples,
                        config_path = config_path,
                        prior_path = prior_path,
                        exe_path = self.msbayes_exe_path,
                        model_index = model_index,
                        report_parameters = self.report_parameters,
                        seed = seed,
                        schema = 'abctoolbox',
                        stat_patterns = self.stat_patterns,
                        summary_worker = sum_worker,
                        write_header_file = False,
                        tag = model_index)
                prior_workers.append(worker)
        if prior_workers:
            yield prior_workers

    def _parse_observed_paths(self):
        for obs_idx, obs_path in self.observed_stats_paths.iteritems():
            obs_file, close = fileio.process_file_arg(obs_path)
            header = parsing.parse_header(obs_file, seek = False)
            all_stat_indices = parsing.get_stat_indices(header,
                    stat_patterns=parsing.ALL_STAT_PATTERNS)
            for i, line in enumerate(obs_file):
                self.num_observed += 1
                self.n_observed_map[obs_idx] += 1
                tmp_obs_path = self.get_observed_path(obs_idx, i)
                l = line.strip().split()
                with open(tmp_obs_path, 'w') as out:
                    out.write('{0}\n{1}\n'.format(
                        '\t'.join([header[idx] for idx in all_stat_indices]),
                        '\t'.join([str(l[idx]) for idx in all_stat_indices])))
            if close:
                obs_file.close()
        if ((self.num_observed < 2) and
                (self.num_processors > 1) and 
                (not self.use_previous_priors) and
                (self.num_batches > 1)):
            self.duplicate_rejection_workers = True
        if not self.global_estimate_only:
            self.num_samples_processed = dict(zip(self.models.keys(),
                    [self.num_posterior_samples * \
                            self.num_observed for i in range(len(
                                    self.models.keys()))]))
        else:
            self.num_samples_processed = {
                    'combined': self.num_posterior_samples * self.num_observed}
        if self.duplicate_rejection_workers:
            for model_idx in self.num_samples_processed.iterkeys():
                self.num_samples_processed[model_idx] *= self.num_processors

    def _write_keys(self):
        out, close = fileio.process_file_arg(self.model_key_path, 'w')
        for idx, cfg_path in self.models.iteritems():
            out.write('m{0} = {1}\n'.format(idx, os.path.relpath(cfg_path,
                    os.path.dirname(self.model_key_path))))
        out.close()
        out, close = fileio.process_file_arg(self.data_key_path, 'w')
        for idx, data_path in self.observed_stats_paths.iteritems():
            out.write('d{0} = {1}\n'.format(idx, os.path.relpath(data_path,
                    os.path.dirname(self.data_key_path))))
        out.close()

    def _rejection_worker_iter(self, max_num_workers = 500):
        for observed_idx, n_observed in self.n_observed_map.iteritems():
            if observed_idx < self.start_from_observed:
                continue
            for i in range(self.start_from_simulation, n_observed,
                    max_num_workers):
                rejection_workers = {}
                for j in range(i, i + max_num_workers):
                    if j >= n_observed:
                        break
                    obs_path, post_paths = self.get_observed_posterior_paths(
                            observed_idx,
                            j,
                            per_processor = self.duplicate_rejection_workers)
                    combined_post_paths = None
                    if post_paths.has_key('combined'):
                        combined_post_paths = post_paths.pop('combined')
                    if self.global_estimate_only:
                        post_paths = {'combined': combined_post_paths}
                    for model_idx, p_paths in post_paths.iteritems():
                        for k, p_path in enumerate(p_paths):
                            prior_paths = []
                            if os.path.exists(p_path):
                                temp_post_path = self.temp_fs.get_file_path(
                                        parent = self.old_posterior_temp_dir,
                                        prefix = ('{0}-d{1}-{2}-s{3}-{4}-'
                                                  'post-tmp-'.format(
                                                self.temp_fs.token_id,
                                                observed_idx,
                                                self.model_strings[model_idx],
                                                j, k)),
                                        register = False)
                                shutil.copy(p_path, temp_post_path)
                                prior_paths.append(temp_post_path)
                            rw = workers.EuRejectWorker(
                                    temp_fs = self.temp_fs,
                                    observed_path = obs_path,
                                    prior_paths = prior_paths,
                                    num_posterior_samples = \
                                            self.num_posterior_samples,
                                    num_standardizing_samples = 0,
                                    summary_in_path = \
                                            self.summary_paths[model_idx],
                                    summary_out_path = None,
                                    posterior_path = p_path,
                                    regression_worker = None,
                                    exe_path = self.eureject_exe_path,
                                    keep_temps = self.keep_temps,
                                    tag = model_idx)
                            if rejection_workers.has_key(rw.tag):
                                rejection_workers[rw.tag].append(rw)
                            else:
                                rejection_workers[rw.tag] = [rw]
                yield rejection_workers

    def _merging_rejection_worker_iter(self, max_num_workers = 500):
        for observed_idx, n_observed in self.n_observed_map.iteritems():
            if observed_idx < self.start_from_observed:
                continue
            for i in range(self.start_from_simulation, n_observed,
                    max_num_workers):
                rejection_workers = []
                for j in range(i, i + max_num_workers):
                    if j >= n_observed:
                        break
                    obs_path, post_paths = self.get_observed_posterior_paths(
                            observed_idx,
                            j,
                            per_processor = True)
                    final_post_paths = self.get_posterior_paths(
                            observed_idx,
                            j,
                            per_processor = False)
                    combined_post_paths = None
                    if post_paths.has_key('combined'):
                        combined_post_paths = post_paths.pop('combined')
                    if self.global_estimate_only:
                        post_paths = {'combined': combined_post_paths}
                    for model_idx, p_paths in post_paths.iteritems():
                        prior_paths = []
                        for k, p_path in enumerate(p_paths):
                            if os.path.exists(p_path):
                                prior_paths.append(p_path)
                        if prior_paths:
                            rejection_workers.append(workers.EuRejectWorker(
                                    temp_fs = self.temp_fs,
                                    observed_path = obs_path,
                                    prior_paths = prior_paths,
                                    num_posterior_samples = \
                                            self.num_posterior_samples,
                                    num_standardizing_samples = 0,
                                    summary_in_path = \
                                            self.summary_paths[model_idx],
                                    summary_out_path = None,
                                    posterior_path = final_post_paths[model_idx][0],
                                    regression_worker = None,
                                    exe_path = self.eureject_exe_path,
                                    keep_temps = self.keep_temps,
                                    tag = model_idx))
                if rejection_workers:
                    yield rejection_workers

    def _final_rejection_worker_iter(self, max_num_workers = 500):
        for observed_idx, n_observed in self.n_observed_map.iteritems():
            if observed_idx < self.start_from_observed:
                continue
            for i in range(self.start_from_simulation, n_observed,
                    max_num_workers):
                rejection_workers = []
                for j in range(i, i + max_num_workers):
                    if j >= n_observed:
                        break
                    obs_path, post_paths = self.get_observed_posterior_paths(
                            observed_idx,
                            j,
                            per_processor = False)
                    combined_post_paths = post_paths.pop('combined')
                    prior_paths = []
                    for model_idx, p_paths in post_paths.iteritems():
                        assert os.path.exists(p_paths[0])
                        prior_paths.append(p_paths[0])
                    if prior_paths:
                        rejection_workers.append(workers.EuRejectWorker(
                                temp_fs = self.temp_fs,
                                observed_path = obs_path,
                                prior_paths = prior_paths,
                                num_posterior_samples = \
                                        self.num_posterior_samples,
                                num_standardizing_samples = 0,
                                summary_in_path = \
                                        self.summary_paths['combined'],
                                summary_out_path = None,
                                posterior_path = combined_post_paths[0],
                                regression_worker = None,
                                exe_path = self.eureject_exe_path,
                                keep_temps = self.keep_temps,
                                tag = model_idx))
                if rejection_workers:
                    yield rejection_workers

    def _regression_worker_iter(self, batch_index, max_num_workers = 500):
        for observed_idx, n_observed in self.n_observed_map.iteritems():
            if observed_idx < self.start_from_observed:
                continue
            for i in range(self.start_from_simulation, n_observed,
                    max_num_workers):
                posterior_workers = []
                for j in range(i, i + max_num_workers):
                    if j >= n_observed:
                        break
                    obs_path, post_paths = self.get_observed_posterior_paths(
                            observed_idx,
                            j,
                            per_processor = False)
                    if self.global_estimate_only:
                        post_paths = {'combined': post_paths['combined']}
                    if (not self.global_estimate) and post_paths.has_key(
                            'combined'):
                        post_paths.pop('combined')
                    for model_idx, p_paths in post_paths.iteritems():
                        fname_parts = os.path.basename(p_paths[0]).split('-')
                        if model_idx == 'combined':
                            fname_prefix = self.output_prefix + \
                                    '-'.join(fname_parts[1:5] + \
                                            [str(batch_index)])
                        else:
                            fname_prefix = self.output_prefix + \
                                    '-'.join(fname_parts[1:4] + \
                                            [str(batch_index)])
                        out_prefix = os.path.join(
                                self.model_dirs[observed_idx][model_idx],
                                fname_prefix)
                        post_out_path = out_prefix + '-posterior-sample.txt'
                        if model_idx == 'combined':
                            model_indices = self.models.keys()
                        else:
                            model_indices = [model_idx]

                        posterior_workers.append(workers.PosteriorWorker(
                                temp_fs = self.temp_fs,
                                observed_path = obs_path,
                                posterior_path = p_paths[0],
                                num_taxon_pairs = self.num_taxon_pairs,
                                posterior_out_path = post_out_path,
                                output_prefix = out_prefix,
                                model_indices = model_indices,
                                abctoolbox_exe_path = self.abctoolbox_exe_path,
                                abctoolbox_num_posterior_quantiles = \
                                        self.num_posterior_density_quantiles,
                                omega_threshold = self.omega_threshold,
                                compress = self.compress,
                                keep_temps = self.keep_temps,
                                tag = model_idx))
                if posterior_workers:
                    yield posterior_workers
    
    def _run_rejection_on_provided_priors(self, remove_files = True):
        for k in self.num_samples_processed.iterkeys():
            self.num_samples_processed[k] = 0
        for observed_idx, n_observed in self.n_observed_map.iteritems():
            if observed_idx < self.start_from_observed:
                continue
            for i in range(self.start_from_simulation, n_observed,
                    self.num_processors):
                rejection_workers = []
                global_rejection_workers = []
                for j in range(i, i + self.num_processors):
                    if j >= n_observed:
                        break
                    obs_path, post_paths = self.get_observed_posterior_paths(
                            observed_idx,
                            j,
                            per_processor = False)
                    combined_post_paths = None
                    if post_paths.has_key('combined'):
                        combined_post_paths = post_paths.pop('combined')
                        assert len(combined_post_paths) == 1
                    if self.global_estimate_only:
                        post_paths = {'combined': combined_post_paths}
                        combined_post_paths = None
                    for model_idx, p_paths in post_paths.iteritems():
                        assert len(p_paths) == 1
                        p_path = p_paths[0]
                        fname_parts = os.path.basename(p_path).split('-')
                        if model_idx == 'combined':
                            fname_prefix = self.output_prefix + \
                                    '-'.join(fname_parts[1:5] + \
                                            [str(self.num_prior_batch_iters)])
                        else:
                            fname_prefix = self.output_prefix + \
                                    '-'.join(fname_parts[1:4] + \
                                            [str(self.num_prior_batch_iters)])
                        out_prefix = os.path.join(
                                self.model_dirs[observed_idx][model_idx],
                                fname_prefix)
                        post_out_path = out_prefix + '-posterior-sample.txt'
                        if model_idx == 'combined':
                            model_indices = self.models.keys()
                        else:
                            model_indices = [model_idx]
                        pw = workers.PosteriorWorker(
                                temp_fs = self.temp_fs,
                                observed_path = obs_path,
                                posterior_path = p_path,
                                num_taxon_pairs = self.num_taxon_pairs,
                                posterior_out_path = post_out_path,
                                output_prefix = out_prefix,
                                model_indices = model_indices,
                                abctoolbox_exe_path = self.abctoolbox_exe_path,
                                abctoolbox_num_posterior_quantiles = \
                                        self.num_posterior_density_quantiles,
                                omega_threshold = self.omega_threshold,
                                compress = self.compress,
                                keep_temps = self.keep_temps,
                                tag = model_idx)
                        rw = workers.EuRejectWorker(
                                temp_fs = self.temp_fs,
                                observed_path = obs_path,
                                prior_paths = [self.model_samples[model_idx]],
                                num_posterior_samples = \
                                        self.num_posterior_samples,
                                num_standardizing_samples = 0,
                                summary_in_path = \
                                        self.summary_paths[model_idx],
                                summary_out_path = None,
                                posterior_path = p_path,
                                regression_worker = pw,
                                exe_path = self.eureject_exe_path,
                                keep_temps = self.keep_temps,
                                tag = model_idx)
                        rejection_workers.append(rw)
                    if combined_post_paths and self.global_estimate:
                        p_path = combined_post_paths[0]
                        model_idx = 'combined'
                        fname_parts = os.path.basename(p_path).split('-')
                        fname_prefix = self.output_prefix + \
                                '-'.join(fname_parts[1:5] + \
                                        [str(self.num_prior_batch_iters)])
                        out_prefix = os.path.join(
                                self.model_dirs[observed_idx][model_idx],
                                fname_prefix)
                        post_out_path = out_prefix + '-posterior-sample.txt'
                        model_indices = self.models.keys()
                        prior_paths = [pp[0] for pp in post_paths.itervalues()]
                        pw = workers.PosteriorWorker(
                                temp_fs = self.temp_fs,
                                observed_path = obs_path,
                                posterior_path = p_path,
                                num_taxon_pairs = self.num_taxon_pairs,
                                posterior_out_path = post_out_path,
                                output_prefix = out_prefix,
                                model_indices = model_indices,
                                abctoolbox_exe_path = self.abctoolbox_exe_path,
                                abctoolbox_num_posterior_quantiles = \
                                        self.num_posterior_density_quantiles,
                                omega_threshold = self.omega_threshold,
                                compress = self.compress,
                                keep_temps = self.keep_temps,
                                tag = model_idx)
                        rw = workers.EuRejectWorker(
                                temp_fs = self.temp_fs,
                                observed_path = obs_path,
                                prior_paths = prior_paths,
                                num_posterior_samples = \
                                        self.num_posterior_samples,
                                num_standardizing_samples = 0,
                                summary_in_path = \
                                        self.summary_paths[model_idx],
                                summary_out_path = None,
                                posterior_path = p_path,
                                regression_worker = pw,
                                exe_path = self.eureject_exe_path,
                                keep_temps = self.keep_temps,
                                tag = model_idx)
                        global_rejection_workers.append(rw)
                if rejection_workers:
                    rejection_workers = self._run_workers(rejection_workers)
                    for rw in rejection_workers:
                        self.num_samples_processed[rw.tag] += \
                                rw.num_processed
                if global_rejection_workers:
                    global_rejection_workers = self._run_workers(
                            global_rejection_workers)
                if remove_files:
                    for rw in rejection_workers:
                        self.temp_fs.remove_file(rw.posterior_path)
                    for rw in global_rejection_workers:
                        self.temp_fs.remove_file(rw.posterior_path)

    def _run_workers(self, worker_iter, queue_max = 500):
        return manager.Manager.run_workers(
            workers = worker_iter,
            num_processors = self.num_processors,
            queue_max = queue_max)

    def _run_prior_summary_workers(self):
        sum_workers = self._run_workers(list(self.prior_summary_worker_iter))
        prior_paths = {}
        summary_workers = []
        for w in sum_workers:
            self.num_samples_generated += w.sample_size
            self.num_samples_summarized += w.summary_worker.num_summarized
            if prior_paths.has_key(w.tag):
                prior_paths[w.tag].append(str(w.prior_path))
            else:
                prior_paths[w.tag] = [str(w.prior_path)]
            summary_workers.append(copy.deepcopy(w.summary_worker))
            w.purge()
        return prior_paths, summary_workers

    def _run_prior_workers(self, prior_worker_batch):
        pworkers = self._run_workers(prior_worker_batch)
        prior_paths = {}
        for w in pworkers:
            self.num_samples_generated += w.sample_size
            assert w.summary_worker is None
            if prior_paths.has_key(w.tag):
                prior_paths[w.tag].append(str(w.prior_path))
            else:
                prior_paths[w.tag] = [str(w.prior_path)]
            w.purge()
        return prior_paths

    def _run_rejection_workers(self, prior_paths, remove_files = True):
        if self.global_estimate_only:
            p_paths = []
            for path_list in prior_paths.itervalues():
                p_paths.extend(path_list)
            prior_paths = {'combined': p_paths}
        for i, rej_worker_batch in enumerate(self._rejection_worker_iter()):
            if DUMP_DEBUG_INFO:
                with open('pymsbayes-debug.log', 'a') as out:
                    out.write('########################################')
                    out.write('########################################\n')
                    out.write('ABCTeam._run_rejection_workers iter {0}:\n'.format(i))
                    out.write('MEMORY DUMP:\n')
                    dump_memory_info(out)
                    out.write('OPEN FILES ({0}):\n'.format(
                            len(fileio.FileStream.open_files)))
                    out.write('{0}\n'.format('\n'.join(
                            fileio.FileStream.open_files)))
                    out.write('REGISTERED TEMP DIRS ({0}):\n'.format(len(self.temp_fs.dirs)))
                    out.write('{0}\n'.format('\n'.join(self.temp_fs.dirs)))
                    out.write('REGISTERED TEMP FILES ({0}):\n'.format(len(self.temp_fs.files)))
                    out.write('{0}\n'.format('\n'.join(self.temp_fs.files)))
                    out.write('REGISTERED PRIOR TEMP DIRS ({0}):\n'.format(len(self.prior_temp_fs.dirs)))
                    out.write('{0}\n'.format('\n'.join(self.prior_temp_fs.dirs)))
                    out.write('REGISTERED PRIOR TEMP FILES ({0}):\n'.format(len(self.prior_temp_fs.files)))
                    out.write('{0}\n'.format('\n'.join(self.prior_temp_fs.files)))
                    out.write('\n')
            to_run = []
            for model_idx, path_list in prior_paths.iteritems():
                if not self.duplicate_rejection_workers:
                    for rej_worker in rej_worker_batch[model_idx]:
                        rej_worker.prior_paths.extend(path_list)
                        to_run.append(rej_worker)
                else:
                    path_sub_lists = list(functions.list_splitter(path_list,
                            len(rej_worker_batch[model_idx])))
                    for j in range(len(path_sub_lists)):
                        rej_worker_batch[model_idx][j].prior_paths.extend(path_sub_lists[j])
                        to_run.append(rej_worker_batch[model_idx][j])
            to_run = self._run_workers(to_run)
            for rej_worker in to_run:
                self.num_samples_processed[rej_worker.tag] += \
                        (rej_worker.num_processed - self.num_posterior_samples)
            self._purge_old_posterior_temp_dir()
        if remove_files:
            for path_list in prior_paths.itervalues():
                for p in path_list:
                    self.prior_temp_fs.remove_file(p)

    def _purge_old_posterior_temp_dir(self):
        self.temp_fs.clear_dir(self.old_posterior_temp_dir)
        assert os.listdir(self.old_posterior_temp_dir) == []

    def _run_merging_rejection_workers(self, remove_files = False):
        if not self.duplicate_rejection_workers:
            return
        for i, rej_worker_batch in enumerate(self._merging_rejection_worker_iter()):
            rej_worker_batch = self._run_workers(rej_worker_batch)
            if remove_files:
                for rej_worker in rej_worker_batch:
                    for p in rej_worker.prior_paths:
                        self.temp_fs.remove_file(p)
    
    def _run_final_rejections(self):
        if not self.global_estimate:
            return
        if self.global_estimate_only:
            return
        if len(self.models) < 2:
            return
        for i, rej_worker_batch in enumerate(self._final_rejection_worker_iter()):
            rej_worker_batch = self._run_workers(rej_worker_batch)

    def _run_regression_workers(self, batch_index, remove_files = False):
        for i, reg_worker_batch in enumerate(self._regression_worker_iter(
                batch_index)):
            reg_worker_batch = self._run_workers(reg_worker_batch)
            if remove_files:
                for reg_worker in reg_worker_batch:
                    self.temp_fs.remove_file(reg_worker.posterior_path)

    def _merge_summaries(self, summary_workers):
        sum_workers = dict(zip(self.models.keys(),
                [[] for i in range(len(self.models.keys()))]))
        sum_mergers = []
        summaries = {}
        for sw in summary_workers:
            sum_workers[sw.tag].append(sw)
        for model_idx, s_workers in sum_workers.iteritems():
            summary_merger = workers.EuRejectSummaryMerger(s_workers,
                    tag = model_idx)
            sum_mergers.append(summary_merger)
        sum_mergers = self._run_workers(sum_mergers)
        for sm in sum_mergers:
            summaries[sm.tag] = sm.sample_sum_collection
        if len(self.models) > 1:
            summaries['combined'] = stats.SampleSummaryCollection.merge(
                    summaries.itervalues())
        for k, summary in summaries.iteritems():
            summary.write(self.summary_paths[k])

    def _write_prior_samples(self, prior_path_dict, remove_files = True):
        compresslevel = None
        if self.compress:
            compresslevel = 9
        for model_idx, prior_paths in prior_path_dict.iteritems():
            workers.merge_prior_files(
                    paths = prior_paths,
                    dest_path = self._get_prior_path(model_idx),
                    append = True,
                    compresslevel = compresslevel)
        if remove_files:
            for path_list in prior_path_dict.itervalues():
                for p in path_list:
                    self.prior_temp_fs.remove_file(p)

    def run(self):
        self._parse_observed_paths()
        if self.generate_prior_samples_only:
            _LOG.info('Generating prior samples only...')
            self._generate_prior_samples()
        else:
            _LOG.info('Running full analysis...')
            self._run_full_analysis()

    def _process_prior_summary_workers(self,
            run_rejection = True,
            write_priors = False):
        _LOG.info('Running prior summary workers...')
        prior_paths, summary_workers = self._run_prior_summary_workers()

        _LOG.info('Merging summaries...')
        self._merge_summaries(summary_workers)

        if run_rejection:
            _LOG.info('Running rejection on summary worker priors...')
            self._run_rejection_workers(prior_paths)
        if write_priors:
            _LOG.info('Writing prior samples...')
            self._write_prior_samples(prior_paths)

        _LOG.info('Number of samples generated: {0}'.format(
            self.num_samples_generated))
        _LOG.info('Number of samples summarized: {0}'.format(
            self.num_samples_summarized))

    def _generate_prior_samples(self):
        self._process_prior_summary_workers(run_rejection = False,
                write_priors = True)

        for i, prior_worker_batch in enumerate(self.prior_worker_iter):
            _LOG.info('Running prior worker batch {0} of {1}...'.format(
                    (i + 1), self.num_prior_batch_iters))
            prior_paths = self._run_prior_workers(prior_worker_batch)

            _LOG.info('Writing prior samples...')
            self._write_prior_samples(prior_paths)
            self.iter_count += 1

        self.finished = True

    def _run_full_analysis(self):
        if self.use_previous_priors:
            _LOG.info('Running full analyses on provided priors...')
            self._run_rejection_on_provided_priors()
            self.iter_count += 1
            self.finished = True
            return

        else:
            self._process_prior_summary_workers(run_rejection = True,
                    write_priors = False)

            for i, prior_worker_batch in enumerate(self.prior_worker_iter):
                if self.reporting_indices and (i == self.reporting_indices[0]):
                    self.reporting_indices.pop(0)
                    _LOG.info('Reporting results at iteration {0} of {1}'
                            '...'.format((i + 1), self.num_prior_batch_iters))
                    self._run_merging_rejection_workers(remove_files = False)
                    self._run_final_rejections()
                    self._run_regression_workers(i+1, remove_files = False)

                _LOG.info('Running prior worker batch {0} of {1}...'.format(
                        (i + 1), self.num_prior_batch_iters))
                prior_paths = self._run_prior_workers(prior_worker_batch)

                _LOG.info('Running rejection on prior batch {0} of {1}'
                        '...'.format((i + 1), self.num_prior_batch_iters))
                self._run_rejection_workers(prior_paths)
                self.iter_count += 1

        _LOG.info('Merging rejection team posteriors...')
        self._run_merging_rejection_workers(remove_files = True)

        _LOG.info('Running final rejections...')
        self._run_final_rejections()

        _LOG.info('Number of samples processed:')
        for k, n in self.num_samples_processed.iteritems():
            _LOG.info('    {0}: {1}'.format(k, n))

        _LOG.info('Running final regressions...')
        self._run_regression_workers(self.num_prior_batch_iters,
                remove_files = True)
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
        self.observed_path = fileio.expand_path(observed_path)
        self.output_dir = fileio.expand_path(output_dir)
        self.num_taxon_pairs = int(num_taxon_pairs)
        if model_indices:
            model_indices = list(set(model_indices))
        self.model_indices = model_indices
        self.num_posterior_samples = int(num_posterior_samples)
        self.num_posterior_density_quantiles = int(
                num_posterior_density_quantiles)
        self.eureject_exe_path = eureject_exe_path
        self.abctoolbox_exe_path = abctoolbox_exe_path
        self.abctoolbox_bandwidth = abctoolbox_bandwidth
        self.omega_threshold = omega_threshold
        self.compress = compress
        self.prior_paths = [p for p in prior_paths]
        if summary_in_path:
            summary_in_path = fileio.expand_path(summary_in_path)
        self.summary_in_path = summary_in_path
        self.keep_temps = bool(keep_temps)
        self.keep_priors = True
        self.posterior_path = None
        self.run_regression = run_regression
        self.final_regression = False
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
        self.num_samples_processed = 0
        self.posterior_size = 0
        self.num_rejection_runs = 0
        self.num_regression_runs = 0

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
        self.num_rejection_runs += 1
        p_paths = []
        num_added = 0
        old_posterior_path = None
        prior_paths = copy.deepcopy(self.prior_paths)
        for i in range(len(self.prior_paths)):
            path = self.prior_paths.pop(0)
            p_paths.append(path)
        if self.posterior_path:
            old_posterior_path = str(self.posterior_path)
            p_paths.append(old_posterior_path)
            num_added = self.posterior_size
        self.posterior_path = self.temp_fs.get_file_path(
                prefix = 'posterior-{0}-{1}-{2}-'.format(self.name, 
                        self.num_rejection_runs, self.index))
        rw = workers.EuRejectWorker(
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
        self.posterior_size = rw.num_retained
        if rw.num_processed:
            self.num_samples_processed = self.num_samples_processed + \
                    rw.num_processed - num_added
        if not self.keep_temps:
            if old_posterior_path:
                os.remove(old_posterior_path)
        if not self.keep_priors:
            for pp in prior_paths:
                os.remove(pp)
        assert len(self.prior_paths) == 0

    def _run_regression_workers(self):
        self.num_regression_runs += 1
        out_prefix = self.output_prefix + '-' + str(self.num_rejection_runs)
        post_path = out_prefix + '-posterior-sample.txt'
        pw = workers.PosteriorWorker(
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
                        self.num_posterior_density_quantiles,
                omega_threshold = self.omega_threshold,
                compress = self.compress,
                keep_temps = self.keep_temps,
                tag = self.tag)
        pw.start()
        if self.final_regression:
            os.remove(self.posterior_path)
        self.posterior_path = post_path
        self.div_model_results_path = pw.div_model_results_path
        self.model_results_path = pw.model_results_path
        self.psi_results_path = pw.psi_results_path
        self.omega_results_path = pw.omega_results_path
        self.posterior_summary_path = pw.posterior_summary_path
        self.regress_summary_path = pw.regress_summary_path
        self.regress_posterior_path = pw.regress_posterior_path

class ModelProbabilityEstimatorTeam(object):
    count = 0
    def __init__(self,
            config_paths,
            num_samples = 1000,
            omega_threshold = 0.01,
            num_processors = 4,
            rng = None):
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.rng = rng
        if not rng:
            self.rng = GLOBAL_RNG
        self.np = num_processors
        self.omega_threshold = omega_threshold
        self.configs = dict(zip(config_paths,
            [config.MsBayesConfig(c) for c in config_paths]))
        self.psi_summaries = {}
        self.psi_probs = {}
        self.shared_div_summaries = {}
        self.shared_div_probs = {}
        for p, c in self.configs.iteritems():
            self.psi_summaries[p] = dict(zip([i + 1 for i in range(c.npairs)],
                    [stats.SampleSummary() for i in range(c.npairs)]))
            self.psi_probs[p] = dict(zip([i + 1 for i in range(c.npairs)],
                    [None for i in range(c.npairs)]))
            self.shared_div_summaries[p] = dict(zip(
                    [i + 1 for i in range(1, c.npairs)],
                    [stats.SampleSummary() for i in range(1, c.npairs)]))
            self.shared_div_probs[p] = dict(zip(
                    [i + 1 for i in range(1, c.npairs)],
                    [None for i in range(1, c.npairs)]))
        self.omega_summaries = dict(zip(self.configs.iterkeys(),
                [stats.SampleSummary() for c in self.configs.iterkeys()]))
        self.omega_probs = dict(zip(self.configs.iterkeys(),
                [None for c in self.configs.iterkeys()]))
        self.num_samples = num_samples
        if self.np > self.num_samples:
            self.np = self.num_samples
        self.batch_size, self.remainder = functions.long_division(
                self.num_samples,
                self.np)

    def start(self):
        eworkers = []
        for path, cfg in self.configs.iteritems():
            for i in range(self.np):
                sample_size = self.batch_size
                if i == (self.np - 1):
                    sample_size += self.remainder
                w = workers.ModelProbabilityEstimator(
                        config = cfg,
                        num_samples = sample_size,
                        omega_threshold = self.omega_threshold,
                        tag = path)
                eworkers.append(w)
        _LOG.info('{0}: Generating samples...'.format(self.name))
        eworkers = manager.Manager.run_workers(
                workers = eworkers,
                num_processors = self.np)
        _LOG.info('{0}: Done!'.format(self.name))
        _LOG.info('{0}: Summarizing results...'.format(self.name))
        for w in eworkers:
            for k in self.psi_summaries[w.tag].iterkeys():
                self.psi_summaries[w.tag][k].update(w.psi_summary[k])
            for k in self.shared_div_summaries[w.tag].iterkeys():
                self.shared_div_summaries[w.tag][k].update(
                        w.shared_div_summary[k])
            self.omega_summaries[w.tag].update(w.omega_summary)
        for path in self.configs.iterkeys():
            total = 0.0
            for k, s in self.psi_summaries[path].iteritems():
                self.psi_probs[path][k] = s.mean
                total += s.mean
            if not probability.almost_equal(total, 1.0, places=7):
                raise Exception('Error in estimating probabilities of the '
                        'number of divergences for {0}. The total probability '
                        'summed to {1}'.format(path, total))
            for k, s in self.shared_div_summaries[path].iteritems():
                self.shared_div_probs[path][k] = s.mean
            if not probability.almost_equal(
                    self.shared_div_probs[path][self.configs[path].npairs],
                    self.psi_probs[path][1],
                    places = 7):
                raise Exception('Error in estimating probabilities of the '
                        'shared divergences for {0}. The probability of all'
                        'taxa co-diverging was {1}, but the probability of '
                        'one divergence event was {2}'.format(path,
                                self.shared_div_probs[path][
                                        self.configs[path].npairs],
                                self.psi_probs[path][1]))
            self.omega_probs[path] = self.omega_summaries[path].mean

class DivModelSimulatorTeam(object):
    count = 0
    def __init__(self,
            config_paths,
            num_samples = 1000,
            num_processors = 4,
            rng = None):
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.rng = rng
        if not rng:
            self.rng = GLOBAL_RNG
        self.np = num_processors
        self.configs = dict(zip(config_paths,
            [config.MsBayesConfig(c) for c in config_paths]))
        self.div_models = dict(zip(config_paths,
            [stats.PartitionCollection() for c in config_paths]))
        self.num_samples = num_samples
        if self.np > self.num_samples:
            self.np = self.num_samples
        self.batch_size, self.remainder = functions.long_division(
                self.num_samples,
                self.np)

    def start(self):
        sworkers = []
        for path, cfg in self.configs.iteritems():
            for i in range(self.np):
                sample_size = self.batch_size
                if i == (self.np - 1):
                    sample_size += self.remainder
                w = workers.DivModelSimulator(
                        config = cfg,
                        num_samples = sample_size,
                        tag = path)
                sworkers.append(w)
        _LOG.info('{0}: Generating samples...'.format(self.name))
        sworkers = manager.Manager.run_workers(
                workers = sworkers,
                num_processors = self.np)
        _LOG.info('{0}: Done!'.format(self.name))
        _LOG.info('{0}: Summarizing results...'.format(self.name))
        for w in sworkers:
            self.div_models[w.tag].add_iter(
                    w.div_models.partitions.itervalues())
