#! /usr/bin/env python

import unittest
import os
import sys
import random
import re
import glob

from pymsbayes.teams import *
from pymsbayes.workers import (MsBayesWorker, ABCToolBoxRejectWorker,
        MsRejectWorker, merge_priors, merge_prior_files)
from pymsbayes.test.support import package_paths
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.test import TestLevel, test_enabled
from pymsbayes.utils import GLOBAL_RNG
from pymsbayes.utils.parsing import parse_parameters, add_div_model_column
from pymsbayes.utils.stats import IntegerPartitionCollection
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class RejectionTeamTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')

    def tearDown(self):
        self.tear_down()

    def test_one_prior_path(self):
        obs_worker = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        obs_worker.start()
        prior_worker = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 100,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        prior_worker.start()

        post_path = self.get_test_path(prefix='posterior')
        summary_worker = EuRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_paths = [prior_worker.prior_path],
                num_posterior_samples = 50,
                num_standardizing_samples = 100,
                summary_in_path = None,
                summary_out_path = None,
                posterior_path = post_path,
                regression_worker = None,
                stderr_path = None,
                keep_temps = True,
                tag = 'test')
        summary_worker.start()

        output_dir = self.get_test_subdir(prefix = 'rejection-team-test')
        rt = RejectionTeam(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                output_dir = output_dir,
                output_prefix = self.test_id,
                num_taxon_pairs = 4,
                prior_paths = [prior_worker.prior_path],
                model_indices = None,
                summary_in_path = summary_worker.summary_out_path,
                num_posterior_samples = 50,
                num_posterior_density_quantiles = 100,
                run_regression = False,
                abctoolbox_bandwidth = None,
                omega_threshold = 0.01,
                compress = False,
                keep_temps = False,
                index = 1,
                tag = 'test')

        self.assertEqual(rt.posterior_path, None)
        rt.start()
        self.assertTrue(os.path.exists(rt.posterior_path))
        self.assertEqual(self.get_number_of_lines(rt.posterior_path), 51)
        self.assertEqual(self.get_number_of_header_lines(rt.posterior_path), 1)
        self.assertEqual(rt.num_samples_processed, 100)

        self.assertSameUnsortedFiles([rt.posterior_path,
                summary_worker.posterior_path])
        
        rt.run_regression = True
        self.assertEqual(rt.div_model_results_path, None)
        self.assertEqual(rt.psi_results_path, None)
        self.assertEqual(rt.omega_results_path, None)
        self.assertEqual(rt.posterior_summary_path, None)
        self.assertEqual(rt.regress_summary_path, None)
        self.assertEqual(rt.regress_posterior_path, None)
        self.assertEqual(rt.model_results_path, None)
        rt.start()
        self.assertTrue(os.path.isfile(rt.div_model_results_path))
        self.assertTrue(os.path.isfile(rt.psi_results_path))
        self.assertTrue(os.path.isfile(rt.omega_results_path))
        self.assertTrue(os.path.isfile(rt.posterior_summary_path))
        self.assertTrue(os.path.isfile(rt.regress_summary_path))
        self.assertTrue(os.path.isfile(rt.regress_posterior_path))
        self.assertFalse(os.path.isfile(rt.model_results_path))

        exp_reg_summary_path = self.get_test_path(prefix='expected-summary')
        exp_reg_post_path = self.get_test_path(prefix='expected-posterior')
        self.assertTrue(os.path.isfile(obs_worker.prior_stats_path))
        self.assertTrue(os.path.isfile(rt.posterior_path))
        reg_worker = ABCToolBoxRegressWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                posterior_path = rt.posterior_path,
                regress_summary_path = exp_reg_summary_path,
                regress_posterior_path = exp_reg_post_path,
                bandwidth = None,
                num_posterior_samples = 50,
                num_posterior_quantiles = 100,
                compress = False,
                tag = 'test')
        reg_worker.start()
        self.assertSameFiles([rt.regress_summary_path, exp_reg_summary_path])
        self.assertSameFiles([rt.regress_posterior_path, exp_reg_post_path])

    def test_mult_prior_paths(self):
        obs_worker = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        obs_worker.start()
        prior_workers = []
        for i in range(2):
            prior_worker = MsBayesWorker(
                    temp_fs = self.temp_fs,
                    sample_size = 50,
                    config_path = self.cfg_path,
                    schema = 'abctoolbox',
                    write_stats_file = True)
            prior_worker.start()
            prior_workers.append(prior_worker)

        post_path = self.get_test_path(prefix='posterior')
        summary_worker = EuRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_paths = [pw.prior_path for pw in prior_workers],
                num_posterior_samples = 50,
                num_standardizing_samples = 100,
                summary_in_path = None,
                summary_out_path = None,
                posterior_path = post_path,
                regression_worker = None,
                stderr_path = None,
                keep_temps = True,
                tag = 'test')
        summary_worker.start()

        output_dir = self.get_test_subdir(prefix = 'rejection-team-test')
        rt = RejectionTeam(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                output_dir = output_dir,
                output_prefix = self.test_id,
                num_taxon_pairs = 4,
                prior_paths = [pw.prior_path for pw in prior_workers],
                model_indices = None,
                summary_in_path = summary_worker.summary_out_path,
                num_posterior_samples = 50,
                num_posterior_density_quantiles = 100,
                run_regression = False,
                abctoolbox_bandwidth = None,
                omega_threshold = 0.01,
                compress = False,
                keep_temps = False,
                index = 1,
                tag = 'test')

        self.assertEqual(rt.posterior_path, None)
        rt.start()
        self.assertTrue(os.path.exists(rt.posterior_path))
        self.assertEqual(self.get_number_of_lines(rt.posterior_path), 51)
        self.assertEqual(self.get_number_of_header_lines(rt.posterior_path), 1)
        self.assertEqual(rt.num_samples_processed, 100)

        self.assertSameUnsortedFiles([rt.posterior_path,
                summary_worker.posterior_path])
        
        rt.run_regression = True
        self.assertEqual(rt.div_model_results_path, None)
        self.assertEqual(rt.psi_results_path, None)
        self.assertEqual(rt.omega_results_path, None)
        self.assertEqual(rt.posterior_summary_path, None)
        self.assertEqual(rt.regress_summary_path, None)
        self.assertEqual(rt.regress_posterior_path, None)
        self.assertEqual(rt.model_results_path, None)
        rt.start()
        self.assertTrue(os.path.isfile(rt.div_model_results_path))
        self.assertTrue(os.path.isfile(rt.psi_results_path))
        self.assertTrue(os.path.isfile(rt.omega_results_path))
        self.assertTrue(os.path.isfile(rt.posterior_summary_path))
        self.assertTrue(os.path.isfile(rt.regress_summary_path))
        self.assertTrue(os.path.isfile(rt.regress_posterior_path))
        self.assertFalse(os.path.isfile(rt.model_results_path))

        exp_reg_summary_path = self.get_test_path(prefix='expected-summary')
        exp_reg_post_path = self.get_test_path(prefix='expected-posterior')
        self.assertTrue(os.path.isfile(obs_worker.prior_stats_path))
        self.assertTrue(os.path.isfile(rt.posterior_path))
        reg_worker = ABCToolBoxRegressWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                posterior_path = rt.posterior_path,
                regress_summary_path = exp_reg_summary_path,
                regress_posterior_path = exp_reg_post_path,
                bandwidth = None,
                num_posterior_samples = 50,
                num_posterior_quantiles = 100,
                compress = False,
                tag = 'test')
        reg_worker.start()
        self.assertSameFiles([rt.regress_summary_path, exp_reg_summary_path])
        self.assertSameFiles([rt.regress_posterior_path, exp_reg_post_path])

    def test_mult_prior_paths_and_models(self):
        obs_worker = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        obs_worker.start()
        prior_workers = []
        for i in range(2):
            prior_worker = MsBayesWorker(
                    temp_fs = self.temp_fs,
                    sample_size = 50,
                    config_path = self.cfg_path,
                    schema = 'abctoolbox',
                    model_index = i+1,
                    write_stats_file = True)
            prior_worker.start()
            prior_workers.append(prior_worker)

        post_path = self.get_test_path(prefix='posterior')
        summary_worker = EuRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_paths = [pw.prior_path for pw in prior_workers],
                num_posterior_samples = 50,
                num_standardizing_samples = 100,
                summary_in_path = None,
                summary_out_path = None,
                posterior_path = post_path,
                regression_worker = None,
                stderr_path = None,
                keep_temps = True,
                tag = 'test')
        summary_worker.start()

        output_dir = self.get_test_subdir(prefix = 'rejection-team-test')
        rt = RejectionTeam(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                output_dir = output_dir,
                output_prefix = self.test_id,
                num_taxon_pairs = 4,
                prior_paths = [pw.prior_path for pw in prior_workers],
                model_indices = None,
                summary_in_path = summary_worker.summary_out_path,
                num_posterior_samples = 50,
                num_posterior_density_quantiles = 100,
                run_regression = False,
                abctoolbox_bandwidth = None,
                omega_threshold = 0.01,
                compress = False,
                keep_temps = False,
                index = 1,
                tag = 'test')

        self.assertEqual(rt.posterior_path, None)
        rt.start()
        self.assertTrue(os.path.exists(rt.posterior_path))
        self.assertEqual(self.get_number_of_lines(rt.posterior_path), 51)
        self.assertEqual(self.get_number_of_header_lines(rt.posterior_path), 1)
        self.assertEqual(rt.num_samples_processed, 100)

        self.assertSameUnsortedFiles([rt.posterior_path,
                summary_worker.posterior_path])
        
        rt.run_regression = True
        self.assertEqual(rt.div_model_results_path, None)
        self.assertEqual(rt.psi_results_path, None)
        self.assertEqual(rt.omega_results_path, None)
        self.assertEqual(rt.posterior_summary_path, None)
        self.assertEqual(rt.regress_summary_path, None)
        self.assertEqual(rt.regress_posterior_path, None)
        self.assertEqual(rt.model_results_path, None)
        rt.start()
        self.assertTrue(os.path.isfile(rt.div_model_results_path))
        self.assertTrue(os.path.isfile(rt.psi_results_path))
        self.assertTrue(os.path.isfile(rt.omega_results_path))
        self.assertTrue(os.path.isfile(rt.posterior_summary_path))
        self.assertTrue(os.path.isfile(rt.regress_summary_path))
        self.assertTrue(os.path.isfile(rt.regress_posterior_path))
        self.assertTrue(os.path.isfile(rt.model_results_path))

        exp_reg_summary_path = self.get_test_path(prefix='expected-summary')
        exp_reg_post_path = self.get_test_path(prefix='expected-posterior')
        self.assertTrue(os.path.isfile(obs_worker.prior_stats_path))
        self.assertTrue(os.path.isfile(rt.posterior_path))
        reg_worker = ABCToolBoxRegressWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                posterior_path = rt.posterior_path,
                regress_summary_path = exp_reg_summary_path,
                regress_posterior_path = exp_reg_post_path,
                bandwidth = None,
                num_posterior_samples = 50,
                num_posterior_quantiles = 100,
                compress = False,
                tag = 'test')
        reg_worker.start()
        self.assertSameFiles([rt.regress_summary_path, exp_reg_summary_path])
        self.assertSameFiles([rt.regress_posterior_path, exp_reg_post_path])

class ABCTeamTestCase(PyMsBayesTestCase):
    BASE_DIR_PATTERN = re.compile(r'pymsbayes-output[-]*\d*')
    DATA_DIR_PATTERN = re.compile(r'd\d+')
    MODEL_DIR_PATTERN = re.compile(r'm\d+')
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')
        self.cfg_path2 = package_paths.data_path('4pairs_1locus_maxt5.cfg')
        self.seed = GLOBAL_RNG.randint(1, 999999999)
        self.rng = random.Random()
        self.rng.seed(self.seed)
        self.output_dir = self.get_test_subdir(prefix='abc-team-test-')

    def tearDown(self):
        for base, dirs, files in os.walk(self.output_dir):
            for d in dirs:
                if self.BASE_DIR_PATTERN.match(d) or \
                        self.DATA_DIR_PATTERN.match(d) or \
                        self.MODEL_DIR_PATTERN.match(d) or \
                        d == 'prior-stats-summaries':
                    self.temp_fs._register_dir(os.path.join(base, d))
        self.tear_down()

    def get_result_paths(self, abc_team, observed_idx, model_idx, sim_idx,
            compressed = False):
        prefix = '{0}d{1}-{2}-s{3}-{4}'.format(abc_team.output_prefix,
                observed_idx, abc_team.model_strings[model_idx], sim_idx,
                abc_team.iter_count)
        p = os.path.join(abc_team.model_dirs[observed_idx][model_idx],
                prefix)
        paths = {}
        paths['sample'] = p + '-posterior-sample.txt'
        paths['summary'] = p + '-posterior-summary.txt'
        paths['psi'] = p + '-psi-results.txt'
        paths['model'] = p + '-model-results.txt'
        paths['omega'] = p + '-omega-results.txt'
        paths['div'] = p + '-div-model-results.txt'
        paths['glm-summary'] = p + '-glm-posterior-summary.txt'
        paths['glm-density'] = p + '-glm-posterior-density-estimates.txt'
        if compressed:
            paths['sample'] += '.gz'
            paths['glm-summary'] += '.gz'
            paths['glm-density'] += '.gz'
        return paths

    def test_abc_team(self):
        obs_worker = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        obs_worker.start()

        num_taxon_pairs = 4
        num_prior_samples = 2000
        num_processors = 4
        batch_size = 200
        num_standardizing_samples = 800
        num_posterior_samples = 100
        num_posterior_density_quantiles = 100
        num_batches = num_prior_samples / batch_size

        abct = ABCTeam(
                temp_fs = self.temp_fs,
                observed_stats_files = [obs_worker.prior_stats_path],
                num_taxon_pairs = num_taxon_pairs,
                config_paths = [self.cfg_path],
                num_prior_samples = num_prior_samples,
                num_processors = num_processors,
                num_standardizing_samples = num_standardizing_samples,
                num_posterior_samples = num_posterior_samples,
                num_posterior_density_quantiles = num_posterior_density_quantiles,
                batch_size = batch_size,
                output_dir = self.output_dir,
                output_prefix = self.test_id,
                rng = self.rng,
                abctoolbox_bandwidth = None,
                omega_threshold = 0.01,
                compress = False,
                keep_temps = False,
                global_estimate_only = False)
        self.assertFalse(abct.finished)
        abct.run()
        self.assertTrue(abct.finished)
        self.assertEqual(abct.num_samples_generated, 2000)
        self.assertEqual(abct.num_samples_summarized, 800)
        self.assertEqual(abct.num_samples_processed[1], 2000)

        res = self.get_result_paths(abct, 1, 1, 1)

        self.assertTrue(os.path.isfile(res['sample']))
        self.assertTrue(os.path.isfile(res['summary']))
        self.assertTrue(os.path.isfile(res['model']))
        self.assertTrue(os.path.isfile(res['div']))
        self.assertTrue(os.path.isfile(res['psi']))
        self.assertTrue(os.path.isfile(res['omega']))
        self.assertTrue(os.path.isfile(res['glm-summary']))
        self.assertTrue(os.path.isfile(res['glm-density']))

        self.assertTrue(self.get_number_of_lines(res['sample']),
                num_posterior_samples + 1)
        self.assertTrue(self.get_number_of_lines(res['psi']),
                num_taxon_pairs + 1)
        self.assertTrue(self.get_number_of_lines(res['omega']), 2)
        self.assertTrue(self.get_number_of_lines(res['glm-density']),
                num_posterior_density_quantiles)
        self.assertTrue(self.get_number_of_lines(res['glm-summary']), 20)

        MsBayesWorker.count = 1
        self.rng.seed(self.seed)
        prior_workers = []
        for i in range(num_batches):
            prior_workers.append(MsBayesWorker(
                    temp_fs = self.temp_fs,
                    sample_size = batch_size,
                    config_path = self.cfg_path,
                    model_index = 1,
                    seed = get_random_int(self.rng),
                    schema = 'abctoolbox',
                    tag = 1))
        prior_workers = Manager.run_workers(
                workers = prior_workers,
                num_processors = num_processors)
        post_path = self.get_test_path(prefix='sum-post')
        sum_out_path = self.get_test_path(prefix='summary')
        sum_worker = EuRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_paths = [pw.prior_path for pw in prior_workers],
                posterior_path = post_path,
                num_posterior_samples = 0,
                num_standardizing_samples = num_standardizing_samples,
                summary_out_path = sum_out_path,
                tag = 1)
        sum_worker.start()
        rej_worker = EuRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_paths = [pw.prior_path for pw in prior_workers],
                posterior_path = post_path,
                num_posterior_samples = num_posterior_samples,
                num_standardizing_samples = 0,
                summary_in_path = sum_out_path,
                tag = 1)
        rej_worker.start()
        samples = parse_parameters(post_path)
        ipc = IntegerPartitionCollection(samples['taus'])
        div_models_to_indices = {}
        for i, k in enumerate(ipc.iterkeys()):
            div_models_to_indices[k] = i + 1
        new_post_path = self.get_test_path(prefix='new-post')
        add_div_model_column(post_path, new_post_path,
                div_models_to_indices,
                compresslevel = None)
        
        self.assertSameFiles([res['sample'],
                new_post_path])

        post_out_path = self.get_test_path(prefix='pw-post-out')
        o_prefix = os.path.join(self.temp_fs.base_dir, self.test_id)
        pw = PosteriorWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                posterior_path = post_path,
                num_taxon_pairs = num_taxon_pairs,
                posterior_out_path = post_out_path,
                output_prefix = o_prefix,
                model_indices = [1],
                abctoolbox_bandwidth = None,
                abctoolbox_num_posterior_quantiles = \
                        num_posterior_density_quantiles,
                omega_threshold = 0.01,
                compress = False,
                keep_temps = False,
                tag = 1)
        pw.start()

        self.assertSameFiles([pw.posterior_out_path, new_post_path,
                res['sample']])
        self.assertSameFiles([
                res['glm-density'],
                pw.regress_posterior_path])
        self.assertSameFiles([
                res['glm-summary'],
                pw.regress_summary_path])
        self.assertSameFiles([
                res['div'],
                pw.div_model_results_path])
        self.assertSameFiles([
                res['psi'],
                pw.psi_results_path])
        self.assertSameFiles([
                res['omega'],
                pw.omega_results_path])
        self.assertSameFiles([
                res['summary'],
                pw.posterior_summary_path])

    def test_abc_team_prior_generation(self):
        obs_worker = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        obs_worker.start()

        num_taxon_pairs = 4
        num_prior_samples = 2000
        num_processors = 4
        batch_size = 200
        num_standardizing_samples = 800
        num_posterior_samples = 100
        num_posterior_density_quantiles = 100
        num_batches = num_prior_samples / batch_size

        abct = ABCTeam(
                temp_fs = self.temp_fs,
                observed_stats_files = [obs_worker.prior_stats_path],
                num_taxon_pairs = num_taxon_pairs,
                config_paths = [self.cfg_path],
                num_prior_samples = num_prior_samples,
                num_processors = num_processors,
                num_standardizing_samples = num_standardizing_samples,
                num_posterior_samples = num_posterior_samples,
                num_posterior_density_quantiles = num_posterior_density_quantiles,
                batch_size = batch_size,
                output_dir = self.output_dir,
                output_prefix = self.test_id,
                rng = self.rng,
                abctoolbox_bandwidth = None,
                omega_threshold = 0.01,
                compress = True,
                keep_temps = False,
                global_estimate_only = False,
                generate_prior_samples_only = True)
        self.assertFalse(abct.finished)
        abct.run()
        self.assertTrue(abct.finished)
        self.assertEqual(abct.num_samples_generated, 2000)
        self.assertEqual(abct.num_samples_summarized, 800)

        prior_path = abct._get_prior_path(1)

        MsBayesWorker.count = 1
        self.rng.seed(self.seed)
        prior_workers = []
        for i in range(num_batches):
            prior_workers.append(MsBayesWorker(
                    temp_fs = self.temp_fs,
                    sample_size = batch_size,
                    config_path = self.cfg_path,
                    model_index = 1,
                    seed = get_random_int(self.rng),
                    schema = 'abctoolbox',
                    tag = 1))
        prior_workers = Manager.run_workers(
                workers = prior_workers,
                num_processors = num_processors)
        expected_prior_path = self.get_test_path(prefix='expected-prior')
        merge_prior_files(
                paths = [pw.prior_path for pw in prior_workers],
                dest_path = expected_prior_path)

        self.assertSameUnsortedFiles([prior_path, expected_prior_path])

    def test_abc_team_with_reporting(self):
        obs_worker = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        obs_worker.start()

        num_taxon_pairs = 4
        num_prior_samples = 2000
        num_processors = 2
        batch_size = 200
        num_standardizing_samples = 800
        num_posterior_samples = 100
        num_posterior_density_quantiles = 100
        num_batches = num_prior_samples / batch_size
        reporting_frequency = 2

        abct = ABCTeam(
                temp_fs = self.temp_fs,
                observed_stats_files = [obs_worker.prior_stats_path],
                num_taxon_pairs = num_taxon_pairs,
                config_paths = [self.cfg_path],
                num_prior_samples = num_prior_samples,
                num_processors = num_processors,
                num_standardizing_samples = num_standardizing_samples,
                num_posterior_samples = num_posterior_samples,
                num_posterior_density_quantiles = num_posterior_density_quantiles,
                batch_size = batch_size,
                output_dir = self.output_dir,
                output_prefix = self.test_id,
                rng = self.rng,
                abctoolbox_bandwidth = None,
                omega_threshold = 0.01,
                reporting_frequency = reporting_frequency,
                compress = False,
                keep_temps = False,
                global_estimate_only = False)
        self.assertFalse(abct.finished)
        abct.run()
        self.assertTrue(abct.finished)
        self.assertEqual(abct.num_samples_generated, 2000)
        self.assertEqual(abct.num_samples_summarized, 800)
        self.assertEqual(abct.num_samples_processed[1], 2000)

        res = self.get_result_paths(abct, 1, 1, 1)

        self.assertTrue(os.path.isfile(res['sample']))
        self.assertTrue(os.path.isfile(res['summary']))
        self.assertTrue(os.path.isfile(res['model']))
        self.assertTrue(os.path.isfile(res['div']))
        self.assertTrue(os.path.isfile(res['psi']))
        self.assertTrue(os.path.isfile(res['omega']))
        self.assertTrue(os.path.isfile(res['glm-summary']))
        self.assertTrue(os.path.isfile(res['glm-density']))

        self.assertTrue(self.get_number_of_lines(res['sample']),
                num_posterior_samples + 1)
        self.assertTrue(self.get_number_of_lines(res['psi']),
                num_taxon_pairs + 1)
        self.assertTrue(self.get_number_of_lines(res['omega']), 2)
        self.assertTrue(self.get_number_of_lines(res['glm-density']),
                num_posterior_density_quantiles)
        self.assertTrue(self.get_number_of_lines(res['glm-summary']), 20)

        MsBayesWorker.count = 1
        self.rng.seed(self.seed)
        prior_workers = []
        for i in range(num_batches):
            prior_workers.append(MsBayesWorker(
                    temp_fs = self.temp_fs,
                    sample_size = batch_size,
                    config_path = self.cfg_path,
                    model_index = 1,
                    seed = get_random_int(self.rng),
                    schema = 'abctoolbox',
                    tag = 1))
        prior_workers = Manager.run_workers(
                workers = prior_workers,
                num_processors = num_processors)
        post_path = self.get_test_path(prefix='sum-post')
        sum_out_path = self.get_test_path(prefix='summary')
        sum_worker = EuRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_paths = [pw.prior_path for pw in prior_workers],
                posterior_path = post_path,
                num_posterior_samples = 0,
                num_standardizing_samples = num_standardizing_samples,
                summary_out_path = sum_out_path,
                tag = 1)
        sum_worker.start()
        rej_worker = EuRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_paths = [pw.prior_path for pw in prior_workers],
                posterior_path = post_path,
                num_posterior_samples = num_posterior_samples,
                num_standardizing_samples = 0,
                summary_in_path = sum_out_path,
                tag = 1)
        rej_worker.start()
        samples = parse_parameters(post_path)
        ipc = IntegerPartitionCollection(samples['taus'])
        div_models_to_indices = {}
        for i, k in enumerate(ipc.iterkeys()):
            div_models_to_indices[k] = i + 1
        new_post_path = self.get_test_path(prefix='new-post')
        add_div_model_column(post_path, new_post_path,
                div_models_to_indices,
                compresslevel = None)
        
        self.assertSameFiles([new_post_path, res['sample']])

        post_out_path = self.get_test_path(prefix='pw-post-out')
        o_prefix = os.path.join(self.temp_fs.base_dir, self.test_id)
        pw = PosteriorWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                posterior_path = post_path,
                num_taxon_pairs = num_taxon_pairs,
                posterior_out_path = post_out_path,
                output_prefix = o_prefix,
                model_indices = [1],
                abctoolbox_bandwidth = None,
                abctoolbox_num_posterior_quantiles = \
                        num_posterior_density_quantiles,
                omega_threshold = 0.01,
                compress = False,
                keep_temps = False,
                tag = 1)
        pw.start()

        self.assertSameFiles([pw.posterior_out_path, new_post_path,
                res['sample']])
        self.assertSameFiles([
                res['glm-density'],
                pw.regress_posterior_path])
        self.assertSameFiles([
                res['glm-summary'],
                pw.regress_summary_path])
        self.assertSameFiles([
                res['div'],
                pw.div_model_results_path])
        self.assertSameFiles([
                res['psi'],
                pw.psi_results_path])
        self.assertSameFiles([
                res['omega'],
                pw.omega_results_path])
        self.assertSameFiles([
                res['summary'],
                pw.posterior_summary_path])

    def test_abc_team_repeatability(self):
        obs_worker = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        obs_worker.start()

        num_taxon_pairs = 4
        num_prior_samples = 2000
        num_processors = 4
        batch_size = 200
        num_standardizing_samples = 800
        num_posterior_samples = 100
        num_posterior_density_quantiles = 100
        num_batches = num_prior_samples / batch_size

        abct1 = ABCTeam(
                temp_fs = self.temp_fs,
                observed_stats_files = [obs_worker.prior_stats_path],
                num_taxon_pairs = num_taxon_pairs,
                config_paths = [self.cfg_path],
                num_prior_samples = num_prior_samples,
                num_processors = num_processors,
                num_standardizing_samples = num_standardizing_samples,
                num_posterior_samples = num_posterior_samples,
                num_posterior_density_quantiles = num_posterior_density_quantiles,
                batch_size = batch_size,
                output_dir = self.output_dir,
                output_prefix = self.test_id,
                rng = self.rng,
                abctoolbox_bandwidth = None,
                omega_threshold = 0.01,
                compress = False,
                keep_temps = False,
                global_estimate_only = False)
        self.assertFalse(abct1.finished)
        abct1.run()
        self.assertTrue(abct1.finished)
        self.assertEqual(abct1.num_samples_generated, 2000)
        self.assertEqual(abct1.num_samples_summarized, 800)
        self.assertEqual(abct1.num_samples_processed[1], 2000)

        self.rng.seed(self.seed)

        abct2 = ABCTeam(
                temp_fs = self.temp_fs,
                observed_stats_files = [obs_worker.prior_stats_path],
                num_taxon_pairs = num_taxon_pairs,
                config_paths = [self.cfg_path],
                num_prior_samples = num_prior_samples,
                num_processors = num_processors,
                num_standardizing_samples = num_standardizing_samples,
                num_posterior_samples = num_posterior_samples,
                num_posterior_density_quantiles = num_posterior_density_quantiles,
                batch_size = batch_size,
                output_dir = self.output_dir,
                output_prefix = self.test_id,
                rng = self.rng,
                abctoolbox_bandwidth = None,
                omega_threshold = 0.01,
                reporting_frequency = 2,
                compress = False,
                keep_temps = False,
                global_estimate_only = False)
        self.assertFalse(abct2.finished)
        abct2.run()
        self.assertTrue(abct2.finished)

        res1 = self.get_result_paths(abct1, 1, 1, 1)
        res2 = self.get_result_paths(abct2, 1, 1, 1)

        self.assertSameFiles([
                res1['sample'],
                res2['sample'],
                ])
        self.assertSameFiles([
                res1['glm-density'],
                res2['glm-density'],
                ])
        self.assertSameFiles([
                res1['glm-summary'],
                res2['glm-summary'],
                ])
        self.assertSameFiles([
                res1['div'],
                res2['div'],
                ])
        self.assertSameFiles([
                res1['psi'],
                res2['psi'],
                ])
        self.assertSameFiles([
                res1['omega'],
                res2['omega'],
                ])
        self.assertSameFiles([
                res1['summary'],
                res2['summary'],
                ])

    def test_abc_team_multiple_models(self):
        obs_worker = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        obs_worker.start()

        num_taxon_pairs = 4
        num_prior_samples = 2000
        num_processors = 4
        batch_size = 200
        num_standardizing_samples = 800
        num_posterior_samples = 100
        num_posterior_density_quantiles = 100
        num_batches = num_prior_samples / batch_size

        abct = ABCTeam(
                temp_fs = self.temp_fs,
                observed_stats_files = [obs_worker.prior_stats_path],
                num_taxon_pairs = num_taxon_pairs,
                config_paths = [self.cfg_path, self.cfg_path2],
                num_prior_samples = num_prior_samples,
                num_processors = num_processors,
                num_standardizing_samples = num_standardizing_samples,
                num_posterior_samples = num_posterior_samples,
                num_posterior_density_quantiles = num_posterior_density_quantiles,
                batch_size = batch_size,
                output_dir = self.output_dir,
                output_prefix = self.test_id,
                rng = self.rng,
                abctoolbox_bandwidth = None,
                omega_threshold = 0.01,
                compress = False,
                keep_temps = False,
                global_estimate_only = False)
        self.assertFalse(abct.finished)
        abct.run()
        self.assertTrue(abct.finished)
        self.assertEqual(abct.num_samples_generated, 4000)
        self.assertEqual(abct.num_samples_summarized, 1600)
        self.assertEqual(abct.num_samples_processed[1], 2000)
        self.assertEqual(abct.num_samples_processed[2], 2000)

        for i in [1, 2, 'combined']:
            res = self.get_result_paths(abct, 1, i, 1)
            self.assertTrue(os.path.isfile(res['sample']))
            self.assertTrue(os.path.isfile(res['summary']))
            self.assertTrue(os.path.isfile(res['model']))
            self.assertTrue(os.path.isfile(res['div']))
            self.assertTrue(os.path.isfile(res['psi']))
            self.assertTrue(os.path.isfile(res['omega']))
            self.assertTrue(os.path.isfile(res['glm-summary']))
            self.assertTrue(os.path.isfile(res['glm-density']))

            self.assertTrue(self.get_number_of_lines(res['sample']),
                    num_posterior_samples + 1)
            self.assertTrue(self.get_number_of_lines(res['psi']),
                    num_taxon_pairs + 1)
            self.assertTrue(self.get_number_of_lines(res['omega']), 2)
            self.assertTrue(self.get_number_of_lines(res['glm-density']),
                    num_posterior_density_quantiles)
            self.assertTrue(self.get_number_of_lines(res['glm-summary']), 20)

    def test_abc_team_multiple_models_global_vs_model_wise(self):
        obs_worker = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        obs_worker.start()

        num_taxon_pairs = 4
        num_prior_samples = 2000
        num_processors = 4
        batch_size = 200
        num_standardizing_samples = 800
        num_posterior_samples = 100
        num_posterior_density_quantiles = 100
        num_batches = num_prior_samples / batch_size

        abct1 = ABCTeam(
                temp_fs = self.temp_fs,
                observed_stats_files = [obs_worker.prior_stats_path],
                num_taxon_pairs = num_taxon_pairs,
                config_paths = [self.cfg_path, self.cfg_path2],
                num_prior_samples = num_prior_samples,
                num_processors = num_processors,
                num_standardizing_samples = num_standardizing_samples,
                num_posterior_samples = num_posterior_samples,
                num_posterior_density_quantiles = num_posterior_density_quantiles,
                batch_size = batch_size,
                output_dir = self.output_dir,
                output_prefix = self.test_id,
                rng = self.rng,
                abctoolbox_bandwidth = None,
                omega_threshold = 0.01,
                compress = False,
                keep_temps = False,
                global_estimate_only = False)
        self.assertFalse(abct1.finished)
        abct1.run()
        self.assertTrue(abct1.finished)

        self.rng.seed(self.seed)

        abct2 = ABCTeam(
                temp_fs = self.temp_fs,
                observed_stats_files = [obs_worker.prior_stats_path],
                num_taxon_pairs = num_taxon_pairs,
                config_paths = [self.cfg_path, self.cfg_path2],
                num_prior_samples = num_prior_samples,
                num_processors = num_processors,
                num_standardizing_samples = num_standardizing_samples,
                num_posterior_samples = num_posterior_samples,
                num_posterior_density_quantiles = num_posterior_density_quantiles,
                batch_size = batch_size,
                output_dir = self.output_dir,
                output_prefix = self.test_id,
                rng = self.rng,
                abctoolbox_bandwidth = None,
                omega_threshold = 0.01,
                compress = False,
                keep_temps = False,
                global_estimate_only = True)
        self.assertFalse(abct2.finished)
        abct2.run()
        self.assertTrue(abct2.finished)
        self.assertEqual(abct2.num_samples_generated, 4000)
        self.assertEqual(abct2.num_samples_summarized, 1600)
        self.assertEqual(abct2.num_samples_processed['combined'], 4000)

        res1 = self.get_result_paths(abct1, 1, 'combined', 1)
        res2 = self.get_result_paths(abct2, 1, 'combined', 1)

        self.assertSameSamples(files = [
                res1['sample'],
                res2['sample']],
                columns_to_ignore = [0],
                header = True,
                places = 5,
                num_mismatches_per_sample = 40,
                num_sample_mismatches = 40)

    def test_abc_team_multiple_models_multiple_obs(self):
        obs_worker1 = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        obs_worker1.start()
        obs_worker2 = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        obs_worker2.start()

        obs_path = self.get_test_path(prefix='obs-sims')
        merge_prior_files([obs_worker1.prior_stats_path,
                obs_worker2.prior_stats_path], obs_path)

        num_taxon_pairs = 4
        num_prior_samples = 2000
        num_processors = 4
        batch_size = 200
        num_standardizing_samples = 800
        num_posterior_samples = 100
        num_posterior_density_quantiles = 100
        num_batches = num_prior_samples / batch_size

        abct = ABCTeam(
                temp_fs = self.temp_fs,
                observed_stats_files = [obs_path],
                num_taxon_pairs = num_taxon_pairs,
                config_paths = [self.cfg_path, self.cfg_path2],
                num_prior_samples = num_prior_samples,
                num_processors = num_processors,
                num_standardizing_samples = num_standardizing_samples,
                num_posterior_samples = num_posterior_samples,
                num_posterior_density_quantiles = num_posterior_density_quantiles,
                batch_size = batch_size,
                output_dir = self.output_dir,
                output_prefix = self.test_id,
                rng = self.rng,
                abctoolbox_bandwidth = None,
                omega_threshold = 0.01,
                compress = False,
                keep_temps = False,
                global_estimate_only = False)
        self.assertFalse(abct.finished)
        abct.run()
        self.assertTrue(abct.finished)
        self.assertEqual(abct.num_samples_generated, 4000)
        self.assertEqual(abct.num_samples_summarized, 1600)
        self.assertEqual(abct.num_samples_processed[1], 4000)
        self.assertEqual(abct.num_samples_processed[2], 4000)

        for i in [1, 2, 'combined']:
            for j in [1, 2]:
                res = self.get_result_paths(abct, 1, i, j)
                self.assertTrue(os.path.isfile(res['sample']))
                self.assertTrue(os.path.isfile(res['summary']))
                self.assertTrue(os.path.isfile(res['model']))
                self.assertTrue(os.path.isfile(res['div']))
                self.assertTrue(os.path.isfile(res['psi']))
                self.assertTrue(os.path.isfile(res['omega']))
                self.assertTrue(os.path.isfile(res['glm-summary']))
                self.assertTrue(os.path.isfile(res['glm-density']))

                self.assertTrue(self.get_number_of_lines(res['sample']),
                        num_posterior_samples + 1)
                self.assertTrue(self.get_number_of_lines(res['psi']),
                        num_taxon_pairs + 1)
                self.assertTrue(self.get_number_of_lines(res['omega']), 2)
                self.assertTrue(self.get_number_of_lines(res['glm-density']),
                        num_posterior_density_quantiles)
                self.assertTrue(self.get_number_of_lines(res['glm-summary']), 20)

    def test_abc_team_multiple_models_multiple_obs_no_global(self):
        obs_worker1 = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        obs_worker1.start()
        obs_worker2 = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        obs_worker2.start()

        obs_path = self.get_test_path(prefix='obs-sims')
        merge_prior_files([obs_worker1.prior_stats_path,
                obs_worker2.prior_stats_path], obs_path)

        num_taxon_pairs = 4
        num_prior_samples = 2000
        num_processors = 4
        batch_size = 200
        num_standardizing_samples = 800
        num_posterior_samples = 100
        num_posterior_density_quantiles = 100
        num_batches = num_prior_samples / batch_size

        abct = ABCTeam(
                temp_fs = self.temp_fs,
                observed_stats_files = [obs_path],
                num_taxon_pairs = num_taxon_pairs,
                config_paths = [self.cfg_path, self.cfg_path2],
                num_prior_samples = num_prior_samples,
                num_processors = num_processors,
                num_standardizing_samples = num_standardizing_samples,
                num_posterior_samples = num_posterior_samples,
                num_posterior_density_quantiles = num_posterior_density_quantiles,
                batch_size = batch_size,
                output_dir = self.output_dir,
                output_prefix = self.test_id,
                rng = self.rng,
                abctoolbox_bandwidth = None,
                omega_threshold = 0.01,
                compress = False,
                keep_temps = False,
                global_estimate_only = False,
                global_estimate = False)
        self.assertFalse(abct.finished)
        abct.run()
        self.assertTrue(abct.finished)
        self.assertEqual(abct.num_samples_generated, 4000)
        self.assertEqual(abct.num_samples_summarized, 1600)
        self.assertEqual(abct.num_samples_processed[1], 4000)
        self.assertEqual(abct.num_samples_processed[2], 4000)

        for i in [1, 2]:
            for j in [1, 2]:
                res = self.get_result_paths(abct, 1, i, j)
                self.assertTrue(os.path.isfile(res['sample']))
                self.assertTrue(os.path.isfile(res['summary']))
                self.assertTrue(os.path.isfile(res['model']))
                self.assertTrue(os.path.isfile(res['div']))
                self.assertTrue(os.path.isfile(res['psi']))
                self.assertTrue(os.path.isfile(res['omega']))
                self.assertTrue(os.path.isfile(res['glm-summary']))
                self.assertTrue(os.path.isfile(res['glm-density']))

                self.assertTrue(self.get_number_of_lines(res['sample']),
                        num_posterior_samples + 1)
                self.assertTrue(self.get_number_of_lines(res['psi']),
                        num_taxon_pairs + 1)
                self.assertTrue(self.get_number_of_lines(res['omega']), 2)
                self.assertTrue(self.get_number_of_lines(res['glm-density']),
                        num_posterior_density_quantiles)
                self.assertTrue(self.get_number_of_lines(res['glm-summary']), 20)

        for obs_idx, d in abct.model_dirs.iteritems():
            self.assertEqual(os.listdir(d['combined']), [])

    def test_abc_team_multiple_models_multiple_separate_obs(self):
        obs_worker1 = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        obs_worker1.start()
        obs_worker2 = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        obs_worker2.start()

        num_taxon_pairs = 4
        num_prior_samples = 2000
        num_processors = 4
        batch_size = 200
        num_standardizing_samples = 800
        num_posterior_samples = 100
        num_posterior_density_quantiles = 100
        num_batches = num_prior_samples / batch_size

        abct = ABCTeam(
                temp_fs = self.temp_fs,
                observed_stats_files = [obs_worker1.prior_stats_path,
                        obs_worker2.prior_stats_path],
                num_taxon_pairs = num_taxon_pairs,
                config_paths = [self.cfg_path, self.cfg_path2],
                num_prior_samples = num_prior_samples,
                num_processors = num_processors,
                num_standardizing_samples = num_standardizing_samples,
                num_posterior_samples = num_posterior_samples,
                num_posterior_density_quantiles = num_posterior_density_quantiles,
                batch_size = batch_size,
                output_dir = self.output_dir,
                output_prefix = self.test_id,
                rng = self.rng,
                abctoolbox_bandwidth = None,
                omega_threshold = 0.01,
                compress = False,
                keep_temps = False,
                global_estimate_only = False)
        self.assertFalse(abct.finished)
        abct.run()
        self.assertTrue(abct.finished)
        self.assertEqual(abct.num_samples_generated, 4000)
        self.assertEqual(abct.num_samples_summarized, 1600)
        self.assertEqual(abct.num_samples_processed[1], 4000)
        self.assertEqual(abct.num_samples_processed[2], 4000)

        base_dir_list = os.listdir(abct.output_dir)
        self.assertTrue('d1' in base_dir_list)
        self.assertTrue('d2' in base_dir_list)
        self.assertTrue('prior-stats-summaries' in base_dir_list)

        for i in [1, 2]:
            for j in [1, 2, 'combined']:
                for k in [1]:
                    res = self.get_result_paths(abct, i, j, k)
                    self.assertTrue(os.path.isfile(res['sample']))
                    self.assertTrue(os.path.isfile(res['summary']))
                    self.assertTrue(os.path.isfile(res['model']))
                    self.assertTrue(os.path.isfile(res['div']))
                    self.assertTrue(os.path.isfile(res['psi']))
                    self.assertTrue(os.path.isfile(res['omega']))
                    self.assertTrue(os.path.isfile(res['glm-density']))

                    self.assertTrue(self.get_number_of_lines(res['sample']),
                            num_posterior_samples + 1)
                    self.assertTrue(self.get_number_of_lines(res['psi']),
                            num_taxon_pairs + 1)
                    self.assertTrue(self.get_number_of_lines(res['omega']), 2)
                    self.assertTrue(self.get_number_of_lines(res['glm-density']),
                            num_posterior_density_quantiles)
                    if os.path.isfile(res['glm-summary']):
                        self.assertTrue(self.get_number_of_lines(
                                res['glm-summary']), 20)

    def test_abc_team_repeatability_multiple_models_multiple_obs(self):
        obs_worker1 = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        obs_worker1.start()
        obs_worker2 = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        obs_worker2.start()

        num_taxon_pairs = 4
        num_prior_samples = 2000
        num_processors = 4
        batch_size = 200
        num_standardizing_samples = 800
        num_posterior_samples = 100
        num_posterior_density_quantiles = 100
        num_batches = num_prior_samples / batch_size

        abct1 = ABCTeam(
                temp_fs = self.temp_fs,
                observed_stats_files = [obs_worker1.prior_stats_path,
                        obs_worker2.prior_stats_path],
                num_taxon_pairs = num_taxon_pairs,
                config_paths = [self.cfg_path, self.cfg_path2],
                num_prior_samples = num_prior_samples,
                num_processors = num_processors,
                num_standardizing_samples = num_standardizing_samples,
                num_posterior_samples = num_posterior_samples,
                num_posterior_density_quantiles = num_posterior_density_quantiles,
                batch_size = batch_size,
                output_dir = self.output_dir,
                output_prefix = self.test_id,
                rng = self.rng,
                abctoolbox_bandwidth = None,
                omega_threshold = 0.01,
                reporting_frequency = None,
                compress = False,
                keep_temps = False,
                global_estimate_only = False)
        self.assertFalse(abct1.finished)
        abct1.run()
        self.assertTrue(abct1.finished)
        self.assertEqual(abct1.num_samples_generated, 4000)
        self.assertEqual(abct1.num_samples_summarized, 1600)
        self.assertEqual(abct1.num_samples_processed[1], 4000)
        self.assertEqual(abct1.num_samples_processed[2], 4000)

        self.rng.seed(self.seed)

        abct2 = ABCTeam(
                temp_fs = self.temp_fs,
                observed_stats_files = [obs_worker1.prior_stats_path,
                        obs_worker2.prior_stats_path],
                num_taxon_pairs = num_taxon_pairs,
                config_paths = [self.cfg_path, self.cfg_path2],
                num_prior_samples = num_prior_samples,
                num_processors = num_processors,
                num_standardizing_samples = num_standardizing_samples,
                num_posterior_samples = num_posterior_samples,
                num_posterior_density_quantiles = num_posterior_density_quantiles,
                batch_size = batch_size,
                output_dir = self.output_dir,
                output_prefix = self.test_id,
                rng = self.rng,
                abctoolbox_bandwidth = None,
                omega_threshold = 0.01,
                reporting_frequency = 1,
                compress = False,
                keep_temps = False,
                global_estimate_only = False)
        self.assertFalse(abct2.finished)
        abct2.run()
        self.assertTrue(abct2.finished)
        self.assertEqual(abct2.num_samples_generated, 4000)
        self.assertEqual(abct2.num_samples_summarized, 1600)
        self.assertEqual(abct2.num_samples_processed[1], 4000)
        self.assertEqual(abct2.num_samples_processed[2], 4000)

        for i in [1, 2]:
            for j in [1, 2, 'combined']:
                for k in [1]:
                    res1 = self.get_result_paths(abct1, i, j, k)
                    res2 = self.get_result_paths(abct2, i, j, k)
                    self.assertSameFiles([res1['sample'], res2['sample']])

    def test_abc_team_multiple_models_multiple_obs_previous_priors(self):
        obs_worker1 = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 2,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        obs_worker1.start()
        obs_worker2 = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 2,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        obs_worker2.start()

        num_taxon_pairs = 4
        num_prior_samples = 2000
        num_processors = 4
        batch_size = 200
        num_standardizing_samples = 800
        num_posterior_samples = 100
        num_posterior_density_quantiles = 100
        num_batches = num_prior_samples / batch_size

        abct1 = ABCTeam(
                temp_fs = self.temp_fs,
                observed_stats_files = [obs_worker1.prior_stats_path,
                        obs_worker2.prior_stats_path],
                num_taxon_pairs = num_taxon_pairs,
                config_paths = [self.cfg_path, self.cfg_path2],
                num_prior_samples = num_prior_samples,
                num_processors = num_processors,
                num_standardizing_samples = num_standardizing_samples,
                num_posterior_samples = num_posterior_samples,
                num_posterior_density_quantiles = num_posterior_density_quantiles,
                batch_size = batch_size,
                output_dir = self.output_dir,
                output_prefix = self.test_id,
                rng = self.rng,
                abctoolbox_bandwidth = None,
                omega_threshold = 0.01,
                reporting_frequency = None,
                compress = False,
                keep_temps = False,
                global_estimate_only = False)
        self.assertFalse(abct1.finished)
        abct1.run()
        self.assertTrue(abct1.finished)
        self.assertEqual(abct1.num_samples_generated, 4000)
        self.assertEqual(abct1.num_samples_summarized, 1600)
        self.assertEqual(abct1.num_samples_processed[1], 8000)
        self.assertEqual(abct1.num_samples_processed[2], 8000)

        self.rng.seed(self.seed)

        abct2 = ABCTeam(
                temp_fs = self.temp_fs,
                observed_stats_files = [obs_worker1.prior_stats_path,
                        obs_worker2.prior_stats_path],
                num_taxon_pairs = num_taxon_pairs,
                config_paths = [self.cfg_path, self.cfg_path2],
                num_prior_samples = num_prior_samples,
                num_processors = num_processors,
                num_standardizing_samples = num_standardizing_samples,
                num_posterior_samples = num_posterior_samples,
                num_posterior_density_quantiles = num_posterior_density_quantiles,
                batch_size = batch_size,
                output_dir = self.output_dir,
                output_prefix = self.test_id,
                rng = self.rng,
                abctoolbox_bandwidth = None,
                omega_threshold = 0.01,
                reporting_frequency = None,
                compress = False,
                keep_temps = False,
                global_estimate_only = False,
                generate_prior_samples_only = True)
        self.assertFalse(abct2.finished)
        abct2.run()
        self.assertTrue(abct2.finished)
        self.assertEqual(abct2.num_samples_generated, 4000)
        self.assertEqual(abct2.num_samples_summarized, 1600)

        abct3 = ABCTeam(
                temp_fs = self.temp_fs,
                observed_stats_files = [obs_worker1.prior_stats_path,
                        obs_worker2.prior_stats_path],
                num_taxon_pairs = num_taxon_pairs,
                previous_prior_dir = abct2.summary_dir,
                num_processors = num_processors,
                num_posterior_samples = num_posterior_samples,
                num_posterior_density_quantiles = num_posterior_density_quantiles,
                output_dir = self.output_dir,
                output_prefix = self.test_id,
                rng = self.rng,
                abctoolbox_bandwidth = None,
                omega_threshold = 0.01,
                reporting_frequency = 1,
                compress = False,
                keep_temps = False,
                global_estimate_only = False,
                generate_prior_samples_only = False)
        self.assertFalse(abct3.finished)
        abct3.run()
        self.assertTrue(abct3.finished)
        self.assertEqual(abct3.num_samples_processed[1], 8000)
        self.assertEqual(abct3.num_samples_processed[2], 8000)

        for i in [1, 2]:
            for j in [1, 2, 'combined']:
                for k in [1, 2]:
                    res1 = self.get_result_paths(abct1, i, j, k)
                    res2 = self.get_result_paths(abct3, i, j, k)
                    self.assertSameFiles([res1['sample'], res2['sample']])

        observed_paths = {
                1: dict(zip(range(1,3),
                        [self.get_test_path(prefix='obs1') for i in range(2)])),
                2: dict(zip(range(1,3),
                        [self.get_test_path(prefix='obs2') for i in range(2)]))}
        with open(obs_worker1.prior_stats_path, 'rU') as obs_stream:
            head = obs_stream.next()
            for i, line in enumerate(obs_stream):
                with open(observed_paths[1][i+1], 'w') as out:
                    out.write(head)
                    out.write(line)
        with open(obs_worker2.prior_stats_path, 'rU') as obs_stream:
            head = obs_stream.next()
            for i, line in enumerate(obs_stream):
                with open(observed_paths[2][i+1], 'w') as out:
                    out.write(head)
                    out.write(line)

        prior_paths = {
                1: glob.glob(os.path.join(abct2.summary_dir,
                        '*m1-prior-sample.txt'))[0],
                2: glob.glob(os.path.join(abct2.summary_dir,
                        '*m2-prior-sample.txt'))[0]}
        summary_paths = {
                1: glob.glob(os.path.join(abct2.summary_dir,
                        '*m1-stat-means-and-std-devs.txt'))[0],
                2: glob.glob(os.path.join(abct2.summary_dir,
                        '*m2-stat-means-and-std-devs.txt'))[0],
                'combined': glob.glob(os.path.join(abct2.summary_dir,
                        '*m12-combined-stat-means-and-std-devs.txt'))[0]}

        tmp_post_paths = {}
        post_paths = {}
        for i in [1, 2]: # observed_file
            tmp_post_paths[i] = {}
            post_paths[i] = {}
            for j in [1, 2, 'combined']: # model
                tmp_post_paths[i][j] = {}
                post_paths[i][j] = {}
                for k in [1, 2]: # simulated dataset
                    tmp_post_paths[i][j][k] = self.get_test_path(
                            prefix='tmp-post-d{0}-m{1}-s{2}-'.format(i,j,k))
                    post_paths[i][j][k] = self.get_test_path(
                            prefix='post-d{0}-m{1}-s{2}-'.format(i,j,k))
        for i in [1, 2]: # observed_file
            for j in [1, 2]: # model
                for k in [1, 2]: # simulated dataset
                    rej_worker = EuRejectWorker(
                            temp_fs = self.temp_fs,
                            observed_path = observed_paths[i][k],
                            prior_paths = [prior_paths[j]],
                            posterior_path = tmp_post_paths[i][j][k],
                            num_posterior_samples = num_posterior_samples,
                            num_standardizing_samples = 0,
                            summary_in_path = summary_paths[j])
                    rej_worker.start()

                    samples = parse_parameters(tmp_post_paths[i][j][k])
                    ipc = IntegerPartitionCollection(samples['taus'])
                    div_models_to_indices = {}
                    for idx, key in enumerate(ipc.iterkeys()):
                        div_models_to_indices[key] = idx + 1
                    add_div_model_column(tmp_post_paths[i][j][k], post_paths[i][j][k],
                            div_models_to_indices,
                            compresslevel = None)

                    res1 = self.get_result_paths(abct1, i, j, k)
                    res2 = self.get_result_paths(abct3, i, j, k)
                    self.assertSameFiles([res1['sample'], res2['sample'],
                            post_paths[i][j][k]])

        for i in [1, 2]: # observed_file
            for j in ['combined']: # model
                for k in [1, 2]:
                    rej_worker = EuRejectWorker(
                            temp_fs = self.temp_fs,
                            observed_path = observed_paths[i][k],
                            prior_paths = [tmp_post_paths[i][1][k], tmp_post_paths[i][2][k]],
                            posterior_path = tmp_post_paths[i][j][k],
                            num_posterior_samples = num_posterior_samples,
                            num_standardizing_samples = 0,
                            summary_in_path = summary_paths['combined'])
                    rej_worker.start()

                    samples = parse_parameters(tmp_post_paths[i][j][k])
                    ipc = IntegerPartitionCollection(samples['taus'])
                    div_models_to_indices = {}
                    for idx, key in enumerate(ipc.iterkeys()):
                        div_models_to_indices[key] = idx + 1
                    add_div_model_column(tmp_post_paths[i][j][k], post_paths[i][j][k],
                            div_models_to_indices,
                            compresslevel = None)

                    res1 = self.get_result_paths(abct1, i, j, k)
                    res2 = self.get_result_paths(abct3, i, j, k)
                    self.assertSameFiles([res1['sample'], res2['sample'],
                            post_paths[i][j][k]])


if __name__ == '__main__':
    unittest.main()

