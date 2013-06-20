#! /usr/bin/env python

import unittest
import os
import sys
import random

from pymsbayes.teams import *
from pymsbayes.workers import (MsBayesWorker, ABCToolBoxRejectWorker,
        MsRejectWorker, merge_priors)
from pymsbayes.test.support import package_paths
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.test import TestLevel, test_enabled, GLOBAL_RNG
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
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')
        self.seed = GLOBAL_RNG.randint(1, 999999999)
        self.rng = random.Random()
        self.rng.seed(self.seed)
        self.output_dir = self.get_test_subdir(prefix='abc-team-test-')

    def tearDown(self):
        self.tear_down()
    def test_abc_team(self):
        obs_worker = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        obs_worker.start()

        abct = ABCTeam(
                temp_fs = self.temp_fs,
                observed_sims_file = obs_worker.prior_stats_path,
                num_taxon_pairs = 4,
                model_indices_to_config_paths = {1: self.cfg_path}
                num_prior_samples = 10000,
                num_processors = 4,
                num_standardizing_samples = 1000,
                num_posterior_samples = 100,
                num_posterior_density_quantiles = 100,
                batch_size = 1000,
                output_dir = self.output_dir,
                output_prefix = self.test_id,
                rng = self.rng,
                abctoolbox_bandwidth = None,
                omega_threshold = 0.01,
                compress = False,
                keep_temps = True,
                global_estimate_only = False)
        self.assertFalse(abct.finished)
        abct.run()
        self.assertTrue(abct.finished)


if __name__ == '__main__':
    unittest.main()

