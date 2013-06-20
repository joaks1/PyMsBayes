#! /usr/bin/env python

import unittest
import os
import sys

from pymsbayes.teams import *
from pymsbayes.workers import (MsBayesWorker, ABCToolBoxRejectWorker,
        MsRejectWorker, merge_priors)
from pymsbayes.test.support import package_paths
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.test import TestLevel, test_enabled
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


    # def test_multiple_prior_workers(self):
    #     obs_worker = MsBayesWorker(
    #             temp_fs = self.temp_fs,
    #             sample_size = 1,
    #             config_path = self.cfg_path,
    #             schema = 'abctoolbox',
    #             write_stats_file = True)
    #     obs_worker.start()
    #     prior_workers = []
    #     for i in range(2):
    #         prior_worker = MsBayesWorker(
    #                 temp_fs = self.temp_fs,
    #                 sample_size = 100,
    #                 config_path = self.cfg_path,
    #                 schema = 'abctoolbox',
    #                 write_stats_file = False)
    #         prior_worker.start()
    #         prior_workers.append(prior_worker)
    #     prior_path = self.get_test_path(prefix='prior-')
    #     header_path = self.get_test_path(prefix='prior-header-')
    #     merge_priors(prior_workers,
    #             prior_path = prior_path,
    #             header_path = header_path,
    #             include_header = True)
    #     rt = RejectionTeam(
    #             temp_fs = self.temp_fs,
    #             prior_workers = prior_workers,
    #             observed_path = obs_worker.prior_stats_path,
    #             num_posterior_samples = 10,
    #             keep_temps = True)
    #     rt.start()
    #     self.assertTrue(os.path.exists(rt.posterior_path))
    #     self.assertEqual(self.get_number_of_lines(rt.posterior_path), 11)
    #     self.assertEqual(self.get_number_of_header_lines(rt.posterior_path), 1)
    #     post_path = self.get_test_path(prefix='post-')
    #     rw = ABCToolBoxRejectWorker(
    #             temp_fs = self.temp_fs,
    #             observed_path = obs_worker.prior_stats_path,
    #             prior_path = prior_path,
    #             num_posterior_samples = 10,
    #             posterior_path = post_path,
    #             regression_worker = None,
    #             max_read_sims = 10000)
    #     rw.start()
    #     self.assertSameSamples(files = [rt.posterior_path, rw.posterior_path],
    #             columns_to_ignore = [],
    #             header = True,
    #             places = 4,
    #             num_mismatches_per_sample = 0,
    #             num_sample_mismatches = 2)

    # def test_multiple_prior_workers_large(self):
    #     if test_enabled(
    #             level = TestLevel.EXHAUSTIVE,
    #             log = _LOG,
    #             module_name = '.'.join([self.__class__.__name__,
    #                     sys._getframe().f_code.co_name])):
    #         obs_worker = MsBayesWorker(
    #                 temp_fs = self.temp_fs,
    #                 sample_size = 1,
    #                 config_path = self.cfg_path,
    #                 schema = 'abctoolbox',
    #                 write_stats_file = True)
    #         obs_worker.start()
    #         prior_workers = []
    #         for i in range(2):
    #             prior_worker = MsBayesWorker(
    #                     temp_fs = self.temp_fs,
    #                     sample_size = 2000,
    #                     config_path = self.cfg_path,
    #                     schema = 'abctoolbox',
    #                     write_stats_file = False)
    #             prior_worker.start()
    #             prior_workers.append(prior_worker)
    #         prior_path = self.get_test_path(prefix='prior-')
    #         header_path = self.get_test_path(prefix='prior-header-')
    #         merge_priors(prior_workers,
    #                 prior_path = prior_path,
    #                 header_path = header_path,
    #                 include_header = True)
    #         rt = RejectionTeam(
    #                 temp_fs = self.temp_fs,
    #                 prior_workers = prior_workers,
    #                 observed_path = obs_worker.prior_stats_path,
    #                 num_posterior_samples = 10,
    #                 keep_temps = True)
    #         rt.start()
    #         self.assertTrue(os.path.exists(rt.posterior_path))
    #         self.assertEqual(self.get_number_of_lines(rt.posterior_path), 11)
    #         self.assertEqual(self.get_number_of_header_lines(rt.posterior_path), 1)
    #         post_path = self.get_test_path(prefix='post-')
    #         rw = ABCToolBoxRejectWorker(
    #                 temp_fs = self.temp_fs,
    #                 observed_path = obs_worker.prior_stats_path,
    #                 prior_path = prior_path,
    #                 num_posterior_samples = 10,
    #                 posterior_path = post_path,
    #                 regression_worker = None,
    #                 max_read_sims = 10000)
    #         rw.start()
    #         self.assertSameSamples(files = [rt.posterior_path, rw.posterior_path],
    #                 columns_to_ignore = [],
    #                 header = True,
    #                 places = 4,
    #                 num_mismatches_per_sample = 0,
    #                 num_sample_mismatches = 2)


if __name__ == '__main__':
    unittest.main()

