#! /usr/bin/env python

import unittest
import os
import sys

from pymsbayes.teams import RejectionTeam, ABCTeam, MsRejectionTeam
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

    def test_one_prior_worker(self):
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
        rt = RejectionTeam(
                temp_fs = self.temp_fs,
                prior_workers = [prior_worker],
                observed_path = obs_worker.prior_stats_path,
                num_posterior_samples = 10)
        rt.start()
        self.assertTrue(os.path.exists(rt.posterior_path))
        self.assertEqual(self.get_number_of_lines(rt.posterior_path), 11)
        self.assertEqual(self.get_number_of_header_lines(rt.posterior_path), 1)
        post_path = self.get_test_path(prefix='post-')
        rw = ABCToolBoxRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_path = prior_worker.prior_path,
                num_posterior_samples = 10,
                posterior_path = post_path,
                regression_worker = None,
                max_read_sims = 1000)
        rw.start()
        self.assertSameUnsortedFiles([rt.posterior_path, rw.posterior_path])

    def test_multiple_prior_workers(self):
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
                    sample_size = 100,
                    config_path = self.cfg_path,
                    schema = 'abctoolbox',
                    write_stats_file = False)
            prior_worker.start()
            prior_workers.append(prior_worker)
        prior_path = self.get_test_path(prefix='prior-')
        header_path = self.get_test_path(prefix='prior-header-')
        merge_priors(prior_workers,
                prior_path = prior_path,
                header_path = header_path,
                include_header = True)
        rt = RejectionTeam(
                temp_fs = self.temp_fs,
                prior_workers = prior_workers,
                observed_path = obs_worker.prior_stats_path,
                num_posterior_samples = 10,
                keep_temps = True)
        rt.start()
        self.assertTrue(os.path.exists(rt.posterior_path))
        self.assertEqual(self.get_number_of_lines(rt.posterior_path), 11)
        self.assertEqual(self.get_number_of_header_lines(rt.posterior_path), 1)
        post_path = self.get_test_path(prefix='post-')
        rw = ABCToolBoxRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_path = prior_path,
                num_posterior_samples = 10,
                posterior_path = post_path,
                regression_worker = None,
                max_read_sims = 10000)
        rw.start()
        self.assertSameSamples(files = [rt.posterior_path, rw.posterior_path],
                columns_to_ignore = [],
                header = True,
                places = 4,
                num_mismatches_per_sample = 0,
                num_sample_mismatches = 2)

    def test_multiple_prior_workers_large(self):
        if test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
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
                        sample_size = 2000,
                        config_path = self.cfg_path,
                        schema = 'abctoolbox',
                        write_stats_file = False)
                prior_worker.start()
                prior_workers.append(prior_worker)
            prior_path = self.get_test_path(prefix='prior-')
            header_path = self.get_test_path(prefix='prior-header-')
            merge_priors(prior_workers,
                    prior_path = prior_path,
                    header_path = header_path,
                    include_header = True)
            rt = RejectionTeam(
                    temp_fs = self.temp_fs,
                    prior_workers = prior_workers,
                    observed_path = obs_worker.prior_stats_path,
                    num_posterior_samples = 10,
                    keep_temps = True)
            rt.start()
            self.assertTrue(os.path.exists(rt.posterior_path))
            self.assertEqual(self.get_number_of_lines(rt.posterior_path), 11)
            self.assertEqual(self.get_number_of_header_lines(rt.posterior_path), 1)
            post_path = self.get_test_path(prefix='post-')
            rw = ABCToolBoxRejectWorker(
                    temp_fs = self.temp_fs,
                    observed_path = obs_worker.prior_stats_path,
                    prior_path = prior_path,
                    num_posterior_samples = 10,
                    posterior_path = post_path,
                    regression_worker = None,
                    max_read_sims = 10000)
            rw.start()
            self.assertSameSamples(files = [rt.posterior_path, rw.posterior_path],
                    columns_to_ignore = [],
                    header = True,
                    places = 4,
                    num_mismatches_per_sample = 0,
                    num_sample_mismatches = 2)

class MsRejectionTeamTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')

    def tearDown(self):
        self.tear_down()

    def test_one_prior_worker(self):
        obs_worker = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'msreject',
                include_header = False)
        obs_worker.start()
        prior_worker = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 100,
                config_path = self.cfg_path,
                schema = 'msreject',
                include_header = False)
        prior_worker.start()
        rt = MsRejectionTeam(
                temp_fs = self.temp_fs,
                prior_workers = [prior_worker],
                observed_path = obs_worker.prior_path,
                num_posterior_samples = 10)
        rt.start()
        self.assertTrue(os.path.exists(rt.posterior_path))
        self.assertEqual(self.get_number_of_lines(rt.posterior_path), 11)
        self.assertEqual(self.get_number_of_header_lines(rt.posterior_path), 1)
        post_path = self.get_test_path(prefix='post-')
        rw = MsRejectWorker(
                header = prior_worker.header,
                observed_path = obs_worker.prior_path,
                prior_path = prior_worker.prior_path,
                tolerance = 0.1,
                posterior_path = post_path,
                regression_worker = None)
        rw.start()
        self.assertSameUnsortedFiles([rt.posterior_path, rw.posterior_path])

    def test_multiple_prior_workers(self):
        obs_worker = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'msreject',
                include_header = False)
        obs_worker.start()
        prior_workers = []
        for i in range(2):
            prior_worker = MsBayesWorker(
                    temp_fs = self.temp_fs,
                    sample_size = 100,
                    config_path = self.cfg_path,
                    schema = 'msreject',
                    include_header = False)
            prior_worker.start()
            prior_workers.append(prior_worker)
        prior_path = self.get_test_path(prefix='prior-')
        header_path = self.get_test_path(prefix='prior-header-')
        merge_priors(prior_workers,
                prior_path = prior_path,
                header_path = header_path,
                include_header = False)
        rt = MsRejectionTeam(
                temp_fs = self.temp_fs,
                prior_workers = prior_workers,
                observed_path = obs_worker.prior_path,
                num_posterior_samples = 10,
                keep_temps = True)
        rt.start()
        self.assertTrue(os.path.exists(rt.posterior_path))
        self.assertEqual(self.get_number_of_lines(rt.posterior_path), 11)
        self.assertEqual(self.get_number_of_header_lines(rt.posterior_path), 1)
        post_path = self.get_test_path(prefix='post-')
        rw = MsRejectWorker(
                header = prior_workers[0].header,
                observed_path = obs_worker.prior_path,
                prior_path = prior_path,
                tolerance = 0.05,
                posterior_path = post_path,
                regression_worker = None)
        rw.start()
        self.assertSameSamples(files = [rt.posterior_path, rw.posterior_path],
                columns_to_ignore = [],
                header = True,
                places = 4,
                num_mismatches_per_sample = 0,
                num_sample_mismatches = 1)

    def test_multiple_prior_workers_large(self):
        if test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            obs_worker = MsBayesWorker(
                    temp_fs = self.temp_fs,
                    sample_size = 1,
                    config_path = self.cfg_path,
                    schema = 'msreject',
                    include_header = False)
            obs_worker.start()
            prior_workers = []
            for i in range(2):
                prior_worker = MsBayesWorker(
                        temp_fs = self.temp_fs,
                        sample_size = 2000,
                        config_path = self.cfg_path,
                        schema = 'msreject',
                        include_header = False)
                prior_worker.start()
                prior_workers.append(prior_worker)
            prior_path = self.get_test_path(prefix='prior-')
            header_path = self.get_test_path(prefix='prior-header-')
            merge_priors(prior_workers,
                    prior_path = prior_path,
                    header_path = header_path,
                    include_header = False)
            rt = MsRejectionTeam(
                    temp_fs = self.temp_fs,
                    prior_workers = prior_workers,
                    observed_path = obs_worker.prior_path,
                    num_posterior_samples = 10,
                    keep_temps = True)
            rt.start()
            self.assertTrue(os.path.exists(rt.posterior_path))
            self.assertEqual(self.get_number_of_lines(rt.posterior_path), 11)
            self.assertEqual(self.get_number_of_header_lines(rt.posterior_path), 1)
            post_path = self.get_test_path(prefix='post-')
            rw = MsRejectWorker(
                    header = prior_workers[0].header,
                    observed_path = obs_worker.prior_path,
                    prior_path = prior_path,
                    tolerance = 0.0025,
                    posterior_path = post_path,
                    regression_worker = None)
            rw.start()
            self.assertSameSamples(files = [rt.posterior_path, rw.posterior_path],
                    columns_to_ignore = [],
                    header = True,
                    places = 4,
                    num_mismatches_per_sample = 0,
                    num_sample_mismatches = 0)

if __name__ == '__main__':
    unittest.main()

