#! /usr/bin/env python

import unittest
import os
import sys

from pymsbayes import workers
from pymsbayes.test.support import package_paths
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.test import TestLevel, test_enabled
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class MsBayesWorkerTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')

    def tearDown(self):
        self.tear_down()

    def test_schema_error(self):
        kwargs = {'temp_fs': self.temp_fs,
                'sample_size': 10,
                'config_path': self.cfg_path,
                'schema': 'bogus'}
        self.assertRaises(ValueError, workers.MsBayesWorker, **kwargs)

    def _assert_success(self, w, num_pairs, sample_size):
        self.assertTrue(w.finished)
        self.assertEqual(0, w.exit_code)
        self.assertTrue(os.path.isdir(w.output_dir))
        self.assertTrue(os.path.isfile(w.prior_path))
        if w.prior_stats_path:
            self.assertTrue(os.path.isfile(w.prior_stats_path))
        self.assertTrue(os.path.isfile(w.header_path))
        self.assertEqual(w.header,
                open(w.header_path, 'rU').read().strip().split())
        dummy_column = True
        if w.schema.startswith('abctoolbox'):
            dummy_column = False
        expected_p_indices, expected_s_indices = self.get_expected_indices(
                num_pairs = num_pairs,
                dummy_column = dummy_column,
                parameters_reported = True)
        self.assertEqual(w.stat_indices, expected_s_indices)
        self.assertEqual(w.parameter_indices, expected_p_indices)
        self.assertTrue(self.prior_file_is_valid(w.prior_path,
                num_of_samples = sample_size,
                num_of_columns = len(expected_p_indices + \
                        expected_s_indices) + int(dummy_column)))
        if w.prior_stats_path:
            self.assertTrue(self.prior_file_is_valid(w.prior_stats_path,
                num_of_samples = sample_size,
                num_of_columns = len(expected_s_indices)))

    def test_schema_msreject(self):
        w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'msreject')
        self.assertIsInstance(w, workers.MsBayesWorker)
        self.assertFalse(w.finished)
        w.start()
        self._assert_success(w, 4, 10)

    def test_schema_abctoolbox(self):
        w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'abctoolbox')
        self.assertIsInstance(w, workers.MsBayesWorker)
        self.assertFalse(w.finished)
        w.start()
        self._assert_success(w, 4, 10)

    def test_schema_abctoolbox_write_stats(self):
        w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        self.assertIsInstance(w, workers.MsBayesWorker)
        self.assertFalse(w.finished)
        w.start()
        self._assert_success(w, 4, 10)

    def test_staging_dir(self):
        staging_dir = self.get_test_subdir(prefix='test-staging-dir-')
        self.assertTrue(os.path.isdir(staging_dir))
        self.assertEqual(os.listdir(staging_dir), [])
        w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'msreject',
                staging_dir = staging_dir)
        self.assertIsInstance(w, workers.MsBayesWorker)
        self.assertFalse(w.finished)
        w.start()
        self._assert_success(w, 4, 10)
        self.assertTrue(os.path.isdir(staging_dir))
        self.assertEqual(os.listdir(staging_dir), [])

    def test_repeatability(self):
        jobs = []
        for i in range(4):
            w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'msreject',
                seed = 1111111)
            jobs.append(w)
        for w in jobs:
            w.start()
        for w in jobs:
            self._assert_success(w, 4, 10)
        self.assertSameFiles([j.prior_path for j in jobs])
        self.assertSameFiles([j.header_path for j in jobs])

class MergePriorTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')

    def tearDown(self):
        self.tear_down()

    def _has_header(self, path):
        first_line = open(path, 'rU').next()
        if workers.HEADER_PATTERN.match(first_line):
            return True
        return False

    def _correct_n_lines(self, path, n):
        nlines = 0
        f = open(path, 'rU')
        for l in f:
            nlines += 1
        if n == nlines:
            return True
        return False

    def test_no_header_explicit_header_path(self):
        jobs = []
        for i in range(4):
            w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'msreject')
            jobs.append(w)
        for w in jobs:
            w.start()
        ppath = self.get_test_path(prefix='merged_prior', create=False)
        hpath = self.get_test_path(prefix='merged_header', create=False)
        workers.merge_priors(workers = jobs,
                prior_path = ppath,
                header_path = hpath,
                include_header = False)
        self.assertTrue(os.path.exists(ppath))
        self.assertTrue(os.path.exists(hpath))
        self.assertFalse(self._has_header(ppath))
        self.assertTrue(self._has_header(hpath))
        self.assertTrue(self._correct_n_lines(ppath, 40))
        self.assertTrue(self._correct_n_lines(hpath, 1))

    def test_no_header_no_explicit_header_path(self):
        jobs = []
        for i in range(4):
            w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'msreject')
            jobs.append(w)
        for w in jobs:
            w.start()
        ppath = self.get_test_path(prefix='merged_prior', create=False)
        hpath = ppath + '.header'
        self.temp_fs._register_file(hpath)
        workers.merge_priors(workers = jobs,
                prior_path = ppath,
                header_path = hpath,
                include_header = False)
        self.assertTrue(os.path.exists(ppath))
        self.assertTrue(os.path.exists(hpath))
        self.assertFalse(self._has_header(ppath))
        self.assertTrue(self._has_header(hpath))
        self.assertTrue(self._correct_n_lines(ppath, 40))
        self.assertTrue(self._correct_n_lines(hpath, 1))

    def test_header_explicit_header_path(self):
        jobs = []
        for i in range(4):
            w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'msreject')
            jobs.append(w)
        for w in jobs:
            w.start()
        ppath = self.get_test_path(prefix='merged_prior', create=False)
        hpath = self.get_test_path(prefix='merged_header', create=False)
        workers.merge_priors(workers = jobs,
                prior_path = ppath,
                header_path = hpath,
                include_header = True)
        self.assertTrue(os.path.exists(ppath))
        self.assertTrue(os.path.exists(hpath))
        self.assertTrue(self._has_header(ppath))
        self.assertTrue(self._has_header(hpath))
        self.assertTrue(self._correct_n_lines(ppath, 41))
        self.assertTrue(self._correct_n_lines(hpath, 1))

    def test_merge_priors_with_headers(self):
        jobs = []
        for i in range(4):
            w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'msreject',
                include_header = True)
            jobs.append(w)
        for w in jobs:
            w.start()
        for w in jobs:
            self.assertEqual(self.get_number_of_lines(w.prior_path), 11)
            self.assertEqual(self.get_number_of_header_lines(w.prior_path), 1)
        ppath = self.get_test_path(prefix='merged_prior', create=False)
        hpath = self.get_test_path(prefix='merged_header', create=False)
        workers.merge_priors(workers = jobs,
                prior_path = ppath,
                header_path = hpath,
                include_header = True)
        self.assertTrue(os.path.exists(ppath))
        self.assertTrue(os.path.exists(hpath))
        self.assertTrue(self._has_header(ppath))
        self.assertTrue(self._has_header(hpath))
        self.assertTrue(self._correct_n_lines(ppath, 41))
        self.assertTrue(self._correct_n_lines(hpath, 1))
        self.assertEqual(self.get_number_of_header_lines(ppath), 1)
        self.assertEqual(self.get_number_of_header_lines(hpath), 1)

class MergePriorFilesTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')

    def tearDown(self):
        self.tear_down()

    def test_merge_prior_files_with_headers(self):
        jobs = []
        for i in range(4):
            w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'msreject',
                include_header = True)
            jobs.append(w)
        for w in jobs:
            w.start()
        for w in jobs:
            self.assertEqual(self.get_number_of_lines(w.prior_path), 11)
            self.assertEqual(self.get_number_of_header_lines(w.prior_path), 1)
        ppath = self.get_test_path(prefix='merged_prior', create=False)
        workers.merge_prior_files(
                paths = [w.prior_path for w in jobs],
                dest_path = ppath)
        self.assertTrue(os.path.exists(ppath))
        self.assertEqual(self.get_number_of_header_lines(ppath), 1)
        self.assertEqual(self.get_number_of_lines(ppath), 41)

    def test_merge_prior_files_no_headers(self):
        jobs = []
        for i in range(4):
            w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'msreject',
                include_header = False)
            jobs.append(w)
        for w in jobs:
            w.start()
        for w in jobs:
            self.assertEqual(self.get_number_of_lines(w.prior_path), 10)
            self.assertEqual(self.get_number_of_header_lines(w.prior_path), 0)
        ppath = self.get_test_path(prefix='merged_prior', create=False)
        workers.merge_prior_files(
                paths = [w.prior_path for w in jobs],
                dest_path = ppath)
        self.assertTrue(os.path.exists(ppath))
        self.assertEqual(self.get_number_of_header_lines(ppath), 0)
        self.assertEqual(self.get_number_of_lines(ppath), 40)

class MsRejectWorkerTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')

    def tearDown(self):
        self.tear_down()

    def test_reject_no_parameters(self):
        self.msbayes_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'msreject',
                report_parameters = False)
        self.msbayes_worker.start()
        self.prior_path = self.msbayes_worker.prior_path
        self.header = self.msbayes_worker.header
        self.observed_path = self.get_test_path(prefix='test-obs',
                create=True)
        self.posterior_path = self.get_test_path(prefix='test-post',
                create=False)
        prior = open(self.prior_path, 'rU')
        obs = open(self.observed_path, 'w')
        observed = prior.next().strip().split()
        for i in (set(range(len(self.header))) - \
                set(self.msbayes_worker.stat_indices)):
            observed[i] = '999999999'
        obs.write('{0}\t\n'.format('\t'.join(observed)))
        obs.close()
        prior.close()

        w = workers.MsRejectWorker(
            header = self.header,
            observed_path = self.observed_path,
            prior_path = self.prior_path,
            tolerance = 0.1,
            posterior_path = self.posterior_path)
        self.assertFalse(w.finished)
        w.start()
        self.assertTrue(w.finished)
        self.assertEqual(w.exit_code, 0)
        self.assertTrue(os.path.isfile(w.posterior_path))
        post = open(self.posterior_path, 'rU')
        post_header = post.next().strip().split()
        posterior = post.next().strip().split()
        post.close()
        self.assertEqual(self.msbayes_worker.header, w.header)
        pr = open(self.prior_path, 'rU')
        expected_post = pr.next().strip().split()
        pr.close()
        self.assertEqual(posterior, expected_post)

    def test_reject_with_parameters(self):
        self.msbayes_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'msreject',
                report_parameters = True)
        self.msbayes_worker.start()
        self.prior_path = self.msbayes_worker.prior_path
        self.header = self.msbayes_worker.header
        self.observed_path = self.get_test_path(prefix='test-obs',
                create=True)
        self.posterior_path = self.get_test_path(prefix='test-post',
                create=False)
        prior = open(self.prior_path, 'rU')
        obs = open(self.observed_path, 'w')
        observed = prior.next().strip().split()
        for i in (set(range(len(self.header))) - \
                set(self.msbayes_worker.stat_indices)):
            observed[i] = '999999999'
        obs.write('{0}\t\n'.format('\t'.join(observed)))
        obs.close()
        prior.close()

        w = workers.MsRejectWorker(
            header = self.header,
            observed_path = self.observed_path,
            prior_path = self.prior_path,
            tolerance = 0.1,
            posterior_path = self.posterior_path)
        self.assertFalse(w.finished)
        w.start()
        self.assertTrue(w.finished)
        self.assertEqual(w.exit_code, 0)
        self.assertTrue(os.path.isfile(w.posterior_path))
        post = open(self.posterior_path, 'rU')
        post_header = post.next().strip().split()
        posterior = post.next().strip().split()
        post.close()
        self.assertEqual(self.msbayes_worker.header, w.header)
        pr = open(self.prior_path, 'rU')
        expected_post = pr.next().strip().split()
        pr.close()
        self.assertEqual(posterior, expected_post)

    def test_reject_with_parameters_posterier_sample(self):
        self.msbayes_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'msreject',
                report_parameters = True)
        self.msbayes_worker.start()
        self.prior_path = self.msbayes_worker.prior_path
        self.header = self.msbayes_worker.header
        self.observed_path = self.get_test_path(prefix='test-obs',
                create=True)
        self.posterior_path = self.get_test_path(prefix='test-post',
                create=False)
        prior = open(self.prior_path, 'rU')
        obs = open(self.observed_path, 'w')
        obs_sample = prior.next()
        observed = obs_sample.strip().split()
        other_sample = prior.next()
        for i in (set(range(len(self.header))) - \
                set(self.msbayes_worker.stat_indices)):
            observed[i] = '999999999'
        obs.write('{0}\t\n'.format('\t'.join(observed)))
        obs.close()
        self.new_prior_path = self.get_test_path(prefix='test-prior',
                create=True)
        new_prior = open(self.new_prior_path, 'w')
        for i in range(2):
            new_prior.write(obs_sample)
        for i in range(10):
            new_prior.write(other_sample)
        for line in prior:
            new_prior.write(line)
        new_prior.close()
        prior.close()

        w = workers.MsRejectWorker(
            header = self.header,
            observed_path = self.observed_path,
            prior_path = self.new_prior_path,
            tolerance = 0.1,
            posterior_path = self.posterior_path)
        self.assertFalse(w.finished)
        w.start()
        self.assertTrue(w.finished)
        self.assertEqual(w.exit_code, 0)
        self.assertTrue(os.path.isfile(w.posterior_path))
        post = open(self.posterior_path, 'rU')
        post_header = post.next().strip().split()
        posterior1 = post.next().strip().split()
        posterior2 = post.next().strip().split()
        post.close()
        self.assertEqual(self.msbayes_worker.header, w.header)
        pr = open(self.prior_path, 'rU')
        expected_post = pr.next().strip().split()
        pr.close()
        self.assertEqual(posterior1, expected_post)
        self.assertEqual(posterior2, expected_post)

class ABCToolBoxRejectWorkerTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')

    def tearDown(self):
        self.tear_down()

    def test_rejection_with_parameters(self):
        prior_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 100,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                report_parameters = True)
        prior_worker.start()
        obs_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True,
                report_parameters = True)
        obs_worker.start()

        post_path = self.get_test_path(prefix='test-posterior-')
        reject_worker = workers.ABCToolBoxRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_path = prior_worker.prior_path,
                num_posterior_samples = 10,
                posterior_path = post_path,
                regression_worker = None,
                exe_path = None,
                stdout_path = None,
                stderr_path = None,
                keep_temps = False,
                max_read_sims = 1000)
        self.assertFalse(reject_worker.finished)
        reject_worker.start()
        self.assertTrue(reject_worker.finished)
        self.assertTrue(os.path.isfile(reject_worker.posterior_path))
        self.assertEqual(self.get_number_of_lines(
                reject_worker.posterior_path), 11)
        self.assertEqual(self.get_number_of_header_lines(
                reject_worker.posterior_path), 1)

class ABCToolBoxRegressWorkerTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')

    def tearDown(self):
        self.tear_down()

    def test_regression_from_prior(self):
        post_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 100,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                report_parameters = True)
        post_worker.start()
        obs_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True,
                report_parameters = True)
        obs_worker.start()

        summary_path = self.get_test_path(prefix='test-summary-')
        adjusted_path = self.get_test_path(prefix='test-adjusted-')
        regress_worker = workers.ABCToolBoxRegressWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                posterior_path = post_worker.prior_path,
                parameter_indices = None,
                summary_path = summary_path,
                adjusted_path = adjusted_path,
                exe_path = None,
                stdout_path = None,
                stderr_path = None,
                keep_temps = False,
                bandwidth = None,
                num_posterior_quantiles = 100)
        self.assertFalse(regress_worker.finished)
        regress_worker.start()
        self.assertTrue(regress_worker.finished)
        self.assertTrue(os.path.isfile(regress_worker.summary_path))
        self.assertTrue(os.path.isfile(regress_worker.adjusted_path))
        self.assertEqual(self.get_number_of_lines(
                regress_worker.adjusted_path), 101)

    def test_regression_from_prior_with_indices(self):
        post_worker1 = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 100,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                report_parameters = True,
                model_index=1)
        post_worker2 = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 100,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                report_parameters = True,
                model_index=2)
        post_worker1.start()
        post_worker2.start()
        obs_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True,
                report_parameters = True)
        obs_worker.start()
        post_path = self.get_test_path(prefix='test-posterior')
        head_path = self.get_test_path(prefix='test-posterior-head')
        workers.merge_priors(workers = [post_worker1, post_worker2],
                prior_path = post_path,
                header_path = head_path,
                include_header=True)
        summary_path = self.get_test_path(prefix='test-summary-')
        adjusted_path = self.get_test_path(prefix='test-adjusted-')
        regress_worker = workers.ABCToolBoxRegressWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                posterior_path = post_path,
                parameter_indices = None,
                summary_path = summary_path,
                adjusted_path = adjusted_path,
                exe_path = None,
                stdout_path = None,
                stderr_path = None,
                keep_temps = False,
                bandwidth = None,
                num_posterior_quantiles = 100)
        self.assertFalse(regress_worker.finished)
        regress_worker.start()
        self.assertTrue(regress_worker.finished)
        self.assertTrue(os.path.isfile(regress_worker.summary_path))
        self.assertTrue(os.path.isfile(regress_worker.adjusted_path))
        self.assertEqual(self.get_number_of_lines(
                regress_worker.adjusted_path), 101)

    def test_regression_from_posterior_with_indices(self):
        prior_worker1 = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 100,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                report_parameters = True,
                model_index=1)
        prior_worker2 = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 100,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                report_parameters = True,
                model_index=2)
        prior_worker1.start()
        prior_worker2.start()
        obs_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True,
                report_parameters = True)
        obs_worker.start()
        prior_path = self.get_test_path(prefix='test-prior')
        head_path = self.get_test_path(prefix='test-prior-head')
        post_path = self.get_test_path(prefix='test-post')
        workers.merge_priors(workers = [prior_worker1, prior_worker2],
                prior_path = prior_path,
                header_path = head_path,
                include_header=True)
        reject_worker = workers.ABCToolBoxRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_path = prior_path,
                num_posterior_samples = 100,
                posterior_path = post_path,
                regression_worker = None,
                exe_path = None,
                stdout_path = None,
                stderr_path = None,
                keep_temps = False,
                max_read_sims = 1000)
        self.assertFalse(reject_worker.finished)
        reject_worker.start()
        self.assertTrue(reject_worker.finished)
        self.assertTrue(os.path.isfile(reject_worker.posterior_path))
        self.assertEqual(self.get_number_of_lines(
                reject_worker.posterior_path), 101)
        self.assertEqual(self.get_number_of_header_lines(
                reject_worker.posterior_path), 1)
        summary_path = self.get_test_path(prefix='test-summary-')
        adjusted_path = self.get_test_path(prefix='test-adjusted-')
        regress_worker = workers.ABCToolBoxRegressWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                posterior_path = post_path,
                parameter_indices = None,
                summary_path = summary_path,
                adjusted_path = adjusted_path,
                exe_path = None,
                stdout_path = None,
                stderr_path = None,
                keep_temps = False,
                num_posterior_samples = 100,
                bandwidth = None,
                num_posterior_quantiles = 100)
        self.assertFalse(regress_worker.finished)
        regress_worker.start()
        self.assertTrue(regress_worker.finished)
        self.assertTrue(os.path.isfile(regress_worker.summary_path))
        self.assertTrue(os.path.isfile(regress_worker.adjusted_path))
        self.assertEqual(self.get_number_of_lines(
                regress_worker.adjusted_path), 101)


class AssembleMsRejectWorkersTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')
        self.results_dir = self.get_test_subdir(
                prefix='test-posteriors')
        self.posterior_prefix = self.test_id + '-posterior'

    def tearDown(self):
        self.tear_down()

    def _get_prior(self, n=100, schema='msreject', include_header=False):
        w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = n,
                config_path = self.cfg_path,
                schema = schema,
                report_parameters = True,
                include_header = include_header)
        w.start()
        return w
    
    def test_msreject_no_regression(self):
        prior_worker = self._get_prior()
        obs_worker = self._get_prior(n=10, include_header=True)
        jobs = workers.assemble_rejection_workers(
                temp_fs = self.temp_fs,
                observed_sims_file = obs_worker.prior_path,
                prior_path = prior_worker.prior_path,
                tolerance = 0.05,
                results_dir = self.results_dir,
                posterior_prefix = self.posterior_prefix,
                regress = False,
                rejection_tool = 'msreject')
        self.assertEqual(len(jobs), 10)
        for j in jobs:
            self.assertFalse(j.finished)
            j.start()
        for j in jobs:
            self.assertTrue(j.finished)
            self.assertTrue(os.path.isfile(j.posterior_path))
            self.assertTrue(
                    abs(self.get_number_of_lines(j.posterior_path) - 6) < 2)

    def test_abctoolbox_no_regression(self):
        prior_worker = self._get_prior(schema='abctoolbox')
        obs_worker = self._get_prior(n=10, schema='abctoolbox')
        jobs = workers.assemble_rejection_workers(
                temp_fs = self.temp_fs,
                observed_sims_file = obs_worker.prior_path,
                prior_path = prior_worker.prior_path,
                num_prior_samples = 100,
                num_posterior_samples = 10,
                results_dir = self.results_dir,
                posterior_prefix = self.posterior_prefix,
                regress = False,
                rejection_tool = 'abctoolbox')
        self.assertEqual(len(jobs), 10)
        for j in jobs:
            self.assertFalse(j.finished)
            j.start()
        for j in jobs:
            self.assertTrue(j.finished)
            self.assertTrue(os.path.isfile(j.posterior_path))
            self.assertTrue(self.get_number_of_lines(j.posterior_path), 11)

    def test_msreject_llr(self):
        prior_worker = self._get_prior(n=200)
        obs_worker = self._get_prior(n=2, include_header=True)
        msreject_workers = workers.assemble_rejection_workers(
                temp_fs = self.temp_fs,
                observed_sims_file = obs_worker.prior_path,
                prior_path = prior_worker.prior_path,
                tolerance = 0.5,
                results_dir = self.results_dir,
                posterior_prefix = self.posterior_prefix,
                regress = True,
                rejection_tool = 'msreject',
                regression_method = 'llr')
        self.assertEqual(len(msreject_workers), 2)
        for j in msreject_workers:
            self.assertFalse(j.finished)
            self.assertFalse(j.regression_worker.finished)
            j.start()
        for j in msreject_workers:
            self.assertTrue(j.finished)
            self.assertTrue(j.regression_worker.finished)
            self.assertTrue(os.path.isfile(j.posterior_path))
            self.assertTrue(
                    abs(self.get_number_of_lines(j.posterior_path) - 101) < 6)
            self.assertTrue(os.path.isfile(j.regression_worker.summary_path))
            self.assertTrue(os.path.isfile(j.regression_worker.adjusted_path))
            self.assertTrue(abs(self.get_number_of_lines(
                    j.regression_worker.adjusted_path) - 101) < 6)

class RegressionWorkerTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')
        self.summary_path = self.get_test_path(prefix='post_summary')
        self.adjusted_path = self.get_test_path(prefix='adjusted_samples')

    def tearDown(self):
        self.tear_down()

    def test_regression_no_models(self):
        posterior_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 100,
                config_path = self.cfg_path,
                schema = 'msreject',
                report_parameters = True,
                include_header = True)
        posterior_worker.start()
        obs_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'msreject',
                report_parameters = True,
                include_header = False)
        obs_worker.start()

        regress_worker = workers.RegressionWorker(
                observed_path = obs_worker.prior_path,
                posterior_path = posterior_worker.prior_path,
                summary_path = self.summary_path,
                adjusted_path = self.adjusted_path,
                tolerance = 1.0)
        self.assertFalse(regress_worker.finished)
        regress_worker.start()
        self.assertTrue(regress_worker.finished)
        self.assertEqual(
                self.get_number_of_lines(regress_worker.adjusted_path),
                101)
        result_header = open(regress_worker.adjusted_path, 
                'rU').next().strip().split()
        cont_params = [regress_worker.header[
                i] for i in regress_worker.continuous_parameter_indices]
        self.assertEqual(sorted(result_header), sorted(cont_params))
        disc_params = [regress_worker.header[
                i] for i in regress_worker.discrete_parameter_indices]
        stats = [regress_worker.header[
                i] for i in regress_worker.stat_indices]
        results = self.parse_python_config(regress_worker.summary_path)
        expected_keys = [x.replace('PRI.', '').replace('.',
                '_') for x in cont_params + disc_params] + ['settings']
        self.assertEqual(sorted(results.keys()),
                sorted(expected_keys))
        self.assertEqual(sorted(results['settings']['stats_used']),
                sorted(stats))
        self.assertEqual(sorted(results['settings']['continuous_parameters']),
                sorted(cont_params))
        if len(disc_params) > 1:
            self.assertEqual(sorted(results['settings']['discrete_parameters']),
                    sorted(disc_params))
        else:
            self.assertEqual(results['settings']['discrete_parameters'],
                    disc_params[0])
        
    def test_regression_no_models(self):
        posterior_worker1 = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 100,
                config_path = self.cfg_path,
                schema = 'msreject',
                report_parameters = True,
                model_index = 1,
                include_header = True)
        posterior_worker1.start()
        posterior_worker2 = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 100,
                config_path = self.cfg_path,
                schema = 'msreject',
                report_parameters = True,
                model_index = 2,
                include_header = True)
        posterior_worker2.start()
        post_path = self.get_test_path(prefix='test-posterior')
        head_path = self.get_test_path(prefix='test-posterior-header')
        workers.merge_priors([posterior_worker1, posterior_worker2],
                prior_path = post_path,
                header_path = head_path,
                include_header = True)
        obs_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'msreject',
                report_parameters = True,
                model_index = 1,
                include_header = False)
        obs_worker.start()

        regress_worker = workers.RegressionWorker(
                observed_path = obs_worker.prior_path,
                posterior_path = post_path,
                summary_path = self.summary_path,
                adjusted_path = self.adjusted_path,
                tolerance = 1.0)
        self.assertFalse(regress_worker.finished)
        regress_worker.start()
        self.assertTrue(regress_worker.finished)
        self.assertEqual(
                self.get_number_of_lines(regress_worker.adjusted_path),
                201)
        result_header = open(regress_worker.adjusted_path, 
                'rU').next().strip().split()
        cont_params = [regress_worker.header[
                i] for i in regress_worker.continuous_parameter_indices]
        self.assertEqual(sorted(result_header), sorted(cont_params))
        disc_params = [regress_worker.header[
                i] for i in regress_worker.discrete_parameter_indices]
        stats = [regress_worker.header[
                i] for i in regress_worker.stat_indices]
        results = self.parse_python_config(regress_worker.summary_path)
        # expected_keys = [x.replace('PRI.', '').replace('.',
        #         '_') for x in cont_params + disc_params] + ['settings']
        expected_keys = cont_params + disc_params + ['settings']
        self.assertEqual(sorted(results.keys()),
                sorted(expected_keys))
        self.assertEqual(sorted(results['settings']['stats_used']),
                sorted(stats))
        self.assertEqual(sorted(results['settings']['continuous_parameters']),
                sorted(cont_params))
        if len(disc_params) > 1:
            self.assertEqual(sorted(results['settings']['discrete_parameters']),
                    sorted(disc_params))
        else:
            self.assertEqual(results['settings']['discrete_parameters'],
                    disc_params[0])

if __name__ == '__main__':
    unittest.main()

