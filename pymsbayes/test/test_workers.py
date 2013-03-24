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
        self.assertTrue(os.path.isfile(w.header_path))
        self.assertEqual(w.header,
                open(w.header_path, 'rU').read().strip().split())
        expected_p_indices, expected_s_indices = self.get_expected_indices(
                num_pairs = num_pairs,
                dummy_column = True,
                parameters_reported = True)
        self.assertEqual(w.stat_indices, expected_s_indices)
        self.assertEqual(w.parameter_indices, expected_p_indices)
        self.assertTrue(self.prior_file_is_valid(w.prior_path,
               num_of_rows = sample_size,
               num_of_columns = len(expected_p_indices + expected_s_indices)+1))

    def test_simple(self):
        w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'msreject')
        self.assertIsInstance(w, workers.MsBayesWorker)
        self.assertFalse(w.finished)
        w.start()
        self._assert_success(w, 4, 10)

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

class PriorMergeTestCase(PyMsBayesTestCase):
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

class AssembleMsRejectWorkersTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')
        self.results_dir = self.get_test_subdir(
                prefix='test-posteriors')
        self.posterior_prefix = self.test_id + '-posterior'

    def tearDown(self):
        self.tear_down()

    def _get_prior(self, n=100, observed=False):
        w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = n,
                config_path = self.cfg_path,
                schema = 'msreject',
                report_parameters = True,
                observed = observed)
        w.start()
        return w
    
    def test_simple(self):
        prior_worker = self._get_prior()
        obs_worker = self._get_prior(n=10, observed=True)
        jobs = workers.assemble_msreject_workers(
                temp_fs = self.temp_fs,
                observed_sims_file = obs_worker.prior_path,
                prior_path = prior_worker.prior_path,
                tolerance = 0.05,
                results_dir = self.results_dir,
                posterior_prefix = self.posterior_prefix)
        for j in jobs:
            self.assertFalse(j.finished)
            j.start()
        for j in jobs:
            self.assertTrue(j.finished)
            self.assertTrue(os.path.isfile(j.posterior_path))
            self.assertEqual(self.get_number_of_lines(j.posterior_path), 6)

if __name__ == '__main__':
    unittest.main()

