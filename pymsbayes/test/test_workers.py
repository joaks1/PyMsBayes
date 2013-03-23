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

if __name__ == '__main__':
    unittest.main()

