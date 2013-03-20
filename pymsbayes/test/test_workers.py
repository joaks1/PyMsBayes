#! /usr/bin/env python

import unittest
import os

from pymsbayes import workers
from pymsbayes.test.support import package_paths
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class MsBayesWorkerTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')

    def tearDown(self):
        self.tear_down()

    def test_tempfs_error(self):
        kwargs = {'sample_size': 10,
                'config_path': self.cfg_path,
                'schema': 'bogus'}
        self.assertRaises(ValueError, workers.MsBayesWorker, **kwargs)

    def test_schema_error(self):
        kwargs = {'temp_fs': self.temp_fs,
                'sample_size': 10,
                'config_path': self.cfg_path,
                'schema': 'bogus'}
        self.assertRaises(ValueError, workers.MsBayesWorker, **kwargs)

    def _assert_success(self, w, num_pairs, sample_size):
        self.assertFalse(w.is_alive())
        self.assertTrue(w.finished)
        self.assertEqual(0, w.exitcode)
        self.assertEqual(0, w.subprocess_exit_code)
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
        w.join()
        w.finish()
        self._assert_success(w, 4, 10)

    def test_simple_batch(self):
        jobs = []
        for i in range(5):
            w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'msreject')
            jobs.append(w)
            self.assertFalse(w.finished)
            self.assertFalse(w.is_alive())
        for w in jobs:
            w.start()
            self.assertFalse(w.finished)
            self.assertTrue(w.is_alive())
        for w in jobs:
            self.assertFalse(w.finished)
            w.join()
            self.assertFalse(w.finished)
            self.assertFalse(w.is_alive())
            w.finish()
            self.assertTrue(w.finished)
            self.assertFalse(w.is_alive())
        for w in jobs:
            self._assert_success(w, 4, 10)


if __name__ == '__main__':
    unittest.main()

