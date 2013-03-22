#! /usr/bin/env python

import unittest
import os
import sys
import multiprocessing

from pymsbayes.manager import Manager
from pymsbayes.workers import MsBayesWorker
from pymsbayes.test.support import package_paths
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.test import TestLevel, test_enabled
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class ManagerTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')
        self.workq = multiprocessing.Queue()
        self.resultq = multiprocessing.Queue()

    def tearDown(self):
        self.tear_down()

    def _load_workers(self, sample_size=10, n=1):
        for i in range(n):
            w = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = sample_size,
                config_path = self.cfg_path,
                schema = 'msreject')
            self.workq.put(w)

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

    def test_simple_serial(self):
        n = 4
        self._load_workers(n=n)
        m = Manager(work_queue = self.workq,
                result_queue = self.resultq)
        self.assertFalse(m.is_alive())
        m.start()
        self.assertTrue(m.is_alive())
        w = self.resultq.get()
        m.join()
        self.assertFalse(m.is_alive())
        self.assertEqual(m.exitcode, 0)
        self._assert_success(w, 4, 10)

    def test_simple_parallel(self):
        n = 10
        sample_size = 10
        nprocessors = 4
        self._load_workers(sample_size, n)
        managers = []
        for i in range(nprocessors):
            m = Manager(work_queue = self.workq,
                result_queue = self.resultq)
            managers.append(m)
            m.start()
        workers = []
        for i in range(n):
            workers.append(self.resultq.get())
        for m in managers:
            m.join()
        for w in workers:
            self._assert_success(w, 4, sample_size)

    def test_prior_validity(self):
        if test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            n = 100
            sample_size = 100
            nprocessors = 4
            self._load_workers(sample_size, n)
            managers = []
            for i in range(nprocessors):
                m = Manager(work_queue = self.workq,
                    result_queue = self.resultq)
                managers.append(m)
                m.start()
            workers = []
            for i in range(n):
                workers.append(self.resultq.get())
            for m in managers:
                m.join()
            for w in workers:
                self._assert_success(w, 4, sample_size)
            self.assertPriorIsValid(workers, 0)

if __name__ == '__main__':
    unittest.main()

