#! /usr/bin/env python

import unittest
import os
import sys
import multiprocessing

from pymsbayes.fileio import open
from pymsbayes.manager import Manager
from pymsbayes.workers import (MsBayesWorker, merge_priors,
        assemble_rejection_workers, MsRejectWorker)
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

    def _load_workers(self, sample_size=10, n=1, include_header=False):
        for i in range(n):
            w = MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = sample_size,
                config_path = self.cfg_path,
                schema = 'msreject',
                include_header = include_header)
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

    def test_rejection_workers(self):
        """
        Real-world scale number of simulations. This test will hopefully
        catch problems with running into limit of too many open files if
        there are file handle 'leaks'.
        """
        if test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            n = 200
            sample_size = 10
            nprocessors = 4
            # get prior
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
            prior_path = self.get_test_path(prefix='test-prior-')
            header_path = self.get_test_path(prefix='test-prior-header-')
            merge_priors(workers=workers, prior_path=prior_path,
                    header_path=header_path)
            # get observed
            self._load_workers(sample_size, n, include_header=True)
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
            obs_path = self.get_test_path(prefix='test-obs-')
            obs_header_path = self.get_test_path(prefix='test-obs-header-')
            merge_priors(workers=workers, prior_path=obs_path,
                    header_path=obs_header_path, include_header=True)
            results_dir = self.get_test_subdir(prefix='test-rejection-results-')
            msreject_workers = assemble_rejection_workers(
                    temp_fs = self.temp_fs,
                    observed_sims_file = obs_path,
                    prior_path = prior_path,
                    tolerance = 0.1,
                    results_dir = results_dir,
                    posterior_prefix = self.test_id + '-posterior-',
                    regress = False,
                    rejection_tool = 'msreject')
            self.assertEqual(len(msreject_workers), 2000)
            for w in msreject_workers:
                self.workq.put(w)
            managers = []
            for i in range(nprocessors):
                m = Manager(work_queue = self.workq,
                    result_queue = self.resultq)
                managers.append(m)
                m.start()
            workers = []
            for i in range(len(msreject_workers)):
                workers.append(self.resultq.get())
            for m in managers:
                m.join()
            for w in workers:
                self.assertTrue(w.finished)
                self.assertEqual(w.exit_code, 0)
                self.assertTrue(os.path.isfile(w.posterior_path))
                self.assertTrue(self.get_number_of_lines(w.posterior_path) > 2)

if __name__ == '__main__':
    unittest.main()

