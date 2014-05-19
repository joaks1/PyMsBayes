#! /usr/bin/env python

"""
Test suite for end users to run to ensure the bundled tools are working.
"""

import unittest
import os
import sys

from pymsbayes import workers, fileio, manager
from pymsbayes.config import MsBayesConfig
from pymsbayes.test.support import package_paths
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.utils import get_tool_path, MSBAYES_SORT_INDEX
from pymsbayes.utils import parsing, stats
# from pymsbayes.test.test_dmc import DmcTestCase
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__, level = "INFO")

class ToolTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('negros_panay.cfg')
        self.new_cfg_path = package_paths.data_path('negros_panay_new.cfg')
        self.sum_stats_path = package_paths.data_path(
                'negros_panay_sum_stats_sort0.txt')
        self.prior_path = package_paths.data_path(
                'negros_panay_new_prior_sample_10.txt')
        self.posterior_path = package_paths.data_path(
                'posterior-sample.txt.gz')

    def tearDown(self):
        self.tear_down()

    def test_obsSumStats(self):
        ss_path = self.get_test_path(prefix='sum-stats')
        ss_worker = workers.ObsSumStatsWorker(
                temp_fs = self.temp_fs,
                config_path = self.cfg_path,
                output_path = ss_path,
                exe_path = None,
                schema = 'abctoolbox',
                stat_patterns = parsing.DEFAULT_STAT_PATTERNS,
                stderr_path = None,
                tag = None)
        self.assertFalse(ss_worker.finished)
        self.assertEqual(ss_worker.exe_path, get_tool_path('obssumstats'))
        _LOG.info('{0}'.format(ss_worker.exe_path))
        ss_worker.start()
        self.assertTrue(ss_worker.finished)
        self.assertSameFiles([ss_path, self.sum_stats_path])

    def _assert_msbayes_success(self, w, num_pairs, sample_size):
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

    def test_msbayes(self):
        w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 2,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True)
        self.assertIsInstance(w, workers.MsBayesWorker)
        self.assertFalse(w.finished)
        self.assertNotEqual(w.exe_path, get_tool_path('dpp-msbayes'))
        self.assertEqual(w.exe_path, get_tool_path('msbayes'))
        _LOG.info('{0}'.format(w.exe_path))
        w.start()
        self._assert_msbayes_success(w, 9, 2)
        self.assertNotEqual(w.exe_path, get_tool_path('dpp-msbayes'))
        self.assertEqual(w.exe_path, get_tool_path('msbayes'))


    def test_dpp_msbayes(self):
        w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 2,
                config_path = self.new_cfg_path,
                schema = 'abctoolbox')
        self.assertIsInstance(w, workers.MsBayesWorker)
        self.assertFalse(w.finished)
        self.assertEqual(w.exe_path, get_tool_path('dpp-msbayes'))
        self.assertNotEqual(w.exe_path, get_tool_path('msbayes'))
        _LOG.info('{0}'.format(w.exe_path))
        w.start()
        self._assert_msbayes_success(w, 9, 2)
        self.assertEqual(w.exe_path, get_tool_path('dpp-msbayes'))
        self.assertNotEqual(w.exe_path, get_tool_path('msbayes'))


    def test_eureject(self):
        post_path = self.get_test_path(prefix='test-posterior-')
        sum_out_path = self.get_test_path(prefix='test-summary-out-')
        reject_worker = workers.EuRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = self.sum_stats_path,
                prior_paths = [self.prior_path],
                num_posterior_samples = 2,
                num_standardizing_samples = 10,
                summary_in_path = None,
                summary_out_path = sum_out_path,
                posterior_path = post_path,
                regression_worker = None,
                exe_path = None,
                stderr_path = None,
                keep_temps = False,
                tag = 'testcase')
        self.assertFalse(reject_worker.finished)
        self.assertEqual(reject_worker.exe_path, get_tool_path('eureject'))
        _LOG.info('{0}'.format(reject_worker.exe_path))
        reject_worker.start()
        self.assertTrue(reject_worker.finished)
        self.assertTrue(os.path.isfile(reject_worker.posterior_path))
        self.assertEqual(self.get_number_of_lines(
                reject_worker.posterior_path), 3)
        self.assertEqual(self.get_number_of_header_lines(
                reject_worker.posterior_path), 1)
        self.assertTrue(os.path.isfile(reject_worker.summary_out_path))
        self.assertEqual(self.get_number_of_lines(
                reject_worker.summary_out_path), 4)
        self.assertEqual(self.get_number_of_header_lines(
                reject_worker.summary_out_path), 1)
        self.assertEqual(reject_worker.num_processed, 10)
        self.assertEqual(reject_worker.num_summarized, 10)
        self.assertEqual(reject_worker.num_retained, 2)
        self.assertEqual(reject_worker.prior_paths,
                reject_worker.rejection_files)
        self.assertEqual(reject_worker.prior_paths,
                reject_worker.standardizing_files)
        self.assertFalse(os.path.exists(reject_worker.stderr_path))
        self.assertFalse(os.path.exists(reject_worker.output_dir))

    def test_abcestimator(self):
        summary_path = self.get_test_path(prefix='test-summary-')
        post_path = self.get_test_path(prefix='test-posterior-')
        with open(post_path, 'w') as out:
            stream, close = fileio.process_file_arg(self.posterior_path)
            for line in stream:
                out.write('{0}'.format(line))
            if close:
                stream.close()
        regress_posterior_path = self.get_test_path(prefix='test-adjusted-')
        regress_worker = workers.ABCToolBoxRegressWorker(
                temp_fs = self.temp_fs,
                observed_path = self.sum_stats_path,
                posterior_path = post_path,
                parameter_indices = None,
                regress_summary_path = summary_path,
                regress_posterior_path = regress_posterior_path,
                exe_path = None,
                stdout_path = None,
                stderr_path = None,
                keep_temps = False,
                bandwidth = None,
                num_posterior_quantiles = 100)
        self.assertFalse(regress_worker.finished)
        self.assertEqual(regress_worker.exe_path, get_tool_path('abcestimator'))
        _LOG.info('{0}'.format(regress_worker.exe_path))
        regress_worker.start()
        self.assertTrue(regress_worker.finished)
        self.assertTrue(os.path.isfile(regress_worker.regress_summary_path))
        self.assertTrue(os.path.isfile(regress_worker.regress_posterior_path))
        self.assertEqual(self.get_number_of_lines(
                regress_worker.regress_posterior_path), 101)


class MultiprocessingTestCase(unittest.TestCase):
    def setUp(self):
        self.new_cfg_path = package_paths.data_path('negros_panay_new.cfg')

    def test_multiprocessing(self):
        jobs = []
        for i in range(10):
            j = workers.DivModelSimulator(
                    config = self.new_cfg_path,
                    num_samples = 2,
                    )
            jobs.append(j)
        jobs = manager.Manager.run_workers(
                jobs,
                num_processors = 4)
        self.assertEqual(len(jobs), 10)
        for j in jobs:
            self.assertTrue(j.finished)
            self.assertIsInstance(j, workers.DivModelSimulator)
            self.assertIsInstance(j.config, MsBayesConfig)
            self.assertEqual(j.npairs, 9)
            self.assertEqual(j.num_samples, 2)
            self.assertIsInstance(j.div_models, stats.PartitionCollection)


if __name__ == '__main__':
    unittest.main()
