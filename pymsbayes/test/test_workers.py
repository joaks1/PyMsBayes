#! /usr/bin/env python

import unittest
import os
import sys
import shutil

from pymsbayes import workers
from pymsbayes.fileio import is_gzipped
from pymsbayes.test.support import package_paths
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.test import TestLevel
from pymsbayes.utils import ToolPathManager, MSBAYES_SORT_INDEX
from pymsbayes.utils.errors import *
from pymsbayes.utils.parsing import *
from pymsbayes.utils.stats import SampleSummaryCollection
from pymsbayes.utils.functions import whereis
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)


class ObsSumStatsWorkerTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('negros_panay.cfg')
        self.new_cfg_path = package_paths.data_path('negros_panay_new.cfg')
        self.expected = package_paths.data_path('negros_panay_sum_stats.txt')

    def tearDown(self):
        self.tear_down()

    def test_negros_panay(self):
        MSBAYES_SORT_INDEX.set_index(7)
        ss_path = self.get_test_path(prefix='sum-stats')
        ss_worker = workers.ObsSumStatsWorker(
                temp_fs = self.temp_fs,
                config_path = self.cfg_path,
                output_path = ss_path,
                exe_path = None,
                schema = 'abctoolbox',
                stat_patterns = DEFAULT_STAT_PATTERNS,
                stderr_path = None,
                tag = None)
        self.assertFalse(ss_worker.finished)
        ss_worker.start()
        self.assertTrue(ss_worker.finished)
        self.assertSameFiles([ss_path, self.expected])
        MSBAYES_SORT_INDEX.reset_default()

    def test_negros_panay_sort0_and_sort11(self):
        # for single locus data, sort 0 and 11 should be the same
        MSBAYES_SORT_INDEX.set_index(0)
        expected = package_paths.data_path('negros_panay_sum_stats_sort0.txt')
        expected2 = package_paths.data_path('negros_panay_sum_stats_sort11.txt')
        ss_path = self.get_test_path(prefix='sum-stats')
        ss_worker = workers.ObsSumStatsWorker(
                temp_fs = self.temp_fs,
                config_path = self.cfg_path,
                output_path = ss_path,
                exe_path = None,
                schema = 'abctoolbox',
                stat_patterns = DEFAULT_STAT_PATTERNS,
                stderr_path = None,
                tag = None)
        self.assertFalse(ss_worker.finished)
        ss_worker.start()
        self.assertTrue(ss_worker.finished)
        self.assertSameFiles([ss_path, expected])
        self.assertSameFiles([ss_path, expected2])
        MSBAYES_SORT_INDEX.reset_default()


class MsBayesWorkerTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')
        self.new_cfg_path = package_paths.data_path('4pairs_1locus_new.cfg')

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

    def test_purge(self):
        w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'abctoolbox')
        self.assertIsInstance(w, workers.MsBayesWorker)
        self.assertFalse(w.finished)
        w.start()
        self._assert_success(w, 4, 10)
        self.assertTrue(os.path.exists(w.output_dir))
        self.assertTrue(w.output_dir in self.temp_fs.dirs)
        w.purge()
        self.assertFalse(os.path.exists(w.output_dir))
        self.assertFalse(w.output_dir in self.temp_fs.dirs)

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

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
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

    def test_new_implementation(self):
        w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.new_cfg_path,
                schema = 'abctoolbox')
        self.assertIsInstance(w, workers.MsBayesWorker)
        self.assertFalse(w.finished)
        self.assertEqual(w.exe_path,
                ToolPathManager.get_tool_path('dpp-msbayes.pl'))
        self.assertNotEqual(w.exe_path,
                ToolPathManager.get_tool_path('msbayes.pl'))
        _LOG.warning('\n\n{0}\n\n'.format(w.exe_path))
        w.start()
        self._assert_success(w, 4, 10)
        self.assertEqual(w.exe_path,
                ToolPathManager.get_tool_path('dpp-msbayes.pl'))
        self.assertNotEqual(w.exe_path,
                ToolPathManager.get_tool_path('msbayes.pl'))


class MergePriorTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')

    def tearDown(self):
        self.tear_down()

    def _has_header(self, path):
        first_line = open(path, 'rU').next()
        if HEADER_PATTERN.match(first_line):
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

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
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
        ppath = self.get_test_path(prefix='merged_prior')
        hpath = self.get_test_path(prefix='merged_header')
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

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
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
        ppath = self.get_test_path(prefix='merged_prior')
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

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
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
        ppath = self.get_test_path(prefix='merged_prior')
        hpath = self.get_test_path(prefix='merged_header')
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

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
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
        ppath = self.get_test_path(prefix='merged_prior')
        hpath = self.get_test_path(prefix='merged_header')
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

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
    def test_merge_error_with_headers(self):
        alt_cfg_path = package_paths.data_path('negros_panay.cfg')
        w1 = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                include_header = True)
        w1.start()
        w2 = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = alt_cfg_path,
                schema = 'abctoolbox',
                include_header = True)
        w2.start()
        ppath = self.get_test_path(prefix='merged_prior')
        self.assertRaises(PriorMergeError, workers.merge_prior_files,
                [w1.prior_path, w2.prior_path], ppath, True)

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
    def test_merge_error_without_headers(self):
        alt_cfg_path = package_paths.data_path('negros_panay.cfg')
        w1 = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'msreject',
                include_header = False)
        w1.start()
        w2 = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = alt_cfg_path,
                schema = 'msreject',
                include_header = False)
        w2.start()
        ppath = self.get_test_path(prefix='merged_prior')
        self.assertRaises(PriorMergeError, workers.merge_prior_files,
                [w1.prior_path, w2.prior_path], ppath, True)

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
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
        ppath = self.get_test_path(prefix='merged_prior')
        workers.merge_prior_files(
                paths = [w.prior_path for w in jobs],
                dest_path = ppath)
        self.assertTrue(os.path.exists(ppath))
        self.assertEqual(self.get_number_of_header_lines(ppath), 1)
        self.assertEqual(self.get_number_of_lines(ppath), 41)

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
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
        ppath = self.get_test_path(prefix='merged_prior')
        workers.merge_prior_files(
                paths = [w.prior_path for w in jobs],
                dest_path = ppath)
        self.assertTrue(os.path.exists(ppath))
        self.assertEqual(self.get_number_of_header_lines(ppath), 0)
        self.assertEqual(self.get_number_of_lines(ppath), 40)

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
    def test_append_with_headers(self):
        jobs = []
        for i in range(4):
            w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                include_header = True)
            jobs.append(w)
        for w in jobs:
            w.start()
        ppath1 = self.get_test_path(prefix='merged_prior')
        ppath2 = self.get_test_path(prefix='merged_prior_append')
        ppath3 = self.get_test_path(prefix='merged_prior_append_last')
        ppath4 = self.get_test_path(prefix='merged_prior_append_each')
        workers.merge_prior_files(
                paths = [w.prior_path for w in jobs],
                dest_path = ppath1,
                append=False)
        self.assertTrue(os.path.exists(ppath1))
        self.assertEqual(self.get_number_of_header_lines(ppath1), 1)
        self.assertEqual(self.get_number_of_lines(ppath1), 41)
        workers.merge_prior_files(
                paths = [w.prior_path for w in jobs],
                dest_path = ppath2,
                append=True)
        shutil.copy(jobs[0].prior_path, ppath3)
        for i in range(1, len(jobs)):
            workers.merge_prior_files(
                    paths = [jobs[i].prior_path],
                    dest_path = ppath3,
                    append = True)
        for w in jobs:
            workers.merge_prior_files(
                    paths = [w.prior_path],
                    dest_path = ppath4,
                    append = True)
        self.assertSameFiles([ppath1, ppath2, ppath3, ppath4])

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
    def test_append_without_headers(self):
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
        ppath1 = self.get_test_path(prefix='merged_prior')
        ppath2 = self.get_test_path(prefix='merged_prior_append')
        ppath3 = self.get_test_path(prefix='merged_prior_append_last')
        ppath4 = self.get_test_path(prefix='merged_prior_append_each')
        workers.merge_prior_files(
                paths = [w.prior_path for w in jobs],
                dest_path = ppath1,
                append=False)
        self.assertTrue(os.path.exists(ppath1))
        self.assertEqual(self.get_number_of_header_lines(ppath1), 0)
        self.assertEqual(self.get_number_of_lines(ppath1), 40)
        workers.merge_prior_files(
                paths = [w.prior_path for w in jobs],
                dest_path = ppath2,
                append=True)
        shutil.copy(jobs[0].prior_path, ppath3)
        for i in range(1, len(jobs)):
            workers.merge_prior_files(
                    paths = [jobs[i].prior_path],
                    dest_path = ppath3,
                    append = True)
        for w in jobs:
            workers.merge_prior_files(
                    paths = [w.prior_path],
                    dest_path = ppath4,
                    append = True)
        self.assertSameFiles([ppath1, ppath2, ppath3, ppath4])

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
    def test_append_error_with_headers(self):
        alt_cfg_path = package_paths.data_path('negros_panay.cfg')
        w1 = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                include_header = True)
        w1.start()
        w2 = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = alt_cfg_path,
                schema = 'abctoolbox',
                include_header = True)
        w2.start()
        ppath = self.get_test_path(prefix='merged_prior')
        shutil.copy(w1.prior_path, ppath)
        self.assertRaises(PriorMergeError, workers.merge_prior_files,
                [w2.prior_path], ppath, True)

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
    def test_append_error_without_headers(self):
        alt_cfg_path = package_paths.data_path('negros_panay.cfg')
        w1 = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'msreject',
                include_header = False)
        w1.start()
        w2 = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = alt_cfg_path,
                schema = 'msreject',
                include_header = False)
        w2.start()
        ppath = self.get_test_path(prefix='merged_prior')
        shutil.copy(w1.prior_path, ppath)
        self.assertRaises(PriorMergeError, workers.merge_prior_files,
                [w2.prior_path], ppath, True)

    def test_append_with_compression(self):
        jobs = []
        for i in range(4):
            w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 10,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                include_header = True)
            jobs.append(w)
        for w in jobs:
            w.start()
        ppath = self.get_test_path(prefix='merged_prior')
        for w in jobs:
            workers.merge_prior_files(
                    paths = [w.prior_path],
                    dest_path = ppath,
                    append = True,
                    compresslevel = 9)
        self.assertTrue(os.path.exists(ppath))
        self.assertTrue(is_gzipped(ppath))
        self.assertEqual(self.get_number_of_header_lines(ppath), 1)
        self.assertEqual(self.get_number_of_lines(ppath), 41)


class MsRejectWorkerTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')

    def tearDown(self):
        self.tear_down()

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
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
        self.observed_path = self.get_test_path(prefix='test-obs')
        self.posterior_path = self.get_test_path(prefix='test-post')
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
        self.assertEqual(self.get_number_of_lines(w.posterior_path), 2)
        post = open(self.posterior_path, 'rU')
        post_header = post.next().strip().split()
        posterior = post.next().strip().split()
        post.close()
        self.assertEqual(self.msbayes_worker.header, w.header)
        pr = open(self.prior_path, 'rU')
        expected_post = pr.next().strip().split()
        pr.close()
        self.assertEqual(posterior, expected_post)

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
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
        self.observed_path = self.get_test_path(prefix='test-obs')
        self.posterior_path = self.get_test_path(prefix='test-post')
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
        self.assertEqual(self.get_number_of_lines(w.posterior_path), 2)
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
        self.observed_path = self.get_test_path(prefix='test-obs')
        self.posterior_path = self.get_test_path(prefix='test-post')
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
        self.new_prior_path = self.get_test_path(prefix='test-prior')
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
        self.assertEqual(self.get_number_of_lines(w.posterior_path), 3)
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

class EuRejectWorkerTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')

    def tearDown(self):
        self.tear_down()

    def test_summary_and_rejection(self):
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
        sum_out_path = self.get_test_path(prefix='test-summary-out-')
        reject_worker = workers.EuRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_paths = [prior_worker.prior_path],
                num_posterior_samples = 10,
                num_standardizing_samples = 100,
                summary_in_path = None,
                summary_out_path = sum_out_path,
                posterior_path = post_path,
                regression_worker = None,
                exe_path = None,
                stderr_path = None,
                keep_temps = False,
                tag = 'testcase')
        self.assertFalse(reject_worker.finished)
        reject_worker.start()
        self.assertTrue(reject_worker.finished)
        self.assertTrue(os.path.isfile(reject_worker.posterior_path))
        self.assertEqual(self.get_number_of_lines(
                reject_worker.posterior_path), 11)
        self.assertEqual(self.get_number_of_header_lines(
                reject_worker.posterior_path), 1)
        self.assertTrue(os.path.isfile(reject_worker.summary_out_path))
        self.assertEqual(self.get_number_of_lines(
                reject_worker.summary_out_path), 4)
        self.assertEqual(self.get_number_of_header_lines(
                reject_worker.summary_out_path), 1)
        self.assertEqual(reject_worker.num_processed, 100)
        self.assertEqual(reject_worker.num_summarized, 100)
        self.assertEqual(reject_worker.num_retained, 10)
        self.assertEqual(reject_worker.prior_paths,
                reject_worker.rejection_files)
        self.assertEqual(reject_worker.prior_paths,
                reject_worker.standardizing_files)
        self.assertFalse(os.path.exists(reject_worker.stderr_path))
        self.assertFalse(os.path.exists(reject_worker.output_dir))

    def test_summary_only(self):
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
        sum_out_path = self.get_test_path(prefix='test-summary-out-')
        reject_worker = workers.EuRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_paths = [prior_worker.prior_path],
                num_posterior_samples = 0,
                num_standardizing_samples = 100,
                summary_in_path = None,
                summary_out_path = sum_out_path,
                posterior_path = post_path,
                regression_worker = None,
                exe_path = None,
                stderr_path = None,
                keep_temps = False,
                tag = 'testcase')
        self.assertFalse(reject_worker.finished)
        reject_worker.start()
        self.assertTrue(reject_worker.finished)
        self.assertTrue(os.path.isfile(reject_worker.posterior_path))
        self.assertEqual(self.get_number_of_lines(
                reject_worker.posterior_path), 0)
        self.assertTrue(os.path.isfile(reject_worker.summary_out_path))
        self.assertEqual(self.get_number_of_lines(
                reject_worker.summary_out_path), 4)
        self.assertEqual(self.get_number_of_header_lines(
                reject_worker.summary_out_path), 1)
        self.assertEqual(reject_worker.num_processed, 0)
        self.assertEqual(reject_worker.num_summarized, 100)
        self.assertEqual(reject_worker.num_retained, 0)
        self.assertEqual(reject_worker.rejection_files, ['None'])
        self.assertEqual(reject_worker.prior_paths,
                reject_worker.standardizing_files)

    def test_rejection_only(self):
        prior_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 50,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                report_parameters = True)
        prior_worker.start()
        prior_worker2 = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 50,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                report_parameters = True)
        prior_worker2.start()
        obs_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True,
                report_parameters = True)
        obs_worker.start()

        post_path = self.get_test_path(prefix='test-posterior-')
        sum_out_path = self.get_test_path(prefix='test-summary-out-')
        sum_rej_worker = workers.EuRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_paths = [prior_worker.prior_path],
                num_posterior_samples = 10,
                num_standardizing_samples = 40,
                summary_in_path = None,
                summary_out_path = sum_out_path,
                posterior_path = post_path,
                regression_worker = None,
                exe_path = None,
                stderr_path = None,
                keep_temps = False,
                tag = 'testcase')
        self.assertFalse(sum_rej_worker.finished)
        sum_rej_worker.start()
        self.assertTrue(sum_rej_worker.finished)
        self.assertTrue(os.path.isfile(sum_rej_worker.posterior_path))
        self.assertEqual(self.get_number_of_lines(
                sum_rej_worker.posterior_path), 11)
        self.assertEqual(self.get_number_of_header_lines(
                sum_rej_worker.posterior_path), 1)
        self.assertTrue(os.path.isfile(sum_rej_worker.summary_out_path))
        self.assertEqual(self.get_number_of_lines(
                sum_rej_worker.summary_out_path), 4)
        self.assertEqual(self.get_number_of_header_lines(
                sum_rej_worker.summary_out_path), 1)
        self.assertEqual(sum_rej_worker.num_processed, 50)
        self.assertEqual(sum_rej_worker.num_summarized, 40)
        self.assertEqual(sum_rej_worker.num_retained, 10)
        self.assertEqual(sum_rej_worker.prior_paths,
                sum_rej_worker.rejection_files)
        self.assertEqual(sum_rej_worker.prior_paths,
                sum_rej_worker.standardizing_files)

        post_path2 = self.get_test_path(prefix='test-posterior-')
        sum_out_path2 = self.get_test_path(prefix='test-summary-out-')
        reject_worker = workers.EuRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_paths = [prior_worker2.prior_path],
                num_posterior_samples = 5,
                num_standardizing_samples = 0,
                summary_in_path = sum_rej_worker.summary_out_path,
                summary_out_path = sum_out_path2,
                posterior_path = post_path2,
                regression_worker = None,
                exe_path = None,
                stderr_path = None,
                keep_temps = False,
                tag = 'testcase')
        self.assertFalse(reject_worker.finished)
        reject_worker.start()
        self.assertTrue(reject_worker.finished)
        self.assertTrue(os.path.isfile(reject_worker.posterior_path))
        self.assertEqual(self.get_number_of_header_lines(
                reject_worker.posterior_path), 1)
        self.assertEqual(self.get_number_of_lines(
                reject_worker.posterior_path), 6)
        self.assertTrue(os.path.isfile(reject_worker.summary_out_path))
        self.assertEqual(self.get_number_of_lines(
                reject_worker.summary_out_path), 4)
        self.assertEqual(self.get_number_of_header_lines(
                reject_worker.summary_out_path), 1)
        self.assertEqual(reject_worker.num_processed, 50)
        self.assertEqual(reject_worker.num_summarized, 0)
        self.assertEqual(reject_worker.num_retained, 5)
        self.assertEqual(reject_worker.prior_paths,
                reject_worker.rejection_files)
        self.assertEqual(reject_worker.standardizing_files, ['None'])

        self.assertSameFiles([reject_worker.summary_out_path,
            sum_rej_worker.summary_out_path])

        post_path3 = self.get_test_path(prefix='test-posterior-')
        sum_out_path3 = self.get_test_path(prefix='test-summary-out-')
        reject_worker2 = workers.EuRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_paths = [prior_worker.prior_path],
                num_posterior_samples = 10,
                num_standardizing_samples = 0,
                summary_in_path = reject_worker.summary_out_path,
                summary_out_path = sum_out_path3,
                posterior_path = post_path3,
                regression_worker = None,
                exe_path = None,
                stderr_path = None,
                keep_temps = False,
                tag = 'testcase')
        self.assertFalse(reject_worker2.finished)
        reject_worker2.start()
        self.assertTrue(reject_worker2.finished)
        self.assertTrue(os.path.isfile(reject_worker2.posterior_path))
        self.assertEqual(self.get_number_of_header_lines(
                reject_worker2.posterior_path), 1)
        self.assertEqual(self.get_number_of_lines(
                reject_worker2.posterior_path), 11)
        self.assertTrue(os.path.isfile(reject_worker2.summary_out_path))
        self.assertEqual(self.get_number_of_lines(
                reject_worker2.summary_out_path), 4)
        self.assertEqual(self.get_number_of_header_lines(
                reject_worker2.summary_out_path), 1)
        self.assertEqual(reject_worker2.num_processed, 50)
        self.assertEqual(reject_worker2.num_summarized, 0)
        self.assertEqual(reject_worker2.num_retained, 10)
        self.assertEqual(reject_worker2.prior_paths,
                reject_worker2.rejection_files)
        self.assertEqual(reject_worker2.standardizing_files, ['None'])

        self.assertSameFiles([reject_worker2.posterior_path,
            sum_rej_worker.posterior_path])
        self.assertSameFiles([reject_worker2.summary_out_path,
            sum_rej_worker.summary_out_path])

class EuRejectSummaryMergerTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')

    def tearDown(self):
        self.tear_down()

    def test_num_summarize_exception(self):
        prior_workers = []
        for i in range(2):
            prior_worker = workers.MsBayesWorker(
                    temp_fs = self.temp_fs,
                    sample_size = 10,
                    config_path = self.cfg_path,
                    schema = 'abctoolbox',
                    report_parameters = True)
            prior_worker.start()
            prior_workers.append(prior_worker)
        obs_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True,
                report_parameters = True)
        obs_worker.start()

        reject_workers = []
        for pw in prior_workers:
            post_path = self.get_test_path(prefix='test-posterior-')
            sum_out_path = self.get_test_path(prefix='test-summary-out-')
            reject_worker = workers.EuRejectWorker(
                    temp_fs = self.temp_fs,
                    observed_path = obs_worker.prior_stats_path,
                    prior_paths = [pw.prior_path],
                    num_posterior_samples = 0,
                    num_standardizing_samples = 11,
                    summary_in_path = None,
                    summary_out_path = sum_out_path,
                    posterior_path = post_path,
                    regression_worker = None,
                    exe_path = None,
                    stderr_path = None,
                    keep_temps = False,
                    tag = 'testcase')
            reject_worker.start()
            reject_workers.append(reject_worker)

        sum_out_path = self.get_test_path(prefix='test-summary-merged-')
        merger = workers.EuRejectSummaryMerger(reject_workers)
        self.assertRaises(Exception, merger.start)

    def test_merge_equal_sample_sizes(self):
        prior_workers = []
        for i in range(3):
            prior_worker = workers.MsBayesWorker(
                    temp_fs = self.temp_fs,
                    sample_size = 10,
                    config_path = self.cfg_path,
                    schema = 'abctoolbox',
                    report_parameters = True)
            prior_worker.start()
            prior_workers.append(prior_worker)
        obs_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True,
                report_parameters = True)
        obs_worker.start()

        reject_workers = []
        for pw in prior_workers:
            post_path = self.get_test_path(prefix='test-posterior-')
            sum_out_path = self.get_test_path(prefix='test-summary-out-')
            reject_worker = workers.EuRejectWorker(
                    temp_fs = self.temp_fs,
                    observed_path = obs_worker.prior_stats_path,
                    prior_paths = [pw.prior_path],
                    num_posterior_samples = 0,
                    num_standardizing_samples = 10,
                    summary_in_path = None,
                    summary_out_path = sum_out_path,
                    posterior_path = post_path,
                    regression_worker = None,
                    exe_path = None,
                    stderr_path = None,
                    keep_temps = False,
                    tag = 'testcase')
            reject_worker.start()
            reject_workers.append(reject_worker)

        post_path = self.get_test_path(prefix='test-posterior-')
        sum_out_path = self.get_test_path(prefix='test-summary-out-')
        reject_worker = workers.EuRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_paths = [pw.prior_path for pw in prior_workers],
                num_posterior_samples = 0,
                num_standardizing_samples = 30,
                summary_in_path = None,
                summary_out_path = sum_out_path,
                posterior_path = post_path,
                regression_worker = None,
                exe_path = None,
                stderr_path = None,
                keep_temps = False,
                tag = 'testcase')
        reject_worker.start()

        sum_out_path = self.get_test_path(prefix='test-summary-merged-')
        merger = workers.EuRejectSummaryMerger(reject_workers)
        self.assertFalse(merger.finished)
        merger.start()
        self.assertTrue(merger.finished)
        merger.write_summary(sum_out_path)
        self.assertEqual(self.get_number_of_lines(
                sum_out_path), 4)
        self.assertEqual(self.get_number_of_header_lines(
                sum_out_path), 1)
        
        result, header = parse_summary_file(sum_out_path)
        expected_result, expected_header = parse_summary_file(
                reject_worker.summary_out_path)
        self.assertEqual(header, expected_header)
        self.assertEqual(sorted(result.keys()), sorted(expected_result.keys()))
        for stat_name in result.keys():
            self.assertAlmostEqual(result[stat_name]['mean'],
                    expected_result[stat_name]['mean'])
            self.assertAlmostEqual(result[stat_name]['std_deviation'],
                    expected_result[stat_name]['std_deviation'])
            self.assertEqual(result[stat_name]['n'],
                    expected_result[stat_name]['n'])

    def test_merge_unequal_sample_sizes(self):
        prior_workers = []
        sizes = [20, 15, 7]
        for i in sizes:
            prior_worker = workers.MsBayesWorker(
                    temp_fs = self.temp_fs,
                    sample_size = i,
                    config_path = self.cfg_path,
                    schema = 'abctoolbox',
                    report_parameters = True)
            prior_worker.start()
            prior_workers.append(prior_worker)
        obs_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                write_stats_file = True,
                report_parameters = True)
        obs_worker.start()

        reject_workers = []
        for pw in prior_workers:
            post_path = self.get_test_path(prefix='test-posterior-')
            sum_out_path = self.get_test_path(prefix='test-summary-out-')
            reject_worker = workers.EuRejectWorker(
                    temp_fs = self.temp_fs,
                    observed_path = obs_worker.prior_stats_path,
                    prior_paths = [pw.prior_path],
                    num_posterior_samples = 0,
                    num_standardizing_samples = pw.sample_size,
                    summary_in_path = None,
                    summary_out_path = sum_out_path,
                    posterior_path = post_path,
                    regression_worker = None,
                    exe_path = None,
                    stderr_path = None,
                    keep_temps = False,
                    tag = 'testcase')
            reject_worker.start()
            reject_workers.append(reject_worker)

        post_path = self.get_test_path(prefix='test-posterior-')
        sum_out_path = self.get_test_path(prefix='test-summary-out-')
        reject_worker = workers.EuRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_paths = [pw.prior_path for pw in prior_workers],
                num_posterior_samples = 0,
                num_standardizing_samples = sum(sizes),
                summary_in_path = None,
                summary_out_path = sum_out_path,
                posterior_path = post_path,
                regression_worker = None,
                exe_path = None,
                stderr_path = None,
                keep_temps = False,
                tag = 'testcase')
        reject_worker.start()

        sum_out_path = self.get_test_path(prefix='test-summary-merged-')
        merger = workers.EuRejectSummaryMerger(reject_workers)
        self.assertFalse(merger.finished)
        merger.start()
        self.assertTrue(merger.finished)
        merger.write_summary(sum_out_path)
        self.assertEqual(self.get_number_of_lines(
                sum_out_path), 4)
        self.assertEqual(self.get_number_of_header_lines(
                sum_out_path), 1)
        
        result, header = parse_summary_file(sum_out_path)
        expected_result, expected_header = parse_summary_file(
                reject_worker.summary_out_path)
        self.assertEqual(header, expected_header)
        self.assertEqual(sorted(result.keys()), sorted(expected_result.keys()))
        for stat_name in result.keys():
            self.assertAlmostEqual(result[stat_name]['mean'],
                    expected_result[stat_name]['mean'])
            self.assertAlmostEqual(result[stat_name]['std_deviation'],
                    expected_result[stat_name]['std_deviation'])
            self.assertEqual(result[stat_name]['n'],
                    expected_result[stat_name]['n'])


class ABCToolBoxRejectWorkerTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')

    def tearDown(self):
        self.tear_down()

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
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

class RejectWorkerComparisonTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')

    def tearDown(self):
        self.tear_down()

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
    def test_compare_rejectors_100(self):
        prior_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 100,
                config_path = self.cfg_path,
                schema = 'msreject',
                include_header = False,
                report_parameters = True)
        prior_worker.start()
        obs_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'msreject',
                include_header = False,
                write_stats_file = True,
                report_parameters = True)
        obs_worker.start()
        headed_prior_path = self.get_test_path(prefix='headed-prior')
        with open(headed_prior_path, 'w') as out:
            out.write('{0}\n'.format('\t'.join(obs_worker.header)))
            with open(prior_worker.prior_path, 'rU') as prior_file:
                for line in prior_file:
                    l = line.strip().split()
                    out.write('{0}\n'.format('\t'.join(l)))

        abc_post_path = self.get_test_path(prefix='abc-posterior-')
        ms_post_path = self.get_test_path(prefix='ms-posterior-')
        eu_post_path = self.get_test_path(prefix='eu-posterior-')
        eu_sum_out_path = self.get_test_path(prefix='eu-sum-out-')
        eu_reject_worker = workers.EuRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_paths = [headed_prior_path],
                num_posterior_samples = 10,
                num_standardizing_samples = 100,
                summary_in_path = None,
                summary_out_path = eu_sum_out_path,
                posterior_path = eu_post_path,
                regression_worker = None,
                exe_path = None,
                stderr_path = None,
                keep_temps = False,
                tag = '')
        eu_reject_worker.start()
        abc_reject_worker = workers.ABCToolBoxRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_path = headed_prior_path,
                num_posterior_samples = 10,
                posterior_path = abc_post_path,
                regression_worker = None,
                exe_path = None,
                stdout_path = None,
                stderr_path = None,
                keep_temps = False,
                max_read_sims = 1000)
        abc_reject_worker.start()
        ms_reject_worker = workers.MsRejectWorker(
                header = prior_worker.header,
                observed_path = obs_worker.prior_path,
                prior_path = prior_worker.prior_path,
                tolerance = 0.1,
                posterior_path = ms_post_path)
        ms_reject_worker.start()
        self.assertSameSamples(files = [abc_reject_worker.posterior_path,
                    ms_reject_worker.posterior_path,
                    eu_reject_worker.posterior_path],
                columns_to_ignore = [0],
                header = True,
                places = 4,
                num_mismatches_per_sample = 1,
                num_sample_mismatches = 1)
        self.assertSameSamples(files = [
                    ms_reject_worker.posterior_path,
                    eu_reject_worker.posterior_path],
                columns_to_ignore = [0],
                header = True,
                places = 5,
                num_mismatches_per_sample = 0,
                num_sample_mismatches = 0)

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
    def test_compare_rejectors_1000(self):
        prior_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1000,
                config_path = self.cfg_path,
                schema = 'msreject',
                include_header = False,
                report_parameters = True)
        prior_worker.start()
        obs_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'msreject',
                include_header = False,
                write_stats_file = True,
                report_parameters = True)
        obs_worker.start()
        headed_prior_path = self.get_test_path(prefix='headed-prior')
        with open(headed_prior_path, 'w') as out:
            out.write('{0}\n'.format('\t'.join(obs_worker.header)))
            with open(prior_worker.prior_path, 'rU') as prior_file:
                for line in prior_file:
                    l = line.strip().split()
                    out.write('{0}\n'.format('\t'.join(l)))

        abc_post_path = self.get_test_path(prefix='abc-posterior-')
        ms_post_path = self.get_test_path(prefix='ms-posterior-')
        eu_post_path = self.get_test_path(prefix='eu-posterior-')
        eu_sum_out_path = self.get_test_path(prefix='eu-sum-out-')
        eu_reject_worker = workers.EuRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_paths = [headed_prior_path],
                num_posterior_samples = 100,
                num_standardizing_samples = 1000,
                summary_in_path = None,
                summary_out_path = eu_sum_out_path,
                posterior_path = eu_post_path,
                regression_worker = None,
                exe_path = None,
                stderr_path = None,
                keep_temps = False,
                tag = '')
        eu_reject_worker.start()
        abc_reject_worker = workers.ABCToolBoxRejectWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                prior_path = headed_prior_path,
                num_posterior_samples = 100,
                posterior_path = abc_post_path,
                regression_worker = None,
                exe_path = None,
                stdout_path = None,
                stderr_path = None,
                keep_temps = False,
                max_read_sims = 1000)
        abc_reject_worker.start()
        ms_reject_worker = workers.MsRejectWorker(
                header = prior_worker.header,
                observed_path = obs_worker.prior_path,
                prior_path = prior_worker.prior_path,
                tolerance = 0.1,
                posterior_path = ms_post_path)
        ms_reject_worker.start()
        self.assertSameSamples(files = [abc_reject_worker.posterior_path,
                    ms_reject_worker.posterior_path,
                    eu_reject_worker.posterior_path],
                columns_to_ignore = [0],
                header = True,
                places = 4,
                num_mismatches_per_sample = 1,
                num_sample_mismatches = 1)
        self.assertSameSamples(files = [
                    ms_reject_worker.posterior_path,
                    eu_reject_worker.posterior_path],
                columns_to_ignore = [0],
                header = True,
                places = 5,
                num_mismatches_per_sample = 0,
                num_sample_mismatches = 0)

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
        regress_posterior_path = self.get_test_path(prefix='test-adjusted-')
        regress_worker = workers.ABCToolBoxRegressWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                posterior_path = post_worker.prior_path,
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
        regress_worker.start()
        self.assertTrue(regress_worker.finished)
        self.assertTrue(os.path.isfile(regress_worker.regress_summary_path))
        self.assertTrue(os.path.isfile(regress_worker.regress_posterior_path))
        self.assertEqual(self.get_number_of_lines(
                regress_worker.regress_posterior_path), 101)

    def test_regression_failure(self):
        post_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
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
        regress_posterior_path = self.get_test_path(prefix='test-adjusted-')
        regress_worker = workers.ABCToolBoxRegressWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                posterior_path = post_worker.prior_path,
                parameter_indices = None,
                regress_summary_path = summary_path,
                regress_posterior_path = regress_posterior_path,
                exe_path = whereis('echo'),
                stdout_path = None,
                stderr_path = None,
                keep_temps = False,
                bandwidth = None,
                num_posterior_quantiles = 100)
        self.assertFalse(regress_worker.finished)
        regress_worker.start()
        self.assertTrue(regress_worker.finished)
        self.assertTrue(regress_worker.failed)
        self.assertEqual(regress_worker.regress_summary_path, None)
        self.assertEqual(regress_worker.regress_posterior_path, None)

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
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
        regress_posterior_path = self.get_test_path(prefix='test-adjusted-')
        regress_worker = workers.ABCToolBoxRegressWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
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
        regress_worker.start()
        self.assertTrue(regress_worker.finished)
        self.assertTrue(os.path.isfile(regress_worker.regress_summary_path))
        self.assertTrue(os.path.isfile(regress_worker.regress_posterior_path))
        self.assertEqual(self.get_number_of_lines(
                regress_worker.regress_posterior_path), 101)

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
        regress_posterior_path = self.get_test_path(prefix='test-adjusted-')
        regress_worker = workers.ABCToolBoxRegressWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                posterior_path = post_path,
                parameter_indices = None,
                regress_summary_path = summary_path,
                regress_posterior_path = regress_posterior_path,
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
        self.assertTrue(os.path.isfile(regress_worker.regress_summary_path))
        self.assertTrue(os.path.isfile(regress_worker.regress_posterior_path))
        self.assertEqual(self.get_number_of_lines(
                regress_worker.regress_posterior_path), 101)


class AssembleMsRejectWorkersTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')
        self.results_dir = self.get_test_subdir(
                prefix='test-posteriors')
        self.posterior_prefix = self.test_id + '-posterior'
        self.expected_continuous_parameters = set(['PRI.omega', 'PRI.E.t',
                'PRI.cv'])
        self.expected_discrete_parameters = set(['PRI.Psi'])
        self.expected_discrete_parameters_with_models = set(['PRI.Psi',
                'PRI.model'])
        self.expected_parameters = set(
                list(self.expected_continuous_parameters) + \
                list(self.expected_discrete_parameters))
        self.expected_parameters_with_models = set(
                list(self.expected_continuous_parameters) + \
                list(self.expected_discrete_parameters_with_models))

    def tearDown(self):
        self.tear_down()

    def _get_prior(self, n=100, schema='msreject',
            model_index = None,
            include_header=False):
        w = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = n,
                config_path = self.cfg_path,
                schema = schema,
                report_parameters = True,
                model_index = model_index,
                include_header = include_header)
        w.start()
        return w
    
    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
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
            self.assertEqual(self.get_number_of_lines(j.posterior_path), 6)

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
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

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
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
            self.assertEqual(self.get_number_of_lines(j.posterior_path), 101)
            self.assertTrue(os.path.isfile(
                    j.regression_worker.regress_summary_path))
            self.assertTrue(os.path.isfile(
                    j.regression_worker.regress_posterior_path))
            self.assertEqual(self.get_number_of_lines(
                    j.regression_worker.regress_posterior_path), 101)
        # check regression results
        parameters = [prior_worker.header[
                i] for i in prior_worker.parameter_indices]
        stats = [prior_worker.header[
                i] for i in prior_worker.stat_indices]
        for j in msreject_workers:
            cont_params = [j.regression_worker.header[
                    i] for i in j.regression_worker.continuous_parameter_indices]
            disc_params = [j.regression_worker.header[
                    i] for i in j.regression_worker.discrete_parameter_indices]
            self.assertEqual(set(cont_params),
                    self.expected_continuous_parameters)
            self.assertEqual(set(disc_params),
                    self.expected_discrete_parameters)
            rej_stats = [j.header[i] for i in j.stat_indices]
            self.assertEqual(set(cont_params + disc_params) - set(parameters), set())
            self.assertEqual(sorted(stats), sorted(rej_stats))
            result_header = open(j.regression_worker.regress_posterior_path, 
                    'rU').next().strip().split()
            self.assertEqual(sorted(result_header), sorted(cont_params))
            results = self.parse_python_config(
                    j.regression_worker.regress_summary_path)
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

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
    def test_msreject_llr_with_models(self):
        prior_workers = []
        for i in range(2):
            prior_workers.append(self._get_prior(n=100, model_index=i))
        prior_path = self.get_test_path()
        header_path = self.get_test_path()
        workers.merge_priors(prior_workers,
                prior_path = prior_path,
                header_path = header_path,
                include_header = False)
        prior_worker = prior_workers[0]
        prior_worker.prior_path = prior_path
        obs_worker = self._get_prior(n=2, model_index=1,
                include_header=True)
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
            self.assertEqual(self.get_number_of_lines(j.posterior_path), 101)
            self.assertTrue(os.path.isfile(
                    j.regression_worker.regress_summary_path))
            self.assertTrue(os.path.isfile(
                    j.regression_worker.regress_posterior_path))
            self.assertEqual(self.get_number_of_lines(
                    j.regression_worker.regress_posterior_path), 101)
        # check regression results
        parameters = [prior_worker.header[
                i] for i in prior_worker.parameter_indices]
        stats = [prior_worker.header[
                i] for i in prior_worker.stat_indices]
        for j in msreject_workers:
            cont_params = [j.regression_worker.header[
                    i] for i in j.regression_worker.continuous_parameter_indices]
            disc_params = [j.regression_worker.header[
                    i] for i in j.regression_worker.discrete_parameter_indices]
            self.assertEqual(set(cont_params),
                    self.expected_continuous_parameters)
            self.assertEqual(set(disc_params),
                    self.expected_discrete_parameters_with_models)
            rej_stats = [j.header[i] for i in j.stat_indices]
            self.assertEqual(set(cont_params + disc_params) - set(parameters), set())
            self.assertEqual(sorted(stats), sorted(rej_stats))
            result_header = open(j.regression_worker.regress_posterior_path, 
                    'rU').next().strip().split()
            self.assertEqual(sorted(result_header), sorted(cont_params))
            results = self.parse_python_config(
                    j.regression_worker.regress_summary_path)
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

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
    def test_msreject_glm(self):
        prior_worker = self._get_prior(n=200)
        obs_worker = self._get_prior(n=2, include_header=True)
        msreject_workers = workers.assemble_rejection_workers(
                temp_fs = self.temp_fs,
                observed_sims_file = obs_worker.prior_path,
                prior_path = prior_worker.prior_path,
                num_prior_samples = 200,
                num_posterior_samples = 100,
                results_dir = self.results_dir,
                posterior_prefix = self.posterior_prefix,
                regress = True,
                rejection_tool = 'msreject',
                regression_method = 'glm',
                num_posterior_quantiles = 50)
        self.assertEqual(len(msreject_workers), 2)
        for j in msreject_workers:
            self.assertFalse(j.finished)
            self.assertFalse(j.regression_worker.finished)
            j.start()
        for j in msreject_workers:
            self.assertTrue(j.finished)
            self.assertTrue(j.regression_worker.finished)
            self.assertTrue(os.path.isfile(j.posterior_path))
            self.assertEqual(self.get_number_of_lines(j.posterior_path), 101)
            self.assertTrue(os.path.isfile(
                    j.regression_worker.regress_summary_path))
            self.assertTrue(os.path.isfile(
                    j.regression_worker.regress_posterior_path))
            self.assertEqual(self.get_number_of_lines(
                    j.regression_worker.regress_posterior_path), 51)
            self.assertEqual(self.get_number_of_lines(
                    j.regression_worker.regress_summary_path), 20)
        # check regression results
        parameters = [prior_worker.header[
                i] for i in prior_worker.parameter_indices]
        stats = [prior_worker.header[
                i] for i in prior_worker.stat_indices]
        for j in msreject_workers:
            reg_params = [j.regression_worker.header[
                    i] for i in j.regression_worker.parameter_indices]
            self.assertEqual(set(reg_params), self.expected_parameters)
            rej_stats = [j.header[i] for i in j.stat_indices]
            self.assertEqual(set(reg_params) - set(parameters), set())
            self.assertEqual(sorted(stats), sorted(rej_stats))
            result_header = set(open(j.regression_worker.regress_posterior_path, 
                    'rU').next().strip().split()) - \
                    set(['number', 'density'])
            expected_header = set(reg_params)
            self.assertEqual(result_header, expected_header)
            summary_header = set(open(j.regression_worker.regress_summary_path,
                    'rU').next().strip().split()) -  set(['what'])
            self.assertEqual(summary_header, expected_header)
            self.assertEqual(sorted(stats),
                    sorted(j.regression_worker.stats_header))

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
    def test_msreject_glm_with_models(self):
        prior_workers = []
        for i in range(2):
            prior_workers.append(self._get_prior(n=100, model_index=i))
        prior_path = self.get_test_path()
        header_path = self.get_test_path()
        workers.merge_priors(prior_workers,
                prior_path = prior_path,
                header_path = header_path,
                include_header = False)
        prior_worker = prior_workers[0]
        prior_worker.prior_path = prior_path
        obs_worker = self._get_prior(n=2, model_index=1,
                include_header=True)
        msreject_workers = workers.assemble_rejection_workers(
                temp_fs = self.temp_fs,
                observed_sims_file = obs_worker.prior_path,
                prior_path = prior_worker.prior_path,
                num_prior_samples = 200,
                num_posterior_samples = 100,
                results_dir = self.results_dir,
                posterior_prefix = self.posterior_prefix,
                regress = True,
                rejection_tool = 'msreject',
                regression_method = 'glm',
                num_posterior_quantiles = 50)
        self.assertEqual(len(msreject_workers), 2)
        for j in msreject_workers:
            self.assertFalse(j.finished)
            self.assertFalse(j.regression_worker.finished)
            j.start()
        for j in msreject_workers:
            self.assertTrue(j.finished)
            self.assertTrue(j.regression_worker.finished)
            self.assertTrue(os.path.isfile(j.posterior_path))
            self.assertEqual(self.get_number_of_lines(j.posterior_path), 101)
            self.assertTrue(os.path.isfile(
                    j.regression_worker.regress_summary_path))
            self.assertTrue(os.path.isfile(
                    j.regression_worker.regress_posterior_path))
            self.assertEqual(self.get_number_of_lines(
                    j.regression_worker.regress_posterior_path), 51)
            self.assertEqual(self.get_number_of_lines(
                    j.regression_worker.regress_summary_path), 20)
        # check regression results
        parameters = [prior_worker.header[
                i] for i in prior_worker.parameter_indices]
        stats = [prior_worker.header[
                i] for i in prior_worker.stat_indices]
        for j in msreject_workers:
            reg_params = [j.regression_worker.header[
                    i] for i in j.regression_worker.parameter_indices]
            self.assertEqual(set(reg_params),
                    self.expected_parameters_with_models)
            rej_stats = [j.header[i] for i in j.stat_indices]
            self.assertEqual(set(reg_params) - set(parameters), set())
            self.assertEqual(sorted(stats), sorted(rej_stats))
            result_header = set(open(j.regression_worker.regress_posterior_path, 
                    'rU').next().strip().split()) - \
                    set(['number', 'density'])
            expected_header = set(reg_params)
            self.assertEqual(result_header, expected_header)
            summary_header = set(open(j.regression_worker.regress_summary_path,
                    'rU').next().strip().split()) -  set(['what'])
            self.assertEqual(summary_header, expected_header)
            self.assertEqual(sorted(stats),
                    sorted(j.regression_worker.stats_header))

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
    def test_abctoolbox_glm(self):
        prior_worker = self._get_prior(n=200, schema='abctoolbox')
        obs_worker = self._get_prior(n=2, schema='abctoolbox')
        jobs = workers.assemble_rejection_workers(
                temp_fs = self.temp_fs,
                observed_sims_file = obs_worker.prior_path,
                prior_path = prior_worker.prior_path,
                num_prior_samples = 200,
                num_posterior_samples = 100,
                results_dir = self.results_dir,
                posterior_prefix = self.posterior_prefix,
                regress = True,
                rejection_tool = 'abctoolbox',
                regression_method = 'glm',
                num_posterior_quantiles = 50)
        self.assertEqual(len(jobs), 2)
        for j in jobs:
            self.assertFalse(j.finished)
            j.start()
        for j in jobs:
            self.assertTrue(j.finished)
            self.assertTrue(os.path.isfile(j.posterior_path))
            self.assertTrue(self.get_number_of_lines(j.posterior_path), 101)
            self.assertTrue(os.path.isfile(
                    j.regression_worker.regress_summary_path))
            self.assertTrue(os.path.isfile(
                    j.regression_worker.regress_posterior_path))
            self.assertEqual(self.get_number_of_lines(
                    j.regression_worker.regress_posterior_path), 51)
            self.assertEqual(self.get_number_of_lines(
                    j.regression_worker.regress_summary_path), 20)
        # check regression results
        parameters = [prior_worker.header[
                i] for i in prior_worker.parameter_indices]
        stats = [prior_worker.header[
                i] for i in prior_worker.stat_indices]
        for j in jobs:
            reg_params = [j.regression_worker.header[
                    i] for i in j.regression_worker.parameter_indices]
            self.assertEqual(set(reg_params), self.expected_parameters)
            rej_params = [j.header[i] for i in j.parameter_indices]
            self.assertEqual(set(reg_params) - set(parameters), set())
            self.assertEqual(sorted(parameters), sorted(rej_params))
            result_header = set(open(j.regression_worker.regress_posterior_path, 
                    'rU').next().strip().split()) - \
                    set(['number', 'density'])
            expected_header = set(reg_params)
            self.assertEqual(result_header, expected_header)
            summary_header = set(open(j.regression_worker.regress_summary_path,
                    'rU').next().strip().split()) -  set(['what'])
            self.assertEqual(summary_header, expected_header)
            self.assertEqual(sorted(stats), sorted(j.stats_header))
            self.assertEqual(sorted(stats),
                    sorted(j.regression_worker.stats_header))

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
    def test_abctoolbox_glm_with_models(self):
        prior_workers = []
        for i in range(2):
            prior_workers.append(self._get_prior(n=100, model_index=i,
                schema='abctoolbox'))
        prior_path = self.get_test_path()
        header_path = self.get_test_path()
        workers.merge_priors(prior_workers,
                prior_path = prior_path,
                header_path = header_path,
                include_header = True)
        prior_worker = prior_workers[0]
        prior_worker.prior_path = prior_path
        obs_worker = self._get_prior(n=2, model_index=1,
                schema='abctoolbox')
        jobs = workers.assemble_rejection_workers(
                temp_fs = self.temp_fs,
                observed_sims_file = obs_worker.prior_path,
                prior_path = prior_worker.prior_path,
                num_prior_samples = 200,
                num_posterior_samples = 100,
                results_dir = self.results_dir,
                posterior_prefix = self.posterior_prefix,
                regress = True,
                rejection_tool = 'abctoolbox',
                regression_method = 'glm',
                num_posterior_quantiles = 50)
        self.assertEqual(len(jobs), 2)
        for j in jobs:
            self.assertFalse(j.finished)
            j.start()
        for j in jobs:
            self.assertTrue(j.finished)
            self.assertTrue(os.path.isfile(j.posterior_path))
            self.assertTrue(self.get_number_of_lines(j.posterior_path), 101)
            self.assertTrue(os.path.isfile(
                    j.regression_worker.regress_summary_path))
            self.assertTrue(os.path.isfile(
                    j.regression_worker.regress_posterior_path))
            self.assertEqual(self.get_number_of_lines(
                    j.regression_worker.regress_posterior_path), 51)
            self.assertEqual(self.get_number_of_lines(
                    j.regression_worker.regress_summary_path), 20)
        # check regression results
        parameters = [prior_worker.header[
                i] for i in prior_worker.parameter_indices]
        stats = [prior_worker.header[
                i] for i in prior_worker.stat_indices]
        for j in jobs:
            reg_params = [j.regression_worker.header[
                    i] for i in j.regression_worker.parameter_indices]
            self.assertEqual(set(reg_params),
                    self.expected_parameters_with_models)
            rej_params = [j.header[i] for i in j.parameter_indices]
            self.assertEqual(set(reg_params) - set(parameters), set())
            self.assertEqual(sorted(parameters), sorted(rej_params))
            result_header = set(open(j.regression_worker.regress_posterior_path, 
                    'rU').next().strip().split()) - \
                    set(['number', 'density'])
            expected_header = set(reg_params)
            self.assertEqual(result_header, expected_header)
            summary_header = set(open(j.regression_worker.regress_summary_path,
                    'rU').next().strip().split()) -  set(['what'])
            self.assertEqual(summary_header, expected_header)
            self.assertEqual(sorted(stats), sorted(j.stats_header))
            self.assertEqual(sorted(stats),
                    sorted(j.regression_worker.stats_header))

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
    def test_abctoolbox_llr(self):
        prior_worker = self._get_prior(n=200, schema='abctoolbox')
        obs_worker = self._get_prior(n=2, schema='abctoolbox')
        jobs = workers.assemble_rejection_workers(
                temp_fs = self.temp_fs,
                observed_sims_file = obs_worker.prior_path,
                prior_path = prior_worker.prior_path,
                num_prior_samples = 200,
                num_posterior_samples = 100,
                results_dir = self.results_dir,
                posterior_prefix = self.posterior_prefix,
                regress = True,
                rejection_tool = 'abctoolbox',
                regression_method = 'llr')
        self.assertEqual(len(jobs), 2)
        for j in jobs:
            self.assertFalse(j.finished)
            j.start()
        for j in jobs:
            self.assertTrue(j.finished)
            self.assertTrue(os.path.isfile(j.posterior_path))
            self.assertTrue(self.get_number_of_lines(j.posterior_path), 101)
            self.assertTrue(os.path.isfile(
                    j.regression_worker.regress_summary_path))
            self.assertTrue(os.path.isfile(
                    j.regression_worker.regress_posterior_path))
            self.assertEqual(self.get_number_of_lines(
                    j.regression_worker.regress_posterior_path), 101)
            self.assertTrue(self.get_number_of_lines(
                    j.regression_worker.regress_summary_path) > 10)
        # check regression results
        parameters = [prior_worker.header[
                i] for i in prior_worker.parameter_indices]
        stats = [prior_worker.header[
                i] for i in prior_worker.stat_indices]
        for j in jobs:
            cont_params = [j.regression_worker.header[
                    i] for i in j.regression_worker.continuous_parameter_indices]
            disc_params = [j.regression_worker.header[
                    i] for i in j.regression_worker.discrete_parameter_indices]
            self.assertEqual(set(cont_params),
                    self.expected_continuous_parameters)
            self.assertEqual(set(disc_params),
                    self.expected_discrete_parameters)
            rej_params = [j.header[i] for i in j.parameter_indices]
            self.assertEqual(set(cont_params + disc_params) - set(parameters), set())
            self.assertEqual(sorted(parameters), sorted(rej_params))
            result_header = open(j.regression_worker.regress_posterior_path, 
                    'rU').next().strip().split()
            self.assertEqual(sorted(result_header), sorted(cont_params))
            results = self.parse_python_config(
                    j.regression_worker.regress_summary_path)
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

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
    def test_abctoolbox_llr_with_models(self):
        prior_workers = []
        for i in range(2):
            prior_workers.append(self._get_prior(n=100, model_index=i,
                schema='abctoolbox'))
        prior_path = self.get_test_path()
        header_path = self.get_test_path()
        workers.merge_priors(prior_workers,
                prior_path = prior_path,
                header_path = header_path,
                include_header = True)
        prior_worker = prior_workers[0]
        prior_worker.prior_path = prior_path
        obs_worker = self._get_prior(n=2, model_index=1, schema='abctoolbox')
        jobs = workers.assemble_rejection_workers(
                temp_fs = self.temp_fs,
                observed_sims_file = obs_worker.prior_path,
                prior_path = prior_worker.prior_path,
                num_prior_samples = 200,
                num_posterior_samples = 100,
                results_dir = self.results_dir,
                posterior_prefix = self.posterior_prefix,
                regress = True,
                rejection_tool = 'abctoolbox',
                regression_method = 'llr')
        self.assertEqual(len(jobs), 2)
        for j in jobs:
            self.assertFalse(j.finished)
            j.start()
        for j in jobs:
            self.assertTrue(j.finished)
            self.assertTrue(os.path.isfile(j.posterior_path))
            self.assertTrue(self.get_number_of_lines(j.posterior_path), 101)
            self.assertTrue(os.path.isfile(
                    j.regression_worker.regress_summary_path))
            self.assertTrue(os.path.isfile(
                    j.regression_worker.regress_posterior_path))
            self.assertEqual(self.get_number_of_lines(
                    j.regression_worker.regress_posterior_path), 101)
            self.assertTrue(self.get_number_of_lines(
                    j.regression_worker.regress_summary_path) > 10)
        # check regression results
        parameters = [prior_worker.header[
                i] for i in prior_worker.parameter_indices]
        stats = [prior_worker.header[
                i] for i in prior_worker.stat_indices]
        for j in jobs:
            cont_params = [j.regression_worker.header[
                    i] for i in j.regression_worker.continuous_parameter_indices]
            disc_params = [j.regression_worker.header[
                    i] for i in j.regression_worker.discrete_parameter_indices]
            self.assertEqual(set(cont_params),
                    self.expected_continuous_parameters)
            self.assertEqual(set(disc_params),
                    self.expected_discrete_parameters_with_models)
            rej_params = [j.header[i] for i in j.parameter_indices]
            self.assertEqual(set(cont_params + disc_params) - set(parameters), set())
            self.assertEqual(sorted(parameters), sorted(rej_params))
            result_header = open(j.regression_worker.regress_posterior_path, 
                    'rU').next().strip().split()
            self.assertEqual(sorted(result_header), sorted(cont_params))
            results = self.parse_python_config(
                    j.regression_worker.regress_summary_path)
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

class RegressionWorkerTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')
        self.summary_path = self.get_test_path(prefix='post_summary')
        self.regress_posterior_path = self.get_test_path(
                prefix='adjusted_samples')

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
                regress_summary_path = self.summary_path,
                regress_posterior_path = self.regress_posterior_path,
                tolerance = 1.0)
        self.assertFalse(regress_worker.finished)
        regress_worker.start()
        self.assertTrue(regress_worker.finished)
        self.assertEqual(
                self.get_number_of_lines(regress_worker.regress_posterior_path),
                101)
        result_header = open(regress_worker.regress_posterior_path, 
                'rU').next().strip().split()
        cont_params = [regress_worker.header[
                i] for i in regress_worker.continuous_parameter_indices]
        self.assertEqual(sorted(result_header), sorted(cont_params))
        disc_params = [regress_worker.header[
                i] for i in regress_worker.discrete_parameter_indices]
        stats = [regress_worker.header[
                i] for i in regress_worker.stat_indices]
        results = self.parse_python_config(regress_worker.regress_summary_path)
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

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
    def test_regression_with_models(self):
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
                regress_summary_path = self.summary_path,
                regress_posterior_path = self.regress_posterior_path,
                tolerance = 1.0)
        self.assertFalse(regress_worker.finished)
        regress_worker.start()
        self.assertTrue(regress_worker.finished)
        self.assertEqual(
                self.get_number_of_lines(regress_worker.regress_posterior_path),
                201)
        result_header = open(regress_worker.regress_posterior_path, 
                'rU').next().strip().split()
        cont_params = [regress_worker.header[
                i] for i in regress_worker.continuous_parameter_indices]
        self.assertEqual(sorted(result_header), sorted(cont_params))
        disc_params = [regress_worker.header[
                i] for i in regress_worker.discrete_parameter_indices]
        stats = [regress_worker.header[
                i] for i in regress_worker.stat_indices]
        results = self.parse_python_config(regress_worker.regress_summary_path)
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

class PosteriorWorkerTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')
        self.output_dir = self.get_test_subdir(prefix='posterior-worker-output')
        self.output_prefix = os.path.join(self.output_dir, self.test_id)

    def tearDown(self):
        self.tear_down()

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
    def test_without_compression(self):
        MSBAYES_SORT_INDEX.set_index(7)
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

        post_out = self.get_test_path(parent = self.output_dir,
                prefix = self.test_id + '-posterior-sample-')
        density_path = self.get_test_path(parent = self.output_dir,
                prefix = self.test_id + '-posterior-density-')
        post_worker = workers.PosteriorWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                posterior_path = prior_worker.prior_path,
                num_taxon_pairs = 4,
                posterior_out_path = post_out,
                output_prefix = self.output_prefix,
                model_indices = None,
                keep_temps = False,
                regress_posterior_path = density_path,
                omega_threshold = 0.01,
                cv_threshold = 0.01,
                abctoolbox_num_posterior_quantiles = 1000,
                compress = False,
                tag = 'test')

        self.assertFalse(post_worker.finished)
        post_worker.start()
        self.assertTrue(post_worker.finished)
        self.assertTrue(os.path.isfile(post_worker.posterior_out_path))
        self.assertFalse(is_gzipped(post_worker.posterior_out_path))
        self.assertTrue(os.path.isfile(post_worker.regress_posterior_path))
        self.assertFalse(is_gzipped(post_worker.regress_posterior_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.regress_posterior_path),
                1001)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.regress_posterior_path),
                1)
        self.assertTrue(os.path.isfile(post_worker.regress_summary_path))
        self.assertFalse(is_gzipped(post_worker.regress_summary_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.regress_summary_path),
                20)
        self.assertTrue(os.path.isfile(post_worker.posterior_summary_path))
        self.assertTrue(os.path.isfile(post_worker.div_model_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.div_model_results_path),
                6)
        self.assertTrue(os.path.isfile(post_worker.psi_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.psi_results_path),
                5)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.psi_results_path),
                1)
        self.assertTrue(os.path.isfile(post_worker.omega_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.omega_results_path),
                2)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.omega_results_path),
                1)
        self.assertTrue(os.path.isfile(post_worker.cv_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.cv_results_path),
                2)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.cv_results_path),
                1)
        self.assertFalse(os.path.isfile(post_worker.model_results_path))
        MSBAYES_SORT_INDEX.reset_default()

    def test_regression_failure(self):
        MSBAYES_SORT_INDEX.set_index(7)
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

        post_out = self.get_test_path(parent = self.output_dir,
                prefix = self.test_id + '-posterior-sample-')
        density_path = self.get_test_path(parent = self.output_dir,
                prefix = self.test_id + '-posterior-density-')
        post_worker = workers.PosteriorWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                posterior_path = prior_worker.prior_path,
                num_taxon_pairs = 4,
                posterior_out_path = post_out,
                output_prefix = self.output_prefix,
                model_indices = None,
                keep_temps = False,
                regress_posterior_path = density_path,
                omega_threshold = 0.01,
                cv_threshold = 0.01,
                abctoolbox_num_posterior_quantiles = -10,
                abctoolbox_bandwidth = -2.0,
                compress = False,
                tag = 'test')

        self.assertFalse(post_worker.finished)
        post_worker.start()
        self.assertTrue(post_worker.finished)
        self.assertTrue(post_worker.regression_failed)
        self.assertTrue(os.path.isfile(post_worker.posterior_out_path))
        self.assertFalse(is_gzipped(post_worker.posterior_out_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.regress_posterior_path), 0)
        self.assertFalse(os.path.exists(post_worker.regress_summary_path))
        self.assertTrue(os.path.isfile(post_worker.div_model_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.div_model_results_path),
                6)
        self.assertTrue(os.path.isfile(post_worker.psi_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.psi_results_path),
                5)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.psi_results_path),
                1)
        self.assertTrue(os.path.isfile(post_worker.omega_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.omega_results_path),
                2)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.omega_results_path),
                1)
        self.assertTrue(os.path.isfile(post_worker.cv_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.cv_results_path),
                2)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.cv_results_path),
                1)
        self.assertFalse(os.path.isfile(post_worker.model_results_path))
        MSBAYES_SORT_INDEX.reset_default()

    def test_without_compression_with_models(self):
        MSBAYES_SORT_INDEX.set_index(7)
        prior_workers = []
        for i in range(4):
            prior_worker = workers.MsBayesWorker(
                    temp_fs = self.temp_fs,
                    sample_size = 25,
                    config_path = self.cfg_path,
                    schema = 'abctoolbox',
                    model_index = i+1,
                    report_parameters = True)
            prior_worker.start()
            prior_workers.append(prior_worker)
        prior_path = self.get_test_path()
        workers.merge_prior_files(
                paths = [pw.prior_path for pw in prior_workers],
                dest_path = prior_path)
        obs_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                model_index = 1,
                write_stats_file = True,
                report_parameters = True)
        obs_worker.start()

        post_out = self.get_test_path(parent = self.output_dir,
                prefix = self.test_id + '-posterior-sample-')
        density_path = self.get_test_path(parent = self.output_dir,
                prefix = self.test_id + '-posterior-density-')
        post_worker = workers.PosteriorWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                posterior_path = prior_path,
                num_taxon_pairs = 4,
                posterior_out_path = post_out,
                output_prefix = self.output_prefix,
                model_indices = [1,2,3,4],
                keep_temps = False,
                regress_posterior_path = density_path,
                omega_threshold = 0.01,
                cv_threshold = 0.01,
                abctoolbox_num_posterior_quantiles = 1000,
                compress = False,
                tag = 'test')

        self.assertFalse(post_worker.finished)
        post_worker.start()
        self.assertTrue(post_worker.finished)
        self.assertTrue(os.path.isfile(post_worker.posterior_out_path))
        self.assertFalse(is_gzipped(post_worker.posterior_out_path))
        self.assertTrue(os.path.isfile(post_worker.regress_posterior_path))
        self.assertFalse(is_gzipped(post_worker.regress_posterior_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.regress_posterior_path),
                1001)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.regress_posterior_path),
                1)
        self.assertTrue(os.path.isfile(post_worker.regress_summary_path))
        self.assertFalse(is_gzipped(post_worker.regress_summary_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.regress_summary_path),
                20)
        self.assertTrue(os.path.isfile(post_worker.posterior_summary_path))
        self.assertTrue(os.path.isfile(post_worker.div_model_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.div_model_results_path),
                6)
        self.assertTrue(os.path.isfile(post_worker.psi_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.psi_results_path),
                5)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.psi_results_path),
                1)
        self.assertTrue(os.path.isfile(post_worker.omega_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.omega_results_path),
                2)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.omega_results_path),
                1)
        self.assertTrue(os.path.isfile(post_worker.cv_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.cv_results_path),
                2)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.cv_results_path),
                1)
        self.assertTrue(os.path.isfile(post_worker.model_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.model_results_path),
                5)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.model_results_path),
                1)
        MSBAYES_SORT_INDEX.reset_default()

    @unittest.skipIf(TestLevel.get_current_level() < TestLevel.EXHAUSTIVE,
            "EXHAUSTIVE test")
    def test_without_compression_with_one_model(self):
        MSBAYES_SORT_INDEX.set_index(7)
        prior_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 100,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                model_index = 1,
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

        post_out = self.get_test_path(parent = self.output_dir,
                prefix = self.test_id + '-posterior-sample-')
        density_path = self.get_test_path(parent = self.output_dir,
                prefix = self.test_id + '-posterior-density-')
        post_worker = workers.PosteriorWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                posterior_path = prior_worker.prior_path,
                num_taxon_pairs = 4,
                posterior_out_path = post_out,
                output_prefix = self.output_prefix,
                model_indices = None,
                keep_temps = False,
                regress_posterior_path = density_path,
                omega_threshold = 0.01,
                cv_threshold = 0.01,
                abctoolbox_num_posterior_quantiles = 1000,
                compress = False,
                tag = 'test')

        self.assertFalse(post_worker.finished)
        post_worker.start()
        self.assertTrue(post_worker.finished)
        self.assertTrue(os.path.isfile(post_worker.posterior_out_path))
        self.assertFalse(is_gzipped(post_worker.posterior_out_path))
        self.assertTrue(os.path.isfile(post_worker.regress_posterior_path))
        self.assertFalse(is_gzipped(post_worker.regress_posterior_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.regress_posterior_path),
                1001)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.regress_posterior_path),
                1)
        self.assertTrue(os.path.isfile(post_worker.regress_summary_path))
        self.assertFalse(is_gzipped(post_worker.regress_summary_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.regress_summary_path),
                20)
        self.assertTrue(os.path.isfile(post_worker.posterior_summary_path))
        self.assertTrue(os.path.isfile(post_worker.div_model_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.div_model_results_path),
                6)
        self.assertTrue(os.path.isfile(post_worker.psi_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.psi_results_path),
                5)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.psi_results_path),
                1)
        self.assertTrue(os.path.isfile(post_worker.omega_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.omega_results_path),
                2)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.omega_results_path),
                1)
        self.assertTrue(os.path.isfile(post_worker.cv_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.cv_results_path),
                2)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.cv_results_path),
                1)
        self.assertTrue(os.path.isfile(post_worker.model_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.model_results_path),
                2)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.model_results_path),
                1)
        MSBAYES_SORT_INDEX.reset_default()

    def test_with_compression(self):
        MSBAYES_SORT_INDEX.set_index(7)
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

        post_out = self.get_test_path(parent = self.output_dir,
                prefix = self.test_id + '-posterior-sample-')
        density_path = self.get_test_path(parent = self.output_dir,
                prefix = self.test_id + '-posterior-density-')
        post_worker = workers.PosteriorWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                posterior_path = prior_worker.prior_path,
                num_taxon_pairs = 4,
                posterior_out_path = post_out,
                output_prefix = self.output_prefix,
                model_indices = None,
                keep_temps = False,
                regress_posterior_path = density_path,
                omega_threshold = 0.01,
                cv_threshold = 0.01,
                abctoolbox_num_posterior_quantiles = 1000,
                compress = True,
                tag = 'test')

        self.assertFalse(post_worker.finished)
        post_worker.start()
        self.assertTrue(post_worker.finished)
        self.assertTrue(os.path.isfile(post_worker.posterior_out_path))
        self.assertTrue(is_gzipped(post_worker.posterior_out_path))
        self.assertTrue(os.path.isfile(post_worker.regress_posterior_path))
        self.assertTrue(is_gzipped(post_worker.regress_posterior_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.regress_posterior_path),
                1001)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.regress_posterior_path),
                1)
        self.assertTrue(os.path.isfile(post_worker.regress_summary_path))
        self.assertTrue(is_gzipped(post_worker.regress_summary_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.regress_summary_path),
                20)
        self.assertTrue(os.path.isfile(post_worker.posterior_summary_path))
        self.assertTrue(os.path.isfile(post_worker.div_model_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.div_model_results_path),
                6)
        self.assertTrue(os.path.isfile(post_worker.psi_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.psi_results_path),
                5)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.psi_results_path),
                1)
        self.assertTrue(os.path.isfile(post_worker.omega_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.omega_results_path),
                2)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.omega_results_path),
                1)
        self.assertTrue(os.path.isfile(post_worker.cv_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.cv_results_path),
                2)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.cv_results_path),
                1)
        self.assertFalse(os.path.isfile(post_worker.model_results_path))
        MSBAYES_SORT_INDEX.reset_default()

    def test_with_models_no_sort(self):
        MSBAYES_SORT_INDEX.set_index(0)
        prior_workers = []
        for i in range(4):
            prior_worker = workers.MsBayesWorker(
                    temp_fs = self.temp_fs,
                    sample_size = 25,
                    config_path = self.cfg_path,
                    schema = 'abctoolbox',
                    model_index = i+1,
                    report_parameters = True)
            prior_worker.start()
            prior_workers.append(prior_worker)
        prior_path = self.get_test_path()
        workers.merge_prior_files(
                paths = [pw.prior_path for pw in prior_workers],
                dest_path = prior_path)
        obs_worker = workers.MsBayesWorker(
                temp_fs = self.temp_fs,
                sample_size = 1,
                config_path = self.cfg_path,
                schema = 'abctoolbox',
                model_index = 1,
                write_stats_file = True,
                report_parameters = True)
        obs_worker.start()

        post_out = self.get_test_path(parent = self.output_dir,
                prefix = self.test_id + '-posterior-sample-')
        density_path = self.get_test_path(parent = self.output_dir,
                prefix = self.test_id + '-posterior-density-')
        post_worker = workers.PosteriorWorker(
                temp_fs = self.temp_fs,
                observed_path = obs_worker.prior_stats_path,
                posterior_path = prior_path,
                num_taxon_pairs = 4,
                posterior_out_path = post_out,
                output_prefix = self.output_prefix,
                model_indices = [1,2,3,4],
                keep_temps = False,
                regress_posterior_path = density_path,
                omega_threshold = 0.01,
                cv_threshold = 0.01,
                abctoolbox_num_posterior_quantiles = 1000,
                compress = False,
                tag = 'test')

        self.assertFalse(post_worker.finished)
        post_worker.start()
        self.assertTrue(post_worker.finished)
        self.assertTrue(os.path.isfile(post_worker.posterior_out_path))
        self.assertFalse(is_gzipped(post_worker.posterior_out_path))
        self.assertTrue(os.path.isfile(post_worker.regress_posterior_path))
        self.assertFalse(is_gzipped(post_worker.regress_posterior_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.regress_posterior_path),
                1001)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.regress_posterior_path),
                1)
        self.assertTrue(os.path.isfile(post_worker.regress_summary_path))
        self.assertFalse(is_gzipped(post_worker.regress_summary_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.regress_summary_path),
                20)
        self.assertTrue(os.path.isfile(post_worker.posterior_summary_path))
        self.assertTrue(os.path.isfile(post_worker.div_model_results_path))
        self.assertTrue(
                self.get_number_of_lines(post_worker.div_model_results_path) > \
                        2)
        self.assertTrue(os.path.isfile(post_worker.psi_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.psi_results_path),
                5)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.psi_results_path),
                1)
        self.assertTrue(os.path.isfile(post_worker.omega_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.omega_results_path),
                2)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.omega_results_path),
                1)
        self.assertTrue(os.path.isfile(post_worker.cv_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.cv_results_path),
                2)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.cv_results_path),
                1)
        self.assertTrue(os.path.isfile(post_worker.model_results_path))
        self.assertEqual(
                self.get_number_of_lines(post_worker.model_results_path),
                5)
        self.assertEqual(
                self.get_number_of_header_lines(post_worker.model_results_path),
                1)
        current_count = post_worker.div_model_summary[0][1]['count']
        for k, d in post_worker.div_model_summary:
            self.assertEqual(len(k.split(',')), 4)
            self.assertTrue(d['count'] <= current_count)
            current_count = d['count']

        MSBAYES_SORT_INDEX.reset_default()

if __name__ == '__main__':
    unittest.main()

