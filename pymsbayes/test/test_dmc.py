#! /usr/bin/env python

import unittest
import os
import sys
import re
import random

from pymsbayes.test.support import package_paths
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.test import TestLevel, test_enabled
from pymsbayes.utils import GLOBAL_RNG
from pymsbayes.fileio import process_file_arg
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class DmcTestCase(PyMsBayesTestCase):
    RESULT_DIR_PATTERN = re.compile(r'pymsbayes-results[-]*\d*')
    BASE_DIR_PATTERN = re.compile(r'pymsbayes-output[-]*\d*')
    DATA_DIR_PATTERN = re.compile(r'd\d+')
    MODEL_DIR_PATTERN = re.compile(r'm\d+')

    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')
        self.cfg_path2 = package_paths.data_path('4pairs_1locus_maxt5.cfg')
        self.seed = GLOBAL_RNG.randint(1, 999999999)
        self.rng = random.Random()
        self.rng.seed(self.seed)
        self.output_dir = self.get_test_subdir(prefix='dmc-test-')
        self.output_prefix = self.temp_fs.token_id

    def tearDown(self):
        for base, dirs, files in os.walk(self.temp_fs.base_dir):
            for d in dirs:
                if self.BASE_DIR_PATTERN.match(d) or \
                        self.RESULT_DIR_PATTERN.match(d) or \
                        self.DATA_DIR_PATTERN.match(d) or \
                        self.MODEL_DIR_PATTERN.match(d) or \
                        d == 'prior-stats-summaries' or \
                        d == 'observed-summary-stats':
                    self.temp_fs._register_dir(os.path.join(base, d))
        self.tear_down()

    def get_result_paths(self, observed_idx, model_str, sim_idx,
            iter_count, compressed = False, output_dir=None):
        if not output_dir:
            output_dir = self.output_dir
        prefix = '{0}d{1}-{2}-s{3}-{4}'.format(self.output_prefix,
                observed_idx, model_str, sim_idx, iter_count)
        p = os.path.join(output_dir, 'pymsbayes-results',
                'pymsbayes-output', 'd' + str(observed_idx), model_str,
                prefix)
        paths = {}
        paths['sample'] = p + '-posterior-sample.txt'
        paths['summary'] = p + '-posterior-summary.txt'
        paths['psi'] = p + '-psi-results.txt'
        paths['model'] = p + '-model-results.txt'
        paths['omega'] = p + '-omega-results.txt'
        paths['div'] = p + '-div-model-results.txt'
        paths['glm-summary'] = p + '-glm-posterior-summary.txt'
        paths['glm-density'] = p + '-glm-posterior-density-estimates.txt'
        if compressed:
            paths['sample'] += '.gz'
            paths['glm-summary'] += '.gz'
            paths['glm-density'] += '.gz'
        paths['prior-dir'] = os.path.join(output_dir, 'pymsbayes-results',
                'pymsbayes-output', 'prior-stats-summaries')
        return paths

    def _has_non_sorted_results(self, div_model_path):
        length = None
        f, close = process_file_arg(div_model_path)
        f.next() # header
        for line in f:
            l = line.strip()
            if l:
                div_model_key = l.split()[0]
                div_model = div_model_key.split(',')
                if not length:
                    length = len(div_model)
                if length != len(div_model):
                    f.close()
                    return False
        f.close()
        return True

    def _exe_dmc(self, args, stdout = None, stderr = None, return_code = 0,
            output_dir = None):
        if not output_dir:
            output_dir = self.output_dir
        args += ['--output-prefix', self.output_prefix,
                 '--output-dir', output_dir]
        self._exe_script(script_name = 'dmc.py', args = args,
                stdout = stdout, stderr = stderr,
                return_code = return_code)

    def test_help(self):
        self._exe_dmc(['-h'], return_code=0)

    def test_no_args(self):
        self._exe_dmc([], return_code=2)

    def test_bogus_arg(self):
        self._exe_dmc(['--cheeseburger'], return_code=2)

    def test_repeatability(self):
        args = ['-o', self.cfg_path,
                '-p', self.cfg_path,
                '-r', 1,
                '-n', 400,
                '--num-posterior-samples', 200,
                '--num-standardizing-samples', 300,
                '-q', 100,
                '--np', 4,
                '--seed', self.seed,
                '--debug']
        self._exe_dmc(args, return_code=0)
        results1 = self.get_result_paths(1, 'm1', 1, 1)
        self.assertTrue(os.path.exists(results1['prior-dir']))
        self.assertTrue(os.path.exists(results1['sample']))

        out_dir1 = self.get_test_subdir(prefix='repeat-')

        args = ['-o', self.cfg_path,
                '-p', self.cfg_path,
                '-r', 1,
                '-n', 400,
                '--num-posterior-samples', 200,
                '--num-standardizing-samples', 300,
                '-q', 100,
                '--np', 4,
                '--seed', self.seed,
                '--debug']
        self._exe_dmc(args, return_code=0, output_dir = out_dir1)
        results2 = self.get_result_paths(1, 'm1', 1, 1, output_dir = out_dir1)
        self.assertTrue(os.path.exists(results2['prior-dir']))
        self.assertTrue(os.path.exists(results2['sample']))

        self.assertNotEqual(results1['prior-dir'], results2['prior-dir'])
        self.assertNotEqual(results1['sample'], results2['sample'])
        self.assertSameFiles([results1['sample'], results2['sample']])

        out_dir2 = self.get_test_subdir(prefix='repeat2-')

        args = ['-o', self.cfg_path,
                '-p', self.cfg_path,
                '-r', 1,
                '-n', 400,
                '--num-posterior-samples', 200,
                '--num-standardizing-samples', 300,
                '-q', 100,
                '--np', 4,
                '--seed', self.seed + 1,
                '--debug']
        self._exe_dmc(args, return_code=0, output_dir = out_dir2)
        results3 = self.get_result_paths(1, 'm1', 1, 1, output_dir = out_dir2)
        self.assertTrue(os.path.exists(results3['prior-dir']))
        self.assertTrue(os.path.exists(results3['sample']))
        self.assertNotEqual(results1['prior-dir'], results3['prior-dir'])
        self.assertNotEqual(results1['sample'], results3['sample'])
        equal, diffs = self.files_equal(results1['sample'], results3['sample'])
        self.assertFalse(equal)

    def test_two_stage(self):
        args = ['-o', self.cfg_path,
                '-p', self.cfg_path,
                '-r', 1,
                '-n', 400,
                '--num-posterior-samples', 200,
                '--num-standardizing-samples', 300,
                '-q', 100,
                '--np', 4,
                '--seed', self.seed,
                '--debug']
        self._exe_dmc(args, return_code=0)
        results1 = self.get_result_paths(1, 'm1', 1, 1)
        self.assertTrue(os.path.exists(results1['prior-dir']))
        self.assertTrue(os.path.exists(results1['sample']))

        out_dir1 = self.get_test_subdir(prefix='repeat-')

        args = ['-o', self.cfg_path,
                '-p', self.cfg_path,
                '-r', 1,
                '-n', 400,
                '--num-posterior-samples', 200,
                '--num-standardizing-samples', 300,
                '-q', 100,
                '--np', 4,
                '--generate-samples-only',
                '--seed', self.seed,
                '--debug']
        self._exe_dmc(args, return_code=0, output_dir = out_dir1)
        results2 = self.get_result_paths(1, 'm1', 1, 1, output_dir = out_dir1)
        self.assertTrue(os.path.exists(results2['prior-dir']))
        self.assertFalse(os.path.exists(results2['sample']))

        out_dir2 = self.get_test_subdir(prefix='repeat2-')

        args = ['-o', self.cfg_path,
                '-p', results2['prior-dir'],
                '-r', 1,
                '-n', 400,
                '--num-posterior-samples', 200,
                '--num-standardizing-samples', 300,
                '-q', 100,
                '--np', 4,
                '--seed', self.seed,
                '--debug']
        self._exe_dmc(args, return_code=0, output_dir = out_dir2)
        results3 = self.get_result_paths(1, 'm1', 1, 1, output_dir = out_dir2)

        self.assertNotEqual(results1['prior-dir'], results3['prior-dir'])
        self.assertNotEqual(results1['sample'], results3['sample'])
        self.assertSameFiles([results1['sample'], results3['sample']])

    def test_sort_index(self):
        args = ['-o', self.cfg_path,
                '-p', self.cfg_path,
                '-r', 1,
                '-n', 400,
                '--num-posterior-samples', 200,
                '--num-standardizing-samples', 300,
                '-q', 100,
                '--np', 4,
                '--sort-index', 7,
                '--seed', self.seed,
                '--debug']
        self._exe_dmc(args, return_code=0)
        results1 = self.get_result_paths(1, 'm1', 1, 1)
        self.assertTrue(os.path.exists(results1['prior-dir']))
        self.assertTrue(os.path.exists(results1['sample']))
        self.assertFalse(self._has_non_sorted_results(results1['div']))

        out_dir1 = self.get_test_subdir(prefix='repeat-')

        args = ['-o', self.cfg_path,
                '-p', self.cfg_path,
                '-r', 1,
                '-n', 400,
                '--num-posterior-samples', 200,
                '--num-standardizing-samples', 300,
                '-q', 100,
                '--np', 4,
                '--sort-index', 0,
                '--seed', self.seed,
                '--debug']
        self._exe_dmc(args, return_code=0, output_dir = out_dir1)
        results2 = self.get_result_paths(1, 'm1', 1, 1, output_dir = out_dir1)
        self.assertTrue(os.path.exists(results2['prior-dir']))
        self.assertTrue(os.path.exists(results2['sample']))
        self.assertTrue(self._has_non_sorted_results(results2['div']))

    def test_no_global_estimate(self):
        args = ['-o', self.cfg_path,
                '-p', self.cfg_path, self.cfg_path2,
                '-r', 1,
                '-n', 400,
                '--num-posterior-samples', 200,
                '--num-standardizing-samples', 300,
                '-q', 100,
                '--np', 4,
                '--seed', self.seed,
                '--debug']
        self._exe_dmc(args, return_code=0)
        results1_1 = self.get_result_paths(1, 'm1', 1, 1)
        results1_2 = self.get_result_paths(1, 'm2', 1, 1)
        results1_c = self.get_result_paths(1, 'm12-combined', 1, 1)
        self.assertTrue(os.path.exists(results1_1['sample']))
        self.assertTrue(os.path.exists(results1_2['sample']))
        self.assertTrue(os.path.exists(results1_c['sample']))

        out_dir1 = self.get_test_subdir(prefix='repeat-')

        args = ['-o', self.cfg_path,
                '-p', self.cfg_path, self.cfg_path2,
                '-r', 1,
                '-n', 400,
                '--num-posterior-samples', 200,
                '--num-standardizing-samples', 300,
                '-q', 100,
                '--np', 4,
                '--no-global-estimate',
                '--seed', self.seed,
                '--debug']
        self._exe_dmc(args, return_code=0, output_dir = out_dir1)
        results2_1 = self.get_result_paths(1, 'm1', 1, 1, output_dir = out_dir1)
        results2_2 = self.get_result_paths(1, 'm2', 1, 1, output_dir = out_dir1)
        results2_c = self.get_result_paths(1, 'm12-combined', 1, 1, output_dir = out_dir1)
        self.assertTrue(os.path.exists(results2_1['sample']))
        self.assertTrue(os.path.exists(results2_2['sample']))
        self.assertFalse(os.path.exists(results2_c['sample']))

        self.assertSameFiles([results1_1['sample'], results2_1['sample']])
        self.assertSameFiles([results1_2['sample'], results2_2['sample']])

if __name__ == '__main__':
    unittest.main()

