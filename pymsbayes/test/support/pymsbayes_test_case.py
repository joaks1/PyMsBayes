#! /usr/bin/env python

import os
import unittest
from cStringIO import StringIO
import re
import random

from pymsbayes.workers import (TAU_PATTERNS,
        MODEL_PATTERNS,
        D_THETA_PATTERNS,
        A_THETA_PATTERNS,
        PSI_PATTERNS,
        MEAN_TAU_PATTERNS,
        OMEGA_PATTERNS,
        HEADER_PATTERN)
from pymsbayes.config import MsBayesConfig
from pymsbayes.utils.stats import SampleSummarizer
from pymsbayes.utils import probability
from pymsbayes.utils.functions import (random_str, process_file_arg,
        get_indices_of_patterns)
from pymsbayes.utils.tempfs import TempFileSystem
from pymsbayes.test.support import package_paths
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class PyMsBayesTestCase(unittest.TestCase):
    
    def set_up(self):
        self.temp_fs = TempFileSystem(
                parent = package_paths.test_path(),
                prefix = 'PyMsBayesTestTemp-')
        self.test_id = 'pymsbayes-' + random_str()

    def tear_down(self):
        self.register_file_system()
        self.temp_fs.purge()

    def get_test_path(self, parent=None, prefix='temp', create=False):
        return self.temp_fs.get_file_path(parent=parent, prefix=prefix,
                create=create)

    def get_test_subdir(self, parent=None, prefix='temp'):
        return self.temp_fs.create_subdir(parent=parent, prefix=prefix)

    def register_file(self, path):
        self.temp_fs._register_file(path)

    def register_dir(self, path):
        self.temp_fs._register_dir(path)

    def register_file_system(self):
        _LOG.debug('registering test file system...')
        for path, dirs, files, in os.walk(self.temp_fs.base_dir):
            for f in files:
                if f.startswith(self.test_id):
                    self.register_file(os.path.join(path, f))
            for d in dirs:
                if d.startswith(self.test_id):
                    self.register_dir(os.path.join(path, d))

    def get_expected_indices(self, num_pairs, dummy_column=True,
            parameters_reported=True):
        num_summary_params = 4
        num_params = 4*num_pairs
        num_default_stats = 4*num_pairs
        start = 0
        if dummy_column:
            start = 1
        param_indices = range(start, start+num_summary_params)
        start += num_summary_params
        if parameters_reported:
            param_indices += range(start, start+num_params)
            start += num_params
        stat_indices = range(start, start+num_default_stats)
        return param_indices, stat_indices
    
    def prior_file_is_valid(self, prior_path, num_of_rows, num_of_columns=None):
        try:
            prior_file = open(prior_path, 'rU')
        except:
            _LOG.error('prior invalid: could not open prior path {0}'.format(
                    prior_path))
            return False
        nrows = 0
        for i, line in enumerate(prior_file):
            nrows += 1
            if not num_of_columns:
                num_of_columns = len(line.strip().split())
            ncols = len(line.strip().split())
            if num_of_columns != ncols:
                _LOG.error('prior invalid: num of columns at line {0} is {1} '
                        'NOT {2}'.format(i+1, ncols, num_of_columns))
                return False
        if num_of_rows != nrows:
            _LOG.error('prior invalid: num of rows is {0} NOT {1}'.format(
                    nrows, num_of_rows))
            return False
        return True

    def get_config_from_msbayes_worker(self, msbayes_worker):
        return MsBayesConfig(msbayes_worker.config_path)

    def get_parameter_summaries_from_msbayes_worker(self, msbayes_worker,
            shuffle_taus=True):
        s = dict(zip(
            [i for i in msbayes_worker.parameter_indices],
            [SampleSummarizer(name=msbayes_worker.header[i]) for i in msbayes_worker.parameter_indices]))
        f = open(msbayes_worker.prior_path, 'rU')
        ncols = None
        for line_idx, row in enumerate(f):
            if not ncols:
                ncols = len(row.strip().split())
            if HEADER_PATTERN.match(row.strip()):
                continue
            r = row.strip().split()
            assert len(r) == ncols
            if shuffle_taus: # because taus are sorted in prior files
                psi_index = get_indices_of_patterns(msbayes_worker.header,
                        PSI_PATTERNS)[0]
                tau_indices = get_indices_of_patterns(msbayes_worker.header,
                        TAU_PATTERNS)
                psi = int(r[psi_index])
                taus = [float(r[i]) for i in tau_indices]
                self.assertEqual(psi, len(set(taus)))
                random.shuffle(taus)
                for n, i in enumerate(tau_indices):
                    s[i].add_sample(taus[n])
                p_set = set(msbayes_worker.parameter_indices) - set(tau_indices)
                p = sorted(list(p_set))
                for i in p:
                    s[i].add_sample(float(r[i]))
            else:
                for i in msbayes_worker.parameter_indices:
                    s[i].add_sample(float(r[i]))
        return s

    def assertPriorIsValid(self, msbayes_worker, places=2):
        self.assertTrue(msbayes_worker.finished)
        param_sums = self.get_parameter_summaries_from_msbayes_worker(
                msbayes_worker)
        for s in param_sums.itervalues():
            self.assertEqual(s.n, msbayes_worker.sample_size)
        cfg = self.get_config_from_msbayes_worker(msbayes_worker)
        psi_indices = get_indices_of_patterns(msbayes_worker.header,
                PSI_PATTERNS)
        self.assertEqual(len(psi_indices), 1)
        model_indices = get_indices_of_patterns(msbayes_worker.header,
                MODEL_PATTERNS)
        if not msbayes_worker.model_index is None:
            self.assertEqual(len(model_indices), 1)
        else:
            self.assertEqual(len(model_indices), 0)
        tau_indices = get_indices_of_patterns(msbayes_worker.header,
                TAU_PATTERNS)
        a_theta_indices = get_indices_of_patterns(msbayes_worker.header,
                A_THETA_PATTERNS)
        d_theta_indices = get_indices_of_patterns(msbayes_worker.header,
                D_THETA_PATTERNS)
        if msbayes_worker.report_parameters:
            self.assertEqual(len(tau_indices), cfg.npairs)
            self.assertEqual(len(a_theta_indices), cfg.npairs)
            self.assertEqual(len(d_theta_indices), 2*cfg.npairs)
        else:
            self.assertEqual(len(tau_indices), 0)
            self.assertEqual(len(a_theta_indices), 0)
            self.assertEqual(len(d_theta_indices), 0)
        for i in psi_indices:
            self.assertSampleIsFromDistribution(param_sums[i], cfg.psi, places)
        for i in tau_indices:
            self.assertSampleIsFromDistribution(param_sums[i], cfg.tau, places)
        for i in a_theta_indices:
            self.assertSampleIsFromDistribution(param_sums[i], cfg.a_theta, places)
        for i in d_theta_indices:
            self.assertSampleIsFromDistribution(param_sums[i], cfg.d_theta, places)
                    
    def assertSampleIsFromDistribution(self, sample_sum, dist, places=2):
        if isinstance(dist, probability.DiscreteUniformDistribution):
            self.assertEqual(sample_sum.minimum, dist.minimum)
            self.assertEqual(sample_sum.maximum, dist.maximum)
        else:
            if dist.minimum != float('-inf') or dist.minimum != float('inf'):
                self.assertAlmostEqual(sample_sum.minimum, dist.minimum, places)
            if dist.maximum != float('-inf') or dist.maximum != float('inf'):
                self.assertAlmostEqual(sample_sum.maximum, dist.maximum, places)
        self.assertAlmostEqual(sample_sum.mean, dist.mean, places)
        self.assertAlmostEqual(sample_sum.variance, dist.variance, places)

    def assertApproxEqual(self, x, y, percent_tol=1e-6):
        eq = (((abs(x-y) / ((abs(x)+abs(y))/2))*100) < percent_tol)
        if not eq:
            _LOG.error('x ({0}) and y ({1}) are not equal'.format(x, y))
        self.assertTrue(eq)

    def files_equal(self, f1, f2):
        equal = True
        diffs = []
        f1, c1 = process_file_arg(f1)
        f2, c2 = process_file_arg(f2)
        line = 0
        f1_end = False
        f2_end = False
        lines_left = True
        while True:
            line += 1
            if f1_end == False:
                try:
                    l1 = f1.next()
                except (StopIteration, EOFError):
                    f1_end = line
                    pass
            if f2_end == False:
                try:
                    l2 = f2.next()
                except (StopIteration, EOFError):
                    f2_end = line
                    pass
            if f1_end != False and f2_end != False:
                break
            if f1_end == False and f2_end == False and l1 != l2:
                diffs.append(line)
                equal = False
        if f1_end != f2_end:
            mn = min([f1_end, f2_end])
            mx = max([f1_end, f2_end])
            diffs.extend(range(mn, mx+1))
            equal = False
        assert len(diffs) == len(set(diffs))
        if c1:
            f1.close()
        if c2:
            f2.close()
        return equal, diffs

    def assertSameFiles(self, files):
        all_equal = True
        diffs = StringIO()
        f1 = files.pop(0)
        for f2 in files:
            equal, diff_list = self.files_equal(f1, f2)
            if not equal:
                all_equal = False
                n1 = f1
                if not isinstance(n1, str):
                    n1 = f1.name
                n2 = f2
                if not isinstance(n2, str):
                    n2 = f2.name
                diffs.write('{0} and {1} differ at lines:\n\t{2}\n'.format(
                        n1, n2, ','.join([str(i) for i in diff_list])))
        if not all_equal:
            _LOG.error('files are not equal:\n{0}\n'.format(diffs.getvalue()))
        self.assertTrue(all_equal)

    def assertSameDistributions(self, d1, d2):
        self.assertEqual(d1.name, d2.name)
        self.assertEqual(str(d1), str(d2))
        self.assertEqual(d1.minimum, d2.minimum)
        self.assertEqual(d1.maximum, d2.maximum)
        self.assertEqual(d1.mean, d2.mean)
        self.assertEqual(d1.variance, d2.variance)

