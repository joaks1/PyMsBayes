#! /usr/bin/env python

import os
import unittest
from cStringIO import StringIO
import re
import random
from configobj import ConfigObj

from pymsbayes.fileio import process_file_arg, open, FileStream
from pymsbayes.utils.parsing import (TAU_PATTERNS,
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
from pymsbayes.utils.functions import random_str, get_indices_of_patterns
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
        self.assertEqual(FileStream.open_files, set())

    def get_test_path(self, parent=None, prefix='temp'):
        return self.temp_fs.get_file_path(parent=parent, prefix=prefix)

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
    
    def prior_file_is_valid(self, prior_path, num_of_samples,
            num_of_columns=None):
        try:
            prior_file = open(prior_path, 'rU')
        except:
            _LOG.error('prior invalid: could not open prior path {0}'.format(
                    prior_path))
            return False
        nrows = 0
        for i, line in enumerate(prior_file):
            if nrows == 0 and HEADER_PATTERN.match(line):
                pass
            else:
                nrows += 1
            if not num_of_columns:
                num_of_columns = len(line.strip().split())
            ncols = len(line.strip().split())
            if num_of_columns != ncols:
                _LOG.error('prior invalid: num of columns at line {0} is {1} '
                        'NOT {2}'.format(i+1, ncols, num_of_columns))
                return False
        prior_file.close()
        if num_of_samples != nrows:
            _LOG.error('prior invalid: num of rows is {0} NOT {1}'.format(
                    nrows, num_of_samples))
            return False
        return True
    
    def get_number_of_lines(self, path):
        f, close = process_file_arg(path)
        count = 0
        for l in f:
            count += 1
        if close:
            f.close()
        return count

    def get_number_of_header_lines(self, path):
        f, close = process_file_arg(path)
        count = 0
        for l in f:
            if HEADER_PATTERN.match(l.strip()):
                count += 1
        if close:
            f.close()
        return count

    def parse_python_config(self, path):
        return ConfigObj(path)

    def get_config_from_msbayes_workers(self, msbayes_workers):
        cfgs = [MsBayesConfig(w.config_path) for w in msbayes_workers]
        self.assertSameConfigs(cfgs)
        return cfgs[0]

    def assertSameConfigs(self, cfgs):
        configs = list(cfgs)
        c1 = configs.pop(0)
        for c2 in cfgs:
            self.assertEqual(c1.npairs, c2.npairs)
            self.assertSameDistributions(c1.psi, c2.psi)
            self.assertSameDistributions(c1.tau, c2.tau)
            self.assertSameDistributions(c1.theta, c2.theta)
            self.assertSameDistributions(c1.a_theta, c2.a_theta)
            self.assertSameDistributions(c1.d_theta, c2.d_theta)
            self.assertSameDistributions(c1.recombination, c2.recombination)
            self.assertSameDistributions(c1.migration, c2.migration)

    def get_parameter_summaries_from_msbayes_workers(self, msbayes_workers,
            shuffle_taus=True):
        msbayes_workers = list(msbayes_workers)
        s = dict(zip(
            [i for i in msbayes_workers[0].parameter_indices],
            [SampleSummarizer(
                tag=msbayes_workers[0].header[i]) for i in msbayes_workers[
                    0].parameter_indices]))
        ncols = None
        header = msbayes_workers[0].header
        pi = msbayes_workers[0].parameter_indices
        for w in msbayes_workers:
            self.assertEqual(w.header, header)
            self.assertEqual(w.parameter_indices, pi)
            f = open(w.prior_path, 'rU')
            for line_idx, row in enumerate(f):
                if not ncols:
                    ncols = len(row.strip().split())
                if HEADER_PATTERN.match(row.strip()):
                    continue
                r = row.strip().split()
                assert len(r) == ncols
                if shuffle_taus: # because taus are sorted in prior files
                    psi_index = get_indices_of_patterns(w.header,
                            PSI_PATTERNS)[0]
                    tau_indices = get_indices_of_patterns(w.header,
                            TAU_PATTERNS)
                    psi = int(r[psi_index])
                    taus = [float(r[i]) for i in tau_indices]
                    self.assertEqual(psi, len(set(taus)))
                    random.shuffle(taus)
                    for n, i in enumerate(tau_indices):
                        s[i].add_sample(taus[n])
                    p_set = set(w.parameter_indices) - set(tau_indices)
                    p = sorted(list(p_set))
                    for i in p:
                        s[i].add_sample(float(r[i]))
                else:
                    for i in w.parameter_indices:
                        s[i].add_sample(float(r[i]))
            f.close()
        return s

    def assertPriorIsPrecise(self, msbayes_workers, places=2):
        msbayes_workers = list(msbayes_workers)
        self.assertWorkersFinished(msbayes_workers)
        param_sums = self.get_parameter_summaries_from_msbayes_workers(
                msbayes_workers)
        sample_size = 0
        for w in msbayes_workers:
            sample_size += w.sample_size
        for s in param_sums.itervalues():
            self.assertEqual(s.n, sample_size)
        cfg = self.get_config_from_msbayes_workers(msbayes_workers)
        psi_indices = get_indices_of_patterns(msbayes_workers[0].header,
                PSI_PATTERNS)
        self.assertEqual(len(psi_indices), 1)
        model_indices = get_indices_of_patterns(msbayes_workers[0].header,
                MODEL_PATTERNS)
        if not msbayes_workers[0].model_index is None:
            self.assertEqual(len(model_indices), 1)
        else:
            self.assertEqual(len(model_indices), 0)
        tau_indices = get_indices_of_patterns(msbayes_workers[0].header,
                TAU_PATTERNS)
        a_theta_indices = get_indices_of_patterns(msbayes_workers[0].header,
                A_THETA_PATTERNS)
        d_theta_indices = get_indices_of_patterns(msbayes_workers[0].header,
                D_THETA_PATTERNS)
        if msbayes_workers[0].report_parameters:
            self.assertEqual(len(tau_indices), cfg.npairs)
            self.assertEqual(len(a_theta_indices), cfg.npairs)
            self.assertEqual(len(d_theta_indices), 2*cfg.npairs)
        else:
            self.assertEqual(len(tau_indices), 0)
            self.assertEqual(len(a_theta_indices), 0)
            self.assertEqual(len(d_theta_indices), 0)
        _LOG.debug('\n{0}\n'.format('\n'.join(
                [str(param_sums[i]) for i in sorted(param_sums.iterkeys())])))
        for i in psi_indices:
            self.assertSampleIsFromDistribution(param_sums[i], cfg.psi,
                    places=places)
        for i in tau_indices:
            self.assertSampleIsFromDistribution(param_sums[i], cfg.tau,
                    places=places)
        for i in a_theta_indices:
            self.assertSampleIsFromDistribution(param_sums[i], cfg.a_theta,
                    places=places)
        for i in d_theta_indices:
            self.assertSampleIsFromDistribution(param_sums[i], cfg.d_theta,
                    mean_adj=cfg.theta.mean,
                    max_adj=cfg.theta.maximum,
                    compare_variance=False,
                    places=places)

    def assertPriorIsAccurate(self, msbayes_workers, places=2):
        msbayes_workers = list(msbayes_workers)
        self.assertWorkersFinished(msbayes_workers)
        pass

    def assertPriorIsValid(self, msbayes_workers, places=2):
        msbayes_workers = list(msbayes_workers)
        self.assertWorkersFinished(msbayes_workers)
        self.assertPriorIsPrecise(msbayes_workers, places=places)
        self.assertPriorIsAccurate(msbayes_workers, places=places)

    def assertWorkersFinished(self, msbayes_workers):
        for w in msbayes_workers:
            self.assertTrue(w.finished)
                    
    def assertSampleIsFromDistribution(self, sample_sum, dist, places=2,
            mean_adj=1,
            max_adj=1,
            compare_variance=True):
        if isinstance(dist, probability.DiscreteUniformDistribution):
            self.assertEqual(sample_sum.minimum, dist.minimum)
            self.assertEqual(sample_sum.maximum, dist.maximum)
        else:
            if dist.minimum != float('-inf') or dist.minimum != float('inf'):
                self.assertAlmostEqual(sample_sum.minimum, dist.minimum, places)
            if dist.maximum != float('-inf') or dist.maximum != float('inf'):
                self.assertAlmostEqual(sample_sum.maximum, dist.maximum*max_adj, places)
        self.assertAlmostEqual(sample_sum.mean, dist.mean*mean_adj, places)
        if compare_variance:
            self.assertAlmostEqual(sample_sum.variance, dist.variance, places)

    def assertApproxEqual(self, x, y, percent_tol=1e-6):
        eq = (((abs(x-y) / ((abs(x)+abs(y))/2))*100) < percent_tol)
        if not eq:
            _LOG.error('x ({0}) and y ({1}) are not equal'.format(x, y))
        self.assertTrue(eq)

    def files_equal(self, f1, f2, exclude_line_endings=False):
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
            if exclude_line_endings:
                l1 = l1.strip()
                l2 = l2.strip()
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

    def assertSameFiles(self, files, exclude_line_endings=False):
        files = list(files)
        all_equal = True
        diffs = StringIO()
        f1 = files.pop(0)
        for f2 in files:
            equal, diff_list = self.files_equal(f1, f2, exclude_line_endings)
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

    def assertSameUnsortedFiles(self, files):
        files = list(files)
        all_equal = True
        diffs = StringIO()
        f1, close = process_file_arg(files.pop(0))
        lines1 = sorted(f1.readlines())
        for f in files:
            f2, close2 = process_file_arg(f)
            lines2 = sorted(f2.readlines())
            if len(lines1) != len(lines2):
                all_equal = False
                diffs.write('{0} ({1}) and {2} ({3}) have different '
                        'number of lines\n'.format(f1.name, len(lines1),
                                f2.name, len(lines2)))
            for i in range(min(len(lines1), len(lines2))):
                if lines1[i].strip().split() != lines2[i].strip().split():
                    all_equal = False
                    diffs.write('{0} and {1} differ at sorted index '
                            '{2}\n'.format(f1.name, f2.name, i))

            if close2:
                f2.close()
        if not all_equal:
            _LOG.error('files are not equal after sorting:\n{0}\n'.format(
                    diffs.getvalue()))
        self.assertTrue(all_equal)
        if close:
            f1.close()

    def same_samples(self, sample1, sample2, places = 4, num_mismatches = 0):
        if len(sample1) != len(sample2):
            return False
        for i in range(len(sample1)):
            if round(float(sample1[i]) - float(sample2[i]), places) != 0:
                if num_mismatches < 1:
                    return False
                num_mismatches -= 1
        return True
            
    def assertSameSamples(self, files, columns_to_ignore = [], header = True,
            places = 5, num_mismatches_per_sample = 0,
            num_sample_mismatches = 0):
        files = list(files)
        all_equal = True
        diffs = StringIO()
        f1, close = process_file_arg(files.pop(0))
        f1_lines = f1.readlines()
        indices = [i for i in range(len(
                f1_lines[0].strip().split())) if i not in columns_to_ignore]
        h1 = []
        if header:
            head = f1_lines.pop(0).strip().split()
            h1 = [head[i] for i in indices]
        lines1 = sorted(f1_lines)
        for f in files:
            f2, close2 = process_file_arg(f)
            f2_lines = f2.readlines()
            h2 = []
            if header:
                head = f2_lines.pop(0).strip().split()
                h2 = [head[i] for i in indices]
                if h1 != h2:
                    all_equal = False
                    diffs.write('{0} and {1} have different headers; not '
                            'comparing further\n'.format(
                                    f1.name, f2.name))
                    continue
            lines2 = sorted(f2_lines)
            if len(lines1) != len(lines2):
                all_equal = False
                diffs.write('{0} ({1}) and {2} ({3}) have different '
                        'number of lines\n'.format(f1.name, len(lines1),
                                f2.name, len(lines2)))
            n_matches = 0
            n_mismatches = 0
            for l1 in lines1:
                found = False
                for l2 in lines2:
                    values1 = l1.strip().split()
                    values2 = l2.strip().split()
                    v1 = [float(values1[x]) for x in indices]
                    v2 = [float(values2[x]) for x in indices]
                    if self.same_samples(v1, v2, places = places,
                            num_mismatches = num_mismatches_per_sample):
                        found = True
                if found:
                    n_matches += 1
                else:
                    n_mismatches += 1
            if n_mismatches > 0:
                if n_mismatches > num_sample_mismatches:
                    all_equal = False
                diffs.write('{0} and {1}\nhave {2} mismatching samples and '
                        'share {3} samples\n'.format(
                                f1.name, f2.name, n_mismatches, n_matches))
            if close2:
                f2.close()
        if diffs.getvalue() != '':
            _LOG.error('files are not equal after sorting:\n{0}\n'.format(
                    diffs.getvalue()))
        self.assertTrue(all_equal)
        if close:
            f1.close()
    
    def assertSameDistributions(self, d1, d2):
        self.assertEqual(d1.name, d2.name)
        self.assertEqual(str(d1), str(d2))
        self.assertEqual(d1.minimum, d2.minimum)
        self.assertEqual(d1.maximum, d2.maximum)
        self.assertEqual(d1.mean, d2.mean)
        self.assertEqual(d1.variance, d2.variance)

    def assertSameIntegerPartitions(self, integer_partitions):
        ips = list(integer_partitions)
        ip1 = ips.pop(0)
        for ip2 in ips:
            self.assertEqual(ip1._initialized, ip2._initialized)
            self.assertEqual(ip1.n, ip2.n)
            self.assertEqual(ip1.key, ip2.key)
            self.assertEqual(ip1.integer_partition, ip2.integer_partition)
            self.assertEqual(ip1._items, ip2._items)

    def assertSamePartitions(self, partitions):
        ps = list(partitions)
        p1 = ps.pop(0)
        for p2 in ps:
            self.assertEqual(p1._initialized, p2._initialized)
            self.assertEqual(p1.n, p2.n)
            self.assertEqual(p1.key, p2.key)
            self.assertEqual(p1.partition, p2.partition)
            self.assertEqual(p1.values, p2.values)

