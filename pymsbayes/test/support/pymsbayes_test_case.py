#! /usr/bin/env python

import os
import unittest

from pymsbayes.utils.functions import random_str
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


