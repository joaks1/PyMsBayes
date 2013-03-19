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
        self.tempfs = TempFileSystem(
                parent = package_paths.test_path(),
                prefix = 'PyMsBayesTestTemp-')
        self.test_id = 'pymsbayes-' + random_str()

    def tear_down(self):
        self.register_file_system()
        self.tempfs.purge()

    def get_test_path(self, parent=None, prefix='temp', create=False):
        return self.tempfs.get_file_path(parent=parent, prefix=prefix,
                create=create)

    def get_test_subdir(self, parent=None, prefix='temp'):
        return self.tempfs.create_subdir(parent=parent, prefix=prefix)

    def register_file(self, path):
        self.tempfs._register_file(path)

    def register_dir(self, path):
        self.tempfs._register_dir(path)

    def register_file_system(self):
        _LOG.debug('registering test file system...')
        for path, dirs, files, in os.walk(self.tempfs.base_dir):
            for f in files:
                if f.startswith(self.test_id):
                    self.register_file(os.path.join(path, f))
            for d in dirs:
                if d.startswith(self.test_id):
                    self.register_dir(os.path.join(path, d))

