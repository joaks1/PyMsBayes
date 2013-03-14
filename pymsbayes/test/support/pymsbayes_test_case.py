#! /usr/bin/env python

import os
import unittest

from pymsbayes.utils.tempfs import TempFileSystem
from pymsbayes.test.support import package_paths
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class PyMsBayesTestCase(unittest.TestCase):
    
    def set_up(self):
        self.tempfs = TempFileSystem(
                parent = package_paths.test_path()
                prefix = 'PyMsBayesTestTemp-')

    def tear_down(self):
        self.tempfs.purge()

    def get_test_path(self, parent=None, prefix='temp', create=False):
        return self.tempfs.get_file_path(parent=parent, prefix=prefix,
                create=create)

    def get_test_subdir(self, parent=None, prefix='temp'):
        return self.tempfs.create_subdir(parent=parent, prefix=prefix)

