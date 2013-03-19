#! /usr/bin/env python

import unittest
import os

from pymsbayes.utils import functions
from pymsbayes.test.support import package_paths
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase

class MkdrTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        p = ['mkdr', 'test', 'dir']
        self.path = os.path.join(self.tempfs.base_dir, *p)
        for i in range(1, len(p)+1):
            self.register_dir(os.path.join(
                    self.tempfs.base_dir,
                    *p[:i]))

    def tearDown(self):
        self.tear_down()
    
    def test_mkdr(self):
        self.assertFalse(os.path.exists(self.path))
        functions.mkdr(self.path)
        self.assertTrue(os.path.exists(self.path))
        self.assertTrue(os.path.isdir(self.path))
        functions.mkdr(self.path)
        self.assertTrue(os.path.exists(self.path))
        self.assertTrue(os.path.isdir(self.path))

class WhereisTestCase(unittest.TestCase):
    def test_whereis(self):
        path = functions.whereis('more')
        self.assertIsNotNone(path)
        self.assertTrue(os.path.exists(path))
        path = functions.whereis('x_bogus_exe_name_x')
        self.assertIsNone(path)

class IsFileTestCase(unittest.TestCase):
    def setUp(self):
        self.file = package_paths.data_path("4pairs_1locus.cfg")
        self.bogus_file = package_paths.data_path("bogusdatafilename")
    
    def test_is_file(self):
        self.assertFalse(functions.is_file(None))
        self.assertFalse(functions.is_file(self.bogus_file))
        self.assertTrue(functions.is_file(self.file))
        
class IsExecutableTestCase(unittest.TestCase):
    def setUp(self):
        self.file = package_paths.data_path("4pairs_1locus.cfg")
        self.bogus_file = package_paths.data_path("bogusdatafilename")
        self.exe = package_paths.bin_path('msbayes.pl')
    
    def test_is_executable(self):
        self.assertFalse(functions.is_executable(None))
        self.assertFalse(functions.is_executable(self.bogus_file))
        self.assertFalse(functions.is_executable(self.file))
        self.assertTrue(functions.is_executable(self.exe))

if __name__ == '__main__':
    unittest.main()
