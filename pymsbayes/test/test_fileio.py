#! /usr/bin/env python

import unittest
import os

from pymsbayes.fileio import open, fopen, process_file_arg, FileStream
from pymsbayes.test.support import package_paths
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase

class OpenTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.test_path = self.get_test_path()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')

    def tearDown(self):
        self.tear_down()

    def test_read(self):
        self.assertEqual(FileStream.open_files, set())
        fs = open(self.cfg_path, 'rU')
        self.assertEqual(FileStream.open_files, set([fs]))
        self.assertEqual(fs.read(), fopen(self.cfg_path, 'rU').read())
        fs.close()
        self.assertEqual(FileStream.open_files, set())

class ProcessFileArgTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.test_path = self.get_test_path()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')

    def tearDown(self):
        self.tear_down()
    
    def test_string_object(self):
        f, close = process_file_arg(self.test_path, 'w')
        self.assertIsInstance(f, file)
        self.assertTrue(close)
        self.assertFalse(f.closed)
        f.close()
        f, close = process_file_arg(self.cfg_path, 'rU')
        self.assertIsInstance(f, file)
        self.assertTrue(close)
        self.assertFalse(f.closed)
        f.close()

    def test_file_object(self):
        f = open(self.cfg_path, 'rU')
        f2, close = process_file_arg(f)
        self.assertIsInstance(f2, file)
        self.assertFalse(close)
        self.assertFalse(f2.closed)
        self.assertFalse(f.closed)
        self.assertEqual(f, f2)
        f.close()
        self.assertTrue(f2.closed)
        self.assertTrue(f.closed)

if __name__ == '__main__':
    unittest.main()

