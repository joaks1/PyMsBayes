#! /usr/bin/env python

import unittest
import os

from pymsbayes.fileio import (open, fopen, process_file_arg, FileStream,
        is_gzipped, GzipFileStream)
from pymsbayes.test.support import package_paths
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

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

    def test_write(self):
        self.assertEqual(FileStream.open_files, set())
        fs = open(self.test_path, 'w')
        self.assertEqual(FileStream.open_files, set([fs]))
        s = 'This\nis\n\na\n\ttest\n'
        fs.write(s)
        fs.close()
        self.assertEqual(FileStream.open_files, set())
        fs = open(self.test_path, 'rU')
        self.assertEqual(FileStream.open_files, set([fs]))
        self.assertEqual(fs.read(), s)
        fs.close()
        self.assertEqual(FileStream.open_files, set())

class IsGzippedTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')
        self.gz_path = package_paths.data_path(
                'abctoolbox_posterior_density_file.txt.gz')

    def tearDown(self):
        self.tear_down()

    def test_simple(self):
        self.assertFalse(is_gzipped(self.cfg_path))
        self.assertTrue(is_gzipped(self.gz_path))

class ProcessFileArgTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.test_path = self.get_test_path()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')
        self.gz_path = package_paths.data_path(
                'abctoolbox_posterior_density_file.txt.gz')
        self.ungz_path = package_paths.data_path(
                'abctoolbox_posterior_density_file.txt')

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

    def test_read_compressed_file(self):
        gzfs, close = process_file_arg(self.gz_path, 'rb')
        out, close_out = process_file_arg(self.test_path, 'w')
        for line in gzfs:
            out.write(line)
        if close_out:
            out.close()
        if close:
            gzfs.close()
        self.assertSameFiles([self.ungz_path, self.test_path], 
                exclude_line_endings=True)

    def test_write_compressed_file(self):
        fs, close = process_file_arg(self.ungz_path, 'rb')
        out, close_out = process_file_arg(self.test_path, 'wb', compresslevel=9)
        for line in fs:
            out.write(line)
        if close_out:
            out.close()
        if close:
            fs.close()
        self.assertTrue(is_gzipped(self.test_path))
        self.assertSameFiles([self.gz_path, self.test_path],
                exclude_line_endings=True)


class GzipFileStreamTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.test_path = self.get_test_path()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')
        self.gz_path = package_paths.data_path(
                'abctoolbox_posterior_density_file.txt.gz')
        self.ungz_path = package_paths.data_path(
                'abctoolbox_posterior_density_file.txt')

    def tearDown(self):
        self.tear_down()

    def test_read(self):
        gzfs = GzipFileStream(self.gz_path, 'rb')
        gzfs_str = gzfs.read()
        with open(self.ungz_path, 'rb') as fs:
            self.assertEqual(gzfs_str, fs.read())
        gzfs.close()

    def test_read_compressed_read_uncompressed(self):
        gzfs = GzipFileStream(self.gz_path, 'rb')
        out = open(self.test_path, 'wb')
        for line in gzfs:
            out.write(line)
        out.close()
        gzfs.close()
        self.assertSameFiles([self.ungz_path, self.test_path],
                exclude_line_endings=True)

    def test_read_compressed_write_compressed(self):
        gzfs = GzipFileStream(self.gz_path, 'rb')
        out = GzipFileStream(self.test_path, 'w')
        for line in gzfs:
            out.write(line)
        out.close()
        gzfs.close()
        self.assertTrue(is_gzipped(self.test_path))
        self.assertSameFiles([self.gz_path, self.test_path], 
                exclude_line_endings=True)

if __name__ == '__main__':
    unittest.main()

