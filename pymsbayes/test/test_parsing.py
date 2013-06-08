#! /usr/bin/env python

import unittest
import os
import sys

from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.utils.parsing import *
from pymsbayes.utils.errors import *
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class ParseSummaryFileTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.summary_path = self.get_test_path(prefix='summary-')
        self.lines = ['stat.1\tstat.2\tstat.3\n',
                      '1.3\t4.78\t5.73\n',
                      '5.6\t7.45\t6.78\n',
                      '10000\t10000\t10000\n']
        self.write_sum_file()
        self.expected = {
                'stat.1': {'mean': 1.3, 'std_deviation': 5.6, 'n': 10000},
                'stat.2': {'mean': 4.78, 'std_deviation': 7.45, 'n': 10000},
                'stat.3': {'mean': 5.73, 'std_deviation': 6.78, 'n': 10000}}
        self.expected_header = ['stat.1', 'stat.2', 'stat.3']

    def tearDown(self):
        self.tear_down()

    def write_sum_file(self):
        out = open(self.summary_path, 'w')
        out.writelines(self.lines)
        out.close()

    def test_line_error(self):
        self.lines = self.lines[:3]
        self.write_sum_file()
        self.assertRaises(SummaryFileParsingError, parse_summary_file,
                self.summary_path)
        self.lines.append('5.6\t7.45\t6.78\n')
        self.lines.append('5.6\t7.45\t6.78\n')
        self.assertRaises(SummaryFileParsingError, parse_summary_file,
                self.summary_path)

    def test_column_error(self):
        self.lines[2] = '5.6\t7.45\t\n'
        self.write_sum_file()
        self.assertRaises(SummaryFileParsingError, parse_summary_file,
                self.summary_path)

    def test_basic(self):
        d, h = parse_summary_file(self.summary_path)
        self.assertEqual(self.expected, d)
        self.assertEqual(self.expected_header, h)

    def test_empty_lines(self):
        self.lines.insert(0, '\n')
        self.lines.insert(0, '\n')
        self.lines.append('\n')
        self.lines.append('\n')
        self.write_sum_file()
        d, h = parse_summary_file(self.summary_path)
        self.assertEqual(self.expected, d)
        self.assertEqual(self.expected_header, h)

if __name__ == '__main__':
    unittest.main()

