#! /usr/bin/env python

import unittest
import os
import sys

from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.test.support import package_paths
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

class ParseParametersTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.posterior_path = self.get_test_path(prefix='summary-')
        self.columns = [('PRI.model', [1, 1, 2]),
                        ('PRI.Psi', [2, 1, 2]),
                        ('PRI.E.t', [0.11, 0.13, 0.12]),
                        ('PRI.omega', [0.01, 0.02, 0.03]),
                        ('PRI.aTheta.1', [0.004, 0.005, 0.007]),
                        ('PRI.aTheta.2', [0.003, 0.002, 0.006]),
                        ('PRI.d1Theta.1', [0.001, 0.009, 0.008]),
                        ('PRI.d1Theta.2', [0.007, 0.003, 0.009]),
                        ('PRI.d2Theta.1', [0.005, 0.005, 0.0001]),
                        ('PRI.d2Theta.2', [0.006, 0.0002, 0.001]),
                        ('PRI.t.1', [0.123, 0.134, 0.343]),
                        ('PRI.t.2', [0.232, 0.134, 0.198]),
                        ]
        self.update_expected()
        self.write_posterior_file()

    def tearDown(self):
        self.tear_down()

    def write_posterior_file(self):
        out = open(self.posterior_path, 'w')
        out.write('{0}\n'.format('\t'.join([k for k, v in self.columns])))
        for i in range(3):
            line = [str(v[i]) for k, v in self.columns if len(v) > i]
            out.write('{0}\n'.format('\t'.join(line)))
        out.close()

    def update_expected(self):
        cols = dict(zip([k for k,v in self.columns],
                [v for k,v in self.columns]))
        self.expected = {}
        if cols.has_key('PRI.E.t'):
            self.expected['mean_tau'] = cols['PRI.E.t']
        if cols.has_key('PRI.omega'):
            self.expected['omega'] = cols['PRI.omega']
        if cols.has_key('PRI.Psi'):
            self.expected['psi'] = cols['PRI.Psi']
        if cols.has_key('PRI.model'):
            self.expected['model'] = cols['PRI.model']
        if cols.has_key('PRI.div.model'):
            self.expected['div_model'] = cols['PRI.div.model']
        tau_keys = [k for k in cols.keys() if k.startswith('PRI.t.')]
        if tau_keys:
            self.expected['taus'] = []
            for i in range(3):
                self.expected['taus'].append([cols[k][i] for k in tau_keys])

    def test_error_extra_psi_column(self):
        self.columns.append(('PRI.Psi', [2, 1, 2]))
        self.write_posterior_file()
        self.assertRaises(ParameterParsingError, parse_parameters,
                self.posterior_path)

    def test_error_missing_cell1(self):
        self.columns[0] = ('PRI.model', [1, '', 2])
        self.write_posterior_file()
        self.assertRaises(ParameterParsingError, parse_parameters,
                self.posterior_path)

    def test_error_missing_cell2(self):
        self.columns[0] = ('PRI.model', [1, 2])
        self.write_posterior_file()
        self.assertRaises(ParameterParsingError, parse_parameters,
                self.posterior_path)

    def test_parse_parameters_simple(self):
        samples = parse_parameters(self.posterior_path)
        self.assertEqual(samples, self.expected)

    def test_parse_parameters_no_model_col(self):
        self.columns.pop(0)
        self.update_expected()
        self.write_posterior_file()
        samples = parse_parameters(self.posterior_path)
        self.assertEqual(samples, self.expected)

    def test_add_div_model_column(self):
        div_model_path = self.get_test_path(prefix = 'posterior-div-models-')
        div_models_to_indices = {'2': 1, '1,1': 2}
        expected_div_model_col = ('PRI.div.model', [2, 1, 2])
        add_div_model_column(self.posterior_path, div_model_path,
                div_models_to_indices)
        self.columns.insert(0, expected_div_model_col)
        self.update_expected()
        samples = parse_parameters(div_model_path)
        self.assertEqual(samples, self.expected)

class ParameterDensityIterTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.pdf_path = package_paths.data_path(
                'abctoolbox_posterior_density_file.txt')

    def test_density_iter(self):
        sums = None
        for i, pd in enumerate(parameter_density_iter(self.pdf_path)):
            if not sums:
                sums = dict(zip([k for k in pd.iterkeys()],
                        [[0.0, 0.0] for k in pd.iterkeys()]))
            for k, val_dens_tup in pd.iteritems():
                sums[k][0] += val_dens_tup[0]
                sums[k][1] += val_dens_tup[1]
        self.assertEqual(i + 1, 10000)
        self.assertAlmostEqual(sums['PRI.Psi'][0], 15000.0000299999)
        self.assertAlmostEqual(sums['PRI.Psi'][1], 10012.60148833639)
        self.assertAlmostEqual(sums['PRI.E.t'][0], 1224.374986974999)
        self.assertAlmostEqual(sums['PRI.E.t'][1], 40869.46144650221)
        self.assertAlmostEqual(sums['PRI.omega'][0], 107.598966652539)
        self.assertAlmostEqual(sums['PRI.omega'][1], 465232.29790955799)

class ParseParameterDensityFileTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.pdf_path = package_paths.data_path(
                'abctoolbox_posterior_density_file.txt')

    def test_parse_parameter_density_file(self):
        pd = parse_parameter_density_file(self.pdf_path)
        sums = dict(zip([k for k in pd.iterkeys()],
                        [[0.0, 0.0] for k in pd.iterkeys()]))
        for k, val_dens_tups in pd.iteritems():
            self.assertEqual(len(val_dens_tups), 10000)
            sums[k][0] += sum([x for x, y in val_dens_tups])
            sums[k][1] += sum([y for x, y in val_dens_tups])
        self.assertAlmostEqual(sums['PRI.Psi'][0], 15000.0000299999)
        self.assertAlmostEqual(sums['PRI.Psi'][1], 10012.60148833639)
        self.assertAlmostEqual(sums['PRI.E.t'][0], 1224.374986974999)
        self.assertAlmostEqual(sums['PRI.E.t'][1], 40869.46144650221)
        self.assertAlmostEqual(sums['PRI.omega'][0], 107.598966652539)
        self.assertAlmostEqual(sums['PRI.omega'][1], 465232.29790955799)

if __name__ == '__main__':
    unittest.main()

