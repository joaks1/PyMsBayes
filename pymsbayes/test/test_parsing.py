#! /usr/bin/env python

import unittest
import os
import sys
import math
from cStringIO import StringIO

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

    def write_posterior_file(self, path=None):
        if not path:
            path = self.posterior_path
        out = open(path, 'w')
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

    def test_strip_div_model_column(self):
        self.write_posterior_file()
        self.columns.append(('PRI.div.model', [2, 1, 2]))
        div_tmp_path = self.get_test_path(prefix = 'div-model-column')
        self.write_posterior_file(path=div_tmp_path)
        stripped_path = self.get_test_path(prefix = 'stripped')
        strip_div_model_column(div_tmp_path, stripped_path)
        self.assertSameFiles([self.posterior_path, stripped_path])

class ParameterDensityIterTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.pdf_path = package_paths.data_path(
                'abctoolbox_posterior_density_file.txt')

    def tearDown(self):
        self.tear_down()

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

    def tearDown(self):
        self.tear_down()

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

class ParseAbctoolboxSummaryFileTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.psf_path = package_paths.data_path(
                'abctoolbox_posterior_summary_file.txt')

    def tearDown(self):
        self.tear_down()

    def test_parse(self):
        summaries = parse_abctoolbox_summary_file(self.psf_path)
        self.assertEqual(len(summaries), 4)
        self.assertAlmostEqual(summaries['PRI.Psi']['mode'], 1.0)
        self.assertAlmostEqual(summaries['PRI.omega']['mode'], 0.000129132)
        self.assertAlmostEqual(summaries['PRI.Psi']['HPD_99_upper_bound'],
                1.07843)
        self.assertAlmostEqual(summaries['PRI.omega']['HPD_99_upper_bound'],
                0.00152915)
        self.assertAlmostEqual(summaries['PRI.E.t']['quantile_99_lower_bound'],
                0.0317952)

class GetDictFromSpreadsheetTestCase(unittest.TestCase):

    def test_tab_with_head(self):
        sep = '\t'
        header = ['PRI.Psi', 'pi.net.1', 'pi.net.2']
        d = dict(zip(header, [['1','5','3', '2', '6'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s1 = StringIO()
        s2 = StringIO()
        s1.write('{0}\n'.format(sep.join(header)))
        s2.write('{0}\n'.format(sep.join(header)))
        for i in range(3):
            s1.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        for i in range(3, len(d.values()[0])):
            s2.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        _LOG.debug('\nTesting get_dict_from_spreadsheet with spreadsheets:\n'
                '{0}\n{1}'.format(s1.getvalue(), s2.getvalue()))
        s1.seek(0)
        s2.seek(0)
        ret = get_dict_from_spreadsheets([s1, s2], sep = sep)
        self.assertEqual(d, ret)

    def test_tab_without_head(self):
        sep = '\t'
        header = ['PRI.Psi', 'pi.net.1', 'pi.net.2']
        d = dict(zip(header, [['1','5','3', '2', '6'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s1 = StringIO()
        s2 = StringIO()
        for i in range(3):
            s1.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        for i in range(3, len(d.values()[0])):
            s2.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        _LOG.debug('\nTesting get_dict_from_spreadsheet with spreadsheets:\n'
                '{0}\n{1}'.format(s1.getvalue(), s2.getvalue()))
        s1.seek(0)
        s2.seek(0)
        ret = get_dict_from_spreadsheets([s1, s2], sep = sep, header = header)
        self.assertEqual(d, ret)

    def test_comma_with_head(self):
        sep = ','
        header = ['PRI.Psi', 'pi.net.1', 'pi.net.2']
        d = dict(zip(header, [['1','5','3', '2', '6'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s1 = StringIO()
        s2 = StringIO()
        s1.write('{0}\n'.format(sep.join(header)))
        s2.write('{0}\n'.format(sep.join(header)))
        for i in range(3):
            s1.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        for i in range(3, len(d.values()[0])):
            s2.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        _LOG.debug('\nTesting get_dict_from_spreadsheet with spreadsheets:\n'
                '{0}\n{1}'.format(s1.getvalue(), s2.getvalue()))
        s1.seek(0)
        s2.seek(0)
        ret = get_dict_from_spreadsheets([s1, s2], sep = sep)
        self.assertEqual(d, ret)

    def test_comma_without_head(self):
        sep = ','
        header = ['PRI.Psi', 'pi.net.1', 'pi.net.2']
        d = dict(zip(header, [['1','5','3', '2', '6'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s1 = StringIO()
        s2 = StringIO()
        for i in range(3):
            s1.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        for i in range(3, len(d.values()[0])):
            s2.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        _LOG.debug('\nTesting get_dict_from_spreadsheet with spreadsheets:\n'
                '{0}\n{1}'.format(s1.getvalue(), s2.getvalue()))
        s1.seek(0)
        s2.seek(0)
        ret = get_dict_from_spreadsheets([s1, s2], sep = sep, header = header)
        self.assertEqual(d, ret)

    def test_error_mismatch_headers(self):
        sep = '\t'
        header = ['PRI.Psi', 'pi.net.1', 'pi.net.2']
        header2 = ['PRI_Psi', 'pi.net.1', 'pi.net.2']
        d = dict(zip(header, [['1','5','3', '2', '6'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s1 = StringIO()
        s2 = StringIO()
        s1.write('{0}\n'.format(sep.join(header)))
        s2.write('{0}\n'.format(sep.join(header2)))
        for i in range(3):
            s1.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        for i in range(3, len(d.values()[0])):
            s2.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        _LOG.debug('\nTesting get_dict_from_spreadsheet with spreadsheets:\n'
                '{0}\n{1}'.format(s1.getvalue(), s2.getvalue()))
        s1.seek(0)
        s2.seek(0)
        self.assertRaises(Exception, get_dict_from_spreadsheets,
                [s1, s2], sep = sep)

    def test_error_missing_header(self):
        sep = '\t'
        header = ['PRI.Psi', 'pi.net.1', 'pi.net.2']
        d = dict(zip(header, [['1','5','3', '2', '6'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s1 = StringIO()
        s2 = StringIO()
        s1.write('{0}\n'.format(sep.join(header)))
        for i in range(3):
            s1.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        for i in range(3, len(d.values()[0])):
            s2.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        _LOG.debug('\nTesting get_dict_from_spreadsheet with spreadsheets:\n'
                '{0}\n{1}'.format(s1.getvalue(), s2.getvalue()))
        s1.seek(0)
        s2.seek(0)
        self.assertRaises(Exception, get_dict_from_spreadsheets,
                [s1, s2], sep = sep)

class GetStatsByTimeTestCase(unittest.TestCase):

    def test_tab_with_head(self):
        sep = '\t'
        header = ['PRI.t.1', 'PRI.t.2', 'pi.1', 'pi.2', 'pi.net.1', 'pi.net.2']
        d = dict(zip(header, [
                ['1.1','5.3','3.1', '2.7', '6.5'],
                ['1.3','5.1','3.3', '2.8', '6.3'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.21', '0.11', '0.21', '0.31', '0.27'],
                ['0.003', '0.0033', '0.0003', '0.01', '0.0031'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        expected = {'PRI.t': [], 'pi': [], 'pi.net': []}
        for i in range(len(d.values()[0])):
            expected['PRI.t'].extend([float(d['PRI.t.1'][i]),
                    float(d['PRI.t.2'][i])])
            expected['pi'].extend([float(d['pi.1'][i]), float(d['pi.2'][i])])
            expected['pi.net'].extend([float(d['pi.net.1'][i]),
                    float(d['pi.net.2'][i])])
        s1 = StringIO()
        s2 = StringIO()
        s1.write('{0}\n'.format(sep.join(header)))
        s2.write('{0}\n'.format(sep.join(header)))
        for i in range(3):
            s1.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        for i in range(3, len(d.values()[0])):
            s2.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        _LOG.debug('\nTesting get_stats_by_time with spreadsheets:\n'
                '{0}\n{1}'.format(s1.getvalue(), s2.getvalue()))
        s1.seek(0)
        s2.seek(0)
        ret = get_stats_by_time([s1, s2], sep = sep)
        self.assertEqual(expected, ret)

    def test_error_missing_stat_column(self):
        sep = '\t'
        header = ['PRI.t.1', 'PRI.t.2', 'pi.1', 'pi.net.1', 'pi.net.2']
        d = dict(zip(header, [
                ['1.1','5.3','3.1', '2.7', '6.5'],
                ['1.3','5.1','3.3', '2.8', '6.3'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.003', '0.0033', '0.0003', '0.01', '0.0031'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s1 = StringIO()
        s2 = StringIO()
        s1.write('{0}\n'.format(sep.join(header)))
        s2.write('{0}\n'.format(sep.join(header)))
        for i in range(3):
            s1.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        for i in range(3, len(d.values()[0])):
            s2.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        _LOG.debug('\nTesting get_stats_by_time with spreadsheets:\n'
                '{0}\n{1}'.format(s1.getvalue(), s2.getvalue()))
        s1.seek(0)
        s2.seek(0)
        self.assertRaises(Exception, get_stats_by_time, [s1, s2], sep = sep)

    def test_error_extra_stat_column(self):
        sep = '\t'
        header = ['PRI.t.1', 'PRI.t.2', 'pi.1', 'pi.2', 'pi.net.1', 'pi.net.2',
                    'pi.net.3']
        d = dict(zip(header, [
                ['1.1','5.3','3.1', '2.7', '6.5'],
                ['1.3','5.1','3.3', '2.8', '6.3'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.21', '0.11', '0.21', '0.31', '0.27'],
                ['0.003', '0.0033', '0.0003', '0.01', '0.0031'],
                ['0.003', '0.0033', '0.0003', '0.01', '0.0031'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s1 = StringIO()
        s2 = StringIO()
        s1.write('{0}\n'.format(sep.join(header)))
        s2.write('{0}\n'.format(sep.join(header)))
        for i in range(3):
            s1.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        for i in range(3, len(d.values()[0])):
            s2.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        _LOG.debug('\nTesting get_stats_by_time with spreadsheets:\n'
                '{0}\n{1}'.format(s1.getvalue(), s2.getvalue()))
        s1.seek(0)
        s2.seek(0)
        self.assertRaises(Exception, get_stats_by_time, [s1, s2], sep = sep)

    def test_missing_tau_column(self):
        sep = '\t'
        header = ['PRI.t.1', 'pi.1', 'pi.2', 'pi.net.1', 'pi.net.2']
        d = dict(zip(header, [
                ['1.3','5.1','3.3', '2.8', '6.3'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.21', '0.11', '0.21', '0.31', '0.27'],
                ['0.003', '0.0033', '0.0003', '0.01', '0.0031'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s1 = StringIO()
        s2 = StringIO()
        s1.write('{0}\n'.format(sep.join(header)))
        s2.write('{0}\n'.format(sep.join(header)))
        for i in range(3):
            s1.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        for i in range(3, len(d.values()[0])):
            s2.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        _LOG.debug('\nTesting get_stats_by_time with spreadsheets:\n'
                '{0}\n{1}'.format(s1.getvalue(), s2.getvalue()))
        s1.seek(0)
        s2.seek(0)
        self.assertRaises(Exception, get_stats_by_time, [s1, s2], sep = sep)

    def test_error_extra_tau_column(self):
        sep = '\t'
        header = ['PRI.t.1', 'PRI.t.2', 'PRI.t.3', 'pi.1', 'pi.2', 'pi.net.1',
                'pi.net.2']
        d = dict(zip(header, [
                ['1.1','5.3','3.1', '2.7', '6.5'],
                ['1.3','5.1','3.3', '2.8', '6.3'],
                ['1.3','5.1','3.3', '2.8', '6.3'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.21', '0.11', '0.21', '0.31', '0.27'],
                ['0.003', '0.0033', '0.0003', '0.01', '0.0031'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s1 = StringIO()
        s2 = StringIO()
        s1.write('{0}\n'.format(sep.join(header)))
        s2.write('{0}\n'.format(sep.join(header)))
        for i in range(3):
            s1.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        for i in range(3, len(d.values()[0])):
            s2.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        _LOG.debug('\nTesting get_stats_by_time with spreadsheets:\n'
                '{0}\n{1}'.format(s1.getvalue(), s2.getvalue()))
        s1.seek(0)
        s2.seek(0)
        self.assertRaises(Exception, get_stats_by_time, [s1, s2], sep = sep)

    def test_error_stat_numbers(self):
        sep = '\t'
        header = ['PRI.t.1', 'PRI.t.2', 'pi.1', 'pi.3', 'pi.net.1', 'pi.net.2']
        d = dict(zip(header, [
                ['1.1','5.3','3.1', '2.7', '6.5'],
                ['1.3','5.1','3.3', '2.8', '6.3'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.21', '0.11', '0.21', '0.31', '0.27'],
                ['0.003', '0.0033', '0.0003', '0.01', '0.0031'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s1 = StringIO()
        s2 = StringIO()
        s1.write('{0}\n'.format(sep.join(header)))
        s2.write('{0}\n'.format(sep.join(header)))
        for i in range(3):
            s1.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        for i in range(3, len(d.values()[0])):
            s2.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        _LOG.debug('\nTesting get_stats_by_time with spreadsheets:\n'
                '{0}\n{1}'.format(s1.getvalue(), s2.getvalue()))
        s1.seek(0)
        s2.seek(0)
        self.assertRaises(Exception, get_stats_by_time, [s1, s2], sep = sep)

    def test_error_tau_numbers(self):
        sep = '\t'
        header = ['PRI.t.1', 'PRI.t.3', 'pi.1', 'pi.2', 'pi.net.1', 'pi.net.2']
        d = dict(zip(header, [
                ['1.1','5.3','3.1', '2.7', '6.5'],
                ['1.3','5.1','3.3', '2.8', '6.3'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.21', '0.11', '0.21', '0.31', '0.27'],
                ['0.003', '0.0033', '0.0003', '0.01', '0.0031'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s1 = StringIO()
        s2 = StringIO()
        s1.write('{0}\n'.format(sep.join(header)))
        s2.write('{0}\n'.format(sep.join(header)))
        for i in range(3):
            s1.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        for i in range(3, len(d.values()[0])):
            s2.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        _LOG.debug('\nTesting get_stats_by_time with spreadsheets:\n'
                '{0}\n{1}'.format(s1.getvalue(), s2.getvalue()))
        s1.seek(0)
        s2.seek(0)
        self.assertRaises(Exception, get_stats_by_time, [s1, s2], sep = sep)

class DictLineIterTestCase(unittest.TestCase):

    def test_tab_with_head(self):
        sep = '\t'
        header = ['PRI.Psi', 'pi.net.1', 'pi.net.2']
        d = dict(zip(header, [['1','5','3', '2', '6'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s = StringIO()
        s.write('{0}\n'.format(sep.join(header)))
        for i in range(len(d.values()[0])):
            s.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        s.seek(0)
        r = StringIO()
        for row in dict_line_iter(d, sep = sep, header = header):
            r.write(row)
        self.assertEqual(s.getvalue(), r.getvalue())

    def test_tab_without_head(self):
        sep = '\t'
        header = ['PRI.Psi', 'pi.net.1', 'pi.net.2']
        header = sorted(header)
        d = dict(zip(header, [['1','5','3', '2', '6'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s = StringIO()
        s.write('{0}\n'.format(sep.join(header)))
        for i in range(len(d.values()[0])):
            s.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        s.seek(0)
        r = StringIO()
        for row in dict_line_iter(d, sep = sep):
            r.write(row)
        self.assertEqual(s.getvalue(), r.getvalue())

class ParseOmegaResultsFileTestCase(unittest.TestCase):
    def test_parse_omega_file(self):
        omega_results_path = os.path.join(package_paths.TEST_DATA_DIR,
                'pymsbayes-results', 'pymsbayes-output', 'd1', 'm3',
                'd1-m3-s1-1-omega-results.txt')
        results = parse_omega_results_file(omega_results_path)
        self.assertEqual(len(results), 3)
        self.assertAlmostEqual(results['threshold'], 0.01)
        self.assertAlmostEqual(results['prob_less'], 0.02)
        self.assertAlmostEqual(results['prob_less_glm'], 0.0235768970431)

    def test_line_number_exception(self):
        s = StringIO()
        s.write('omega_thresh\tprob_less_than\tglm_prob_less_than\n'
                '0.01\t0.912\t0.989123\n'
                '0.01\t0.912\t0.989123\n')
        s.seek(0)
        self.assertRaises(Exception, parse_omega_results_file, s)
        s = StringIO()
        s.write('omega_thresh\tprob_less_than\tglm_prob_less_than\n')
        s.seek(0)
        self.assertRaises(Exception, parse_omega_results_file, s)
    
    def test_bad_header(self):
        s = StringIO()
        s.write('thresh\tprob_less_than\tglm_prob_less_than\n'
                '0.01\t0.912\t0.989123\n')
        s.seek(0)
        self.assertRaises(Exception, parse_omega_results_file, s)

class ParsePsiResultsFileTestCase(unittest.TestCase):
    def test_parse_psi_file(self):
        psi_results_path = os.path.join(package_paths.TEST_DATA_DIR,
                'pymsbayes-results', 'pymsbayes-output', 'd1', 'm3',
                'd1-m3-s1-1-psi-results.txt')
        results = parse_psi_results_file(psi_results_path)
        self.assertIsInstance(results, dict)
        self.assertEqual(len(results), 3)
        self.assertEqual(sorted(results.keys()), list(range(1, 4)))
        self.assertEqual(sorted(results[1].keys()), sorted(['prob', 'prob_glm']))
        self.assertAlmostEqual(results[1]['prob'], 0.0)
        self.assertAlmostEqual(results[2]['prob'], 0.0)
        self.assertAlmostEqual(results[3]['prob'], 1.0)
        for i in range(1, 4):
            self.assertTrue(math.isnan(results[i]['prob_glm']))

    def test_bad_header(self):
        s = StringIO()
        s.write('num_of_dievents\testimated_prob\tglm_adjusted_prob\n'
                '1\t0.12\t0.23\n'
                '2\t0.88\t0.77\n')
        s.seek(0)
        self.assertRaises(Exception, parse_psi_results_file, s)

    def test_bad_value(self):
        s = StringIO()
        s.write('num_of_div_events\testimated_prob\tglm_adjusted_prob\n'
                '1\t0.12\t0.23\n'
                '2\t0.88\tbogus\n')
        s.seek(0)
        self.assertRaises(Exception, parse_psi_results_file, s)

class ParseModelResultsFile(unittest.TestCase):
    def test_parse_model_file(self):
        model_results_path = os.path.join(package_paths.TEST_DATA_DIR,
                'pymsbayes-results', 'pymsbayes-output', 'd1', 'm123-combined',
                'd1-m123-combined-s1-1-model-results.txt')
        results = parse_model_results_file(model_results_path)
        self.assertIsInstance(results, dict)
        self.assertEqual(len(results), 3)
        self.assertEqual(sorted(results.keys()), list(range(1, 4)))
        self.assertEqual(sorted(results[1].keys()), sorted(['prob', 'prob_glm']))
        self.assertAlmostEqual(results[1]['prob'], 0.36)
        self.assertAlmostEqual(results[2]['prob'], 0.32)
        self.assertAlmostEqual(results[3]['prob'], 0.32)
        self.assertAlmostEqual(results[1]['prob_glm'], 0.316788773658)
        self.assertAlmostEqual(results[2]['prob_glm'], 0.470006429464)
        self.assertAlmostEqual(results[3]['prob_glm'], 0.213204796878)

class DMCSimulationResultsTestCase(unittest.TestCase):
    def setUp(self):
        self.info_path = os.path.join(package_paths.TEST_DATA_DIR,
                'pymsbayes-results', 'pymsbayes-info.txt')
        self.prior_configs = {
            1: package_paths.data_path('negros_panay_3pairs_new_dpp.cfg'),
            2: package_paths.data_path('negros_panay_3pairs_new_uniform.cfg'),
            3: package_paths.data_path('negros_panay_3pairs_new_ushaped.cfg'),
        }
        self.observed_configs = self.prior_configs
        self.observed_paths = {
            1: os.path.join(package_paths.TEST_DATA_DIR, 'pymsbayes-results',
                    'observed-summary-stats', 'observed-1.txt'),
            2: os.path.join(package_paths.TEST_DATA_DIR, 'pymsbayes-results',
                    'observed-summary-stats', 'observed-2.txt'),
            3: os.path.join(package_paths.TEST_DATA_DIR, 'pymsbayes-results',
                    'observed-summary-stats', 'observed-3.txt'),
        }

    def test_init(self):
        results = DMCSimulationResults(self.info_path)
        self.assertIsInstance(results, DMCSimulationResults)
        self.assertEqual(results.info_path, self.info_path)
        self.assertEqual(results.results_dir, os.path.dirname(self.info_path))
        self.assertEqual(results.output_dir, os.path.join(
                os.path.dirname(self.info_path), 'pymsbayes-output'))
        self.assertEqual(sorted(results.observed_index_to_config.keys()),
                range(1, 4))
        self.assertEqual(sorted(results.observed_index_to_path.keys()),
                range(1, 4))
        self.assertEqual(sorted(results.observed_config_to_index.values()),
                range(1, 4))
        self.assertEqual(sorted(results.prior_index_to_config.keys()),
                range(1, 4))
        self.assertEqual(sorted(results.prior_config_to_index.values()),
                range(1, 4))
        self.assertEqual(results.prior_index_to_config, self.prior_configs)
        self.assertEqual(results.observed_index_to_config,
                self.observed_configs)
        self.assertEqual(results.observed_index_to_path, self.observed_paths)

    def test_result_path_iter(self):
        results = DMCSimulationResults(self.info_path)
        expected_1_params = {'PRI.t.1': '2.26089234797',
                'PRI.t.2': '0.90622191875',
                'PRI.t.3': '2.04808328054',
                'PRI.d1Theta.1': '0.000191',
                'PRI.d1Theta.2': '0.001532',
                'PRI.d1Theta.3': '0.000721',
                'PRI.d2Theta.1': '0.003364',
                'PRI.d2Theta.2': '0.000612',
                'PRI.d2Theta.3': '0.000010',
                'PRI.aTheta.1': '0.000277',
                'PRI.aTheta.2': '0.001659',
                'PRI.aTheta.3': '0.000411',
                'PRI.model': '1',
                'PRI.Psi': '3',
                'PRI.var.t': '0.530711173422072',
                'PRI.E.t': '1.73839918242',
                'PRI.omega': '0.305287288897178'}
        expected_1_paths = {
                'summary': os.path.join(package_paths.TEST_DATA_DIR,
                        'pymsbayes-results', 'pymsbayes-output', 'd1', 'm1',
                        'd1-m1-s1-1-posterior-summary.txt'),
                'psi': os.path.join(package_paths.TEST_DATA_DIR,
                        'pymsbayes-results', 'pymsbayes-output', 'd1', 'm1',
                        'd1-m1-s1-1-psi-results.txt'),
                'omega': os.path.join(package_paths.TEST_DATA_DIR,
                        'pymsbayes-results', 'pymsbayes-output', 'd1', 'm1',
                        'd1-m1-s1-1-omega-results.txt'),
                'div-model': os.path.join(package_paths.TEST_DATA_DIR,
                        'pymsbayes-results', 'pymsbayes-output', 'd1', 'm1',
                        'd1-m1-s1-1-div-model-results.txt'),
                'model': os.path.join(package_paths.TEST_DATA_DIR,
                        'pymsbayes-results', 'pymsbayes-output', 'd1', 'm1',
                        'd1-m1-s1-1-model-results.txt'),
                }
        rp_iter = results.result_path_iter(1, 1)
        params, paths = rp_iter.next()
        for params2, paths2 in rp_iter:
            pass
        self.assertEqual(paths, expected_1_paths)
        self.assertEqual(params, expected_1_params)
        expected_2_params = {'PRI.t.1': '0.34052413242',
                'PRI.t.2': '0.78091068192',
                'PRI.t.3': '2.60540618445',
                'PRI.d1Theta.1': '0.000209',
                'PRI.d1Theta.2': '0.000539',
                'PRI.d1Theta.3': '0.002246',
                'PRI.d2Theta.1': '0.001421',
                'PRI.d2Theta.2': '0.002391',
                'PRI.d2Theta.3': '0.000322',
                'PRI.aTheta.1': '0.000631',
                'PRI.aTheta.2': '0.000846',
                'PRI.aTheta.3': '0.000169',
                'PRI.model': '1',
                'PRI.Psi': '3',
                'PRI.var.t': '1.44206914355672',
                'PRI.E.t': '1.24228033293',
                'PRI.omega': '1.16082425627355'}
        expected_2_paths = {
                'summary': os.path.join(package_paths.TEST_DATA_DIR,
                        'pymsbayes-results', 'pymsbayes-output', 'd1', 'm1',
                        'd1-m1-s3-1-posterior-summary.txt'),
                'psi': os.path.join(package_paths.TEST_DATA_DIR,
                        'pymsbayes-results', 'pymsbayes-output', 'd1', 'm1',
                        'd1-m1-s3-1-psi-results.txt'),
                'omega': os.path.join(package_paths.TEST_DATA_DIR,
                        'pymsbayes-results', 'pymsbayes-output', 'd1', 'm1',
                        'd1-m1-s3-1-omega-results.txt'),
                'div-model': os.path.join(package_paths.TEST_DATA_DIR,
                        'pymsbayes-results', 'pymsbayes-output', 'd1', 'm1',
                        'd1-m1-s3-1-div-model-results.txt'),
                'model': os.path.join(package_paths.TEST_DATA_DIR,
                        'pymsbayes-results', 'pymsbayes-output', 'd1', 'm1',
                        'd1-m1-s3-1-model-results.txt'),
                }
        self.assertEqual(paths2, expected_2_paths)
        self.assertEqual(params2, expected_2_params)

        expected_params = {'PRI.t.1': '3.62588223655',
                'PRI.t.2': '0.35823047011',
                'PRI.t.3': '1.57901021456',
                'PRI.d1Theta.1': '0.000623',
                'PRI.d1Theta.2': '0.000730',
                'PRI.d1Theta.3': '0.000031',
                'PRI.d2Theta.1': '0.001571',
                'PRI.d2Theta.2': '0.000739',
                'PRI.d2Theta.3': '0.000103',
                'PRI.aTheta.1': '0.000193',
                'PRI.aTheta.2': '0.000413',
                'PRI.aTheta.3': '0.000278',
                'PRI.model': '2',
                'PRI.Psi': '3',
                'PRI.var.t': '2.72625605426388',
                'PRI.E.t': '1.85437430707333',
                'PRI.omega': '1.47017570501535'}
        expected_paths = {
                'summary': os.path.join(package_paths.TEST_DATA_DIR,
                        'pymsbayes-results', 'pymsbayes-output', 'd2', 'm123-combined',
                        'd2-m123-combined-s3-1-posterior-summary.txt'),
                'psi': os.path.join(package_paths.TEST_DATA_DIR,
                        'pymsbayes-results', 'pymsbayes-output', 'd2', 'm123-combined',
                        'd2-m123-combined-s3-1-psi-results.txt'),
                'omega': os.path.join(package_paths.TEST_DATA_DIR,
                        'pymsbayes-results', 'pymsbayes-output', 'd2', 'm123-combined',
                        'd2-m123-combined-s3-1-omega-results.txt'),
                'div-model': os.path.join(package_paths.TEST_DATA_DIR,
                        'pymsbayes-results', 'pymsbayes-output', 'd2', 'm123-combined',
                        'd2-m123-combined-s3-1-div-model-results.txt'),
                'model': os.path.join(package_paths.TEST_DATA_DIR,
                        'pymsbayes-results', 'pymsbayes-output', 'd2', 'm123-combined',
                        'd2-m123-combined-s3-1-model-results.txt'),
                }
        rp_iter = results.result_path_iter(2, '123-combined')
        for params, paths in rp_iter:
            pass
        self.assertEqual(paths, expected_paths)
        self.assertEqual(params, expected_params)

    def test_result_iter(self):
        results = DMCSimulationResults(self.info_path)
        exp = {'mean_tau': {'true': 1.85437430707333,
                            'mode': ((0.827072720697 + 1.1027636276) / 2),
                            'median': 1.22432526401,
                            'mode_glm': 1.26272},
               'omega': {'true': 1.47017570501535,
                         'mode': (0.262383550541 / 2),
                         'median': 0.564675991641,
                         'mode_glm': 0.239378},
               'psi': {'true': 3,
                       'mode': 3,
                       'mode_glm': float('nan'),
                       'probs': {1: {'prob': 0.0, 'prob_glm': float('nan')},
                                 2: {'prob': 0.0, 'prob_glm': float('nan')},
                                 3: {'prob': 1.0, 'prob_glm': float('nan')}}},
               'model': {'true': 2,
                         'mode': 3,
                         'mode_glm': 2.93939,
                         'probs': {1: {'prob': 0.29, 'prob_glm': 0.240650683366},
                                   2: {'prob': 0.31, 'prob_glm': 0.485304802115},
                                   3: {'prob': 0.4, 'prob_glm': 0.274044514518}}}}
        r_iter = results.result_iter(2, '123-combined')
        for results in r_iter:
            pass
        for k in exp['mean_tau'].iterkeys():
            self.assertAlmostEqual(results['mean_tau'][k], exp['mean_tau'][k])
        for k in exp['omega'].iterkeys():
            self.assertAlmostEqual(results['mean_tau'][k], exp['mean_tau'][k])
        for k in ['true', 'mode']:
            self.assertEqual(results['psi'][k], exp['psi'][k])
            self.assertEqual(results['model'][k], exp['model'][k])
        self.assertTrue(math.isnan(exp['psi']['mode_glm']))
        self.assertTrue(math.isnan(results['psi']['mode_glm']))
        self.assertAlmostEqual(results['model']['mode_glm'],
                exp['model']['mode_glm'])
        l = 'probs'
        for k in ['psi', 'model']:
            for i in range(1, 4):
                for p in ['prob', 'prob_glm']:
                    if math.isnan(exp[k][l][i][p]):
                        self.assertTrue(math.isnan(results[k][l][i][p]))
                    else:
                        self.assertAlmostEqual(results[k][l][i][p],
                                exp[k][l][i][p])
    def test_flat_result_iter(self):
        results = DMCSimulationResults(self.info_path)
        exp = {'mean_tau_true': 1.85437430707333,
               'mean_tau_mode': ((0.827072720697 + 1.1027636276) / 2),
               'mean_tau_median': 1.22432526401,
               'mean_tau_mode_glm': 1.26272,
               'omega_true': 1.47017570501535,
               'omega_mode': (0.262383550541 / 2),
               'omega_median': 0.564675991641,
               'omega_mode_glm': 0.239378,
               'psi_true': 3,
               'psi_mode': 3,
               'psi_mode_glm': float('nan'),
               'psi_1_prob': 0.0,
               'psi_1_prob_glm': float('nan'),
               'psi_2_prob': 0.0,
               'psi_2_prob_glm': float('nan'),
               'psi_3_prob': 1.0,
               'psi_3_prob_glm': float('nan'),
               'model_true': 2,
               'model_mode': 3,
               'model_mode_glm': 2.93939,
               'model_1_prob': 0.29,
               'model_1_prob_glm': 0.240650683366,
               'model_2_prob': 0.31,
               'model_2_prob_glm': 0.485304802115,
               'model_3_prob': 0.4,
               'model_3_prob_glm': 0.274044514518}
        r_iter = results.flat_result_iter(2, '123-combined')
        for results in r_iter:
            pass
        for k in exp.iterkeys():
            if math.isnan(exp[k]):
                self.assertTrue(math.isnan(results[k]))
            elif isinstance(exp[k], float):
                self.assertAlmostEqual(exp[k], results[k])
            else:
                self.assertEqual(exp[k], results[k])

    def test_write_result_summaries(self):
        results = DMCSimulationResults(self.info_path)
        results.write_result_summaries()
        prior_keys = results.prior_index_to_config.keys()
        if results.combined_prior_index:
            prior_keys.append(results.combined_prior_index)
        for i in results.observed_index_to_path.iterkeys():
            for j in prior_keys:
                p = os.path.join(results.get_result_dir(i, j),
                        'results.txt.gz')
                self.assertTrue(os.path.exists(p))
        p = os.path.join(results.get_result_dir(2, '123-combined'),
                'results.txt.gz')
        for r in spreadsheet_iter([p]):
            pass
        exp = {'mean_tau_true': 1.85437430707333,
               'mean_tau_mode': ((0.827072720697 + 1.1027636276) / 2),
               'mean_tau_median': 1.22432526401,
               'mean_tau_mode_glm': 1.26272,
               'omega_true': 1.47017570501535,
               'omega_mode': (0.262383550541 / 2),
               'omega_median': 0.564675991641,
               'omega_mode_glm': 0.239378,
               'psi_true': 3,
               'psi_mode': 3,
               'psi_mode_glm': float('nan'),
               'psi_1_prob': 0.0,
               'psi_1_prob_glm': float('nan'),
               'psi_2_prob': 0.0,
               'psi_2_prob_glm': float('nan'),
               'psi_3_prob': 1.0,
               'psi_3_prob_glm': float('nan'),
               'model_true': 2,
               'model_mode': 3,
               'model_mode_glm': 2.93939,
               'model_1_prob': 0.29,
               'model_1_prob_glm': 0.240650683366,
               'model_2_prob': 0.31,
               'model_2_prob_glm': 0.485304802115,
               'model_3_prob': 0.4,
               'model_3_prob_glm': 0.274044514518}
        for k in exp.iterkeys():
            if math.isnan(exp[k]):
                self.assertTrue(math.isnan(float(r[k])))
            elif isinstance(exp[k], float):
                self.assertAlmostEqual(exp[k], float(r[k]))
            else:
                self.assertEqual(exp[k], int(r[k]))
        for i in results.observed_index_to_path.iterkeys():
            for j in prior_keys:
                p = os.path.join(results.get_result_dir(i, j),
                        'results.txt.gz')
                os.remove(p)

if __name__ == '__main__':
    unittest.main()

