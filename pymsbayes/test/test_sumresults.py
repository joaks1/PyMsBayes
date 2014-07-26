#! /usr/bin/env python

import unittest
import os
import sys
import math
from cStringIO import StringIO

import pymsbayes
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.test.support import package_paths
from pymsbayes.utils import parsing
from pymsbayes.utils import sumresults
from pymsbayes.utils.sumresults import *
from pymsbayes.utils.errors import *
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

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
        self.assertEqual(results.num_taxon_pairs, 3)
        self.assertEqual(results.num_sim_reps, 3)

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
                         'mode_glm': 0.239378,
                         'threshold': 0.01,
                         'prob_less': 0.01,
                         'prob_less_glm': 0.028137104656},
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
            self.assertAlmostEqual(results['omega'][k], exp['omega'][k])
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
               'omega_threshold': 0.01,
               'omega_prob_less': 0.01,
               'omega_prob_less_glm': 0.028137104656,
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
        for r in parsing.spreadsheet_iter([p]):
            pass
        exp = {'mean_tau_true': 1.85437430707333,
               'mean_tau_mode': ((0.827072720697 + 1.1027636276) / 2),
               'mean_tau_median': 1.22432526401,
               'mean_tau_mode_glm': 1.26272,
               'omega_true': 1.47017570501535,
               'omega_mode': (0.262383550541 / 2),
               'omega_median': 0.564675991641,
               'omega_mode_glm': 0.239378,
               'omega_threshold': 0.01,
               'omega_prob_less': 0.01,
               'omega_prob_less_glm': 0.028137104656,
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

class ParseDataKeyFile(unittest.TestCase):
    def test_parse_data_key_file(self):
        data_key_path = os.path.join(package_paths.TEST_DATA_DIR,
                'pymsbayes-results', 'pymsbayes-output', 'data-key.txt')
        data_key_dir = os.path.dirname(data_key_path)
        expected = {
                1: os.path.abspath(os.path.join(data_key_dir,
                    '../observed-summary-stats/observed-1.txt')),
                2: os.path.abspath(os.path.join(data_key_dir,
                    '../observed-summary-stats/observed-2.txt')),
                3: os.path.abspath(os.path.join(data_key_dir,
                    '../observed-summary-stats/observed-3.txt'))}
        observed_paths = parse_data_key_file(data_key_path)
        self.assertEqual(observed_paths, expected)

class UnorderedDivergenceModelResultsTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.unordered_path = package_paths.data_path(
                'div-model-results-unordered.txt')
        self.ordered_path = package_paths.data_path(
                os.path.join('posterior-sample.txt.gz'))
        self.error_ordered_path = package_paths.data_path(
                'div-model-results-ordered.txt')

    def tearDown(self):
        self.tear_down()

    def test_unordered_default(self):
        dmr = sumresults.UnorderedDivergenceModelResults(self.unordered_path)
        self.assertEqual(dmr.inclusion_threshold, None)
        self.assertEqual(dmr.n, 10)
        self.assertEqual(len(dmr.models), 10)
        self.assertAlmostEqual(dmr.cumulative_prob, 0.7688)
        m = dmr.models[0]
        self.assertAlmostEqual(m.prob, 0.1242)
        self.assertAlmostEqual(m.glm_prob, 0.0353091998847)
        self.assertEqual(m.int_partition, [3,2,2,1,1])
        expected = [(3, {'median': 1.51219542433,
                         'hpdi_95': (0.00699511563, 5.92942280407)}),
                    (2, {'median': 6.42942517612,
                         'hpdi_95': (0.77383756206, 11.8377290206)}),
                    (2, {'median': 1.33098309672,
                         'hpdi_95': (0.01390793892, 4.21923713971)}),
                    (1, {'median': 6.19116738868,
                         'hpdi_95': (0.46908780651,13.0106869174)}),
                    (1, {'median': 1.70430466557,
                         'hpdi_95': (0.00071452689, 8.43855632832)})]
        for i, (k, d) in enumerate(m.iter_divergences()):
            k_exp, d_exp = expected[i]
            self.assertEqual(k, k_exp)
            self.assertAlmostEqual(d['median'], d_exp['median'])
            self.assertAlmostEqual(d['hpdi_95'][0], d_exp['hpdi_95'][0])
            self.assertAlmostEqual(d['hpdi_95'][1], d_exp['hpdi_95'][1])
        m = dmr.models[1]
        expected = [(2, {'median': 3.38446408604,
                         'hpdi_95': (0.52752450757,10.6079901713)}),
                    (2, {'median': 1.0075068599,
                         'hpdi_95': (0.00675167243,2.92409287271)}),
                    (1, {'median': 9.59616311192,
                         'hpdi_95': (3.39250060621,14.8779562898)}),
                    (1, {'median': 6.31459966875,
                         'hpdi_95': (1.85189189765,11.1637065608)}),
                    (1, {'median': 3.37663824782,
                         'hpdi_95': (0.82221883054,7.26293836489)}),
                    (1, {'median': 1.77483163725,
                         'hpdi_95': (0.23729279096,4.27691977245)}),
                    (1, {'median': 0.735661610655,
                         'hpdi_95': (0.00535012879,2.29339636862)})]
        for i, (k, d) in enumerate(m.iter_divergences()):
            k_exp, d_exp = expected[i]
            self.assertEqual(k, k_exp)
            self.assertAlmostEqual(d['median'], d_exp['median'])
            self.assertAlmostEqual(d['hpdi_95'][0], d_exp['hpdi_95'][0])
            self.assertAlmostEqual(d['hpdi_95'][1], d_exp['hpdi_95'][1])

    def test_unordered_by_num(self):
        dmr = sumresults.UnorderedDivergenceModelResults(self.unordered_path,
                inclusion_threshold = 3)
        self.assertEqual(dmr.inclusion_threshold, 3)
        self.assertEqual(dmr.n, 3)
        self.assertEqual(len(dmr.models), 3)
        self.assertAlmostEqual(dmr.cumulative_prob, 0.3442)
        m = dmr.models[0]
        self.assertAlmostEqual(m.prob, 0.1242)
        self.assertAlmostEqual(m.glm_prob, 0.0353091998847)
        self.assertEqual(m.int_partition, [3,2,2,1,1])
        expected = [(3, {'median': 1.51219542433,
                         'hpdi_95': (0.00699511563, 5.92942280407)}),
                    (2, {'median': 6.42942517612,
                         'hpdi_95': (0.77383756206, 11.8377290206)}),
                    (2, {'median': 1.33098309672,
                         'hpdi_95': (0.01390793892, 4.21923713971)}),
                    (1, {'median': 6.19116738868,
                         'hpdi_95': (0.46908780651,13.0106869174)}),
                    (1, {'median': 1.70430466557,
                         'hpdi_95': (0.00071452689, 8.43855632832)})]
        for i, (k, d) in enumerate(m.iter_divergences()):
            k_exp, d_exp = expected[i]
            self.assertEqual(k, k_exp)
            self.assertAlmostEqual(d['median'], d_exp['median'])
            self.assertAlmostEqual(d['hpdi_95'][0], d_exp['hpdi_95'][0])
            self.assertAlmostEqual(d['hpdi_95'][1], d_exp['hpdi_95'][1])
        m = dmr.models[1]
        expected = [(2, {'median': 3.38446408604,
                         'hpdi_95': (0.52752450757,10.6079901713)}),
                    (2, {'median': 1.0075068599,
                         'hpdi_95': (0.00675167243,2.92409287271)}),
                    (1, {'median': 9.59616311192,
                         'hpdi_95': (3.39250060621,14.8779562898)}),
                    (1, {'median': 6.31459966875,
                         'hpdi_95': (1.85189189765,11.1637065608)}),
                    (1, {'median': 3.37663824782,
                         'hpdi_95': (0.82221883054,7.26293836489)}),
                    (1, {'median': 1.77483163725,
                         'hpdi_95': (0.23729279096,4.27691977245)}),
                    (1, {'median': 0.735661610655,
                         'hpdi_95': (0.00535012879,2.29339636862)})]
        for i, (k, d) in enumerate(m.iter_divergences()):
            k_exp, d_exp = expected[i]
            self.assertEqual(k, k_exp)
            self.assertAlmostEqual(d['median'], d_exp['median'])
            self.assertAlmostEqual(d['hpdi_95'][0], d_exp['hpdi_95'][0])
            self.assertAlmostEqual(d['hpdi_95'][1], d_exp['hpdi_95'][1])

    def test_unordered_by_prob(self):
        dmr = sumresults.UnorderedDivergenceModelResults(self.unordered_path,
                inclusion_threshold = 0.35)
        self.assertEqual(dmr.inclusion_threshold, 0.35)
        self.assertEqual(dmr.n, 4)
        self.assertEqual(len(dmr.models), 4)
        self.assertAlmostEqual(dmr.cumulative_prob, 0.426)
        m = dmr.models[0]
        self.assertAlmostEqual(m.prob, 0.1242)
        self.assertAlmostEqual(m.glm_prob, 0.0353091998847)
        self.assertEqual(m.int_partition, [3,2,2,1,1])
        expected = [(3, {'median': 1.51219542433,
                         'hpdi_95': (0.00699511563, 5.92942280407)}),
                    (2, {'median': 6.42942517612,
                         'hpdi_95': (0.77383756206, 11.8377290206)}),
                    (2, {'median': 1.33098309672,
                         'hpdi_95': (0.01390793892, 4.21923713971)}),
                    (1, {'median': 6.19116738868,
                         'hpdi_95': (0.46908780651,13.0106869174)}),
                    (1, {'median': 1.70430466557,
                         'hpdi_95': (0.00071452689, 8.43855632832)})]
        for i, (k, d) in enumerate(m.iter_divergences()):
            k_exp, d_exp = expected[i]
            self.assertEqual(k, k_exp)
            self.assertAlmostEqual(d['median'], d_exp['median'])
            self.assertAlmostEqual(d['hpdi_95'][0], d_exp['hpdi_95'][0])
            self.assertAlmostEqual(d['hpdi_95'][1], d_exp['hpdi_95'][1])
        m = dmr.models[1]
        expected = [(2, {'median': 3.38446408604,
                         'hpdi_95': (0.52752450757,10.6079901713)}),
                    (2, {'median': 1.0075068599,
                         'hpdi_95': (0.00675167243,2.92409287271)}),
                    (1, {'median': 9.59616311192,
                         'hpdi_95': (3.39250060621,14.8779562898)}),
                    (1, {'median': 6.31459966875,
                         'hpdi_95': (1.85189189765,11.1637065608)}),
                    (1, {'median': 3.37663824782,
                         'hpdi_95': (0.82221883054,7.26293836489)}),
                    (1, {'median': 1.77483163725,
                         'hpdi_95': (0.23729279096,4.27691977245)}),
                    (1, {'median': 0.735661610655,
                         'hpdi_95': (0.00535012879,2.29339636862)})]
        for i, (k, d) in enumerate(m.iter_divergences()):
            k_exp, d_exp = expected[i]
            self.assertEqual(k, k_exp)
            self.assertAlmostEqual(d['median'], d_exp['median'])
            self.assertAlmostEqual(d['hpdi_95'][0], d_exp['hpdi_95'][0])
            self.assertAlmostEqual(d['hpdi_95'][1], d_exp['hpdi_95'][1])

    def test_ordered(self):
        self.assertRaises(Exception, sumresults.UnorderedDivergenceModelResults,
                self.error_ordered_path)
        dmr = sumresults.UnorderedDivergenceModelResults(self.ordered_path,
                inclusion_threshold = 10)
        self.assertEqual(dmr.inclusion_threshold, 10)
        self.assertEqual(dmr.n, 10)
        self.assertEqual(len(dmr.models), 10)
        self.assertTrue((dmr.cumulative_prob > 0.0) and (
                dmr.cumulative_prob < 1.0))
        prev_prob = 1.0
        for m in dmr.models:
            self.assertTrue(m.prob < prev_prob)
            prev_prob = m.prob
            self.assertEqual(m.glm_prob, None)
            prev_med = 999999.9
            prev_i = 99999
            for i, (k, d) in enumerate(m.iter_divergences()):
                if i == prev_i:
                    self.assertTrue(d['median'] < prev_med)
                prev_i = i
                prev_med = d['median']

class OrderedDivergenceModelCollectionTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.ordered_path = package_paths.data_path(
                os.path.join('posterior-sample.txt.gz'))
        post = parsing.parse_parameters(self.ordered_path)
        self.partition_collection = pymsbayes.utils.stats.PartitionCollection(
                post['taus'])

    def tearDown(self):
        self.tear_down()

    def test_ordered(self):
        dmr = sumresults.OrderedDivergenceModelCollection(self.ordered_path)
        self.assertEqual(dmr.inclusion_threshold, None)
        self.assertAlmostEqual(dmr.cumulative_prob, 1.0)
        self.assertAlmostEqual(
                dmr.prob_of_shared_divergence(range(9)),
                self.partition_collection.prob_clustered(range(9)))
        self.assertAlmostEqual(
                dmr.prob_of_shared_divergence([0, 1]),
                self.partition_collection.prob_clustered([0, 1]))
        self.assertAlmostEqual(
                dmr.prob_of_shared_divergence([0, 3, 7]),
                self.partition_collection.prob_clustered([0, 3, 7]))

class OrderedDivergenceModelResultsTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.unordered_path = package_paths.data_path(
                'div-model-results-unordered.txt')
        self.ordered_path = package_paths.data_path(
                'div-model-results-ordered.txt')
        self.post_ordered_path = package_paths.data_path(
                os.path.join('posterior-sample.txt.gz'))

    def tearDown(self):
        self.tear_down()

    def test_ordered_default(self):
        dmr = sumresults.OrderedDivergenceModelResults(self.ordered_path)
        self.assertEqual(dmr.inclusion_threshold, None)
        self.assertEqual(dmr.n, 10)
        self.assertEqual(len(dmr.models), 10)
        self.assertAlmostEqual(dmr.cumulative_prob, (0.0296 + 
                0.015 + 0.0118 + 0.0095 + 0.006 + 0.0046 + 0.0043 + 0.0036 +
                0.0034 + 0.0032))
        m = dmr.models[0]
        self.assertAlmostEqual(m.prob, 0.0296)
        self.assertAlmostEqual(m.glm_prob, 0.00111786398348)
        self.assertEqual(m.partition, [0,0,0,0,0,0,0,0,0])
        expected = [(list(range(9)), {'median': 2.89806067731,
                         'hpdi_95': (1.62602022373,4.24572749138)})
                         ]
        for i, (indices, d) in enumerate(m.iter_divergences()):
            indices_exp, d_exp = expected[i]
            self.assertEqual(indices, indices_exp)
            self.assertAlmostEqual(d['median'], d_exp['median'])
            self.assertAlmostEqual(d['hpdi_95'][0], d_exp['hpdi_95'][0])
            self.assertAlmostEqual(d['hpdi_95'][1], d_exp['hpdi_95'][1])
        m = dmr.models[1]
        expected = [([0], {'median': 2.5934051444,
                          'hpdi_95': (0.00697993267,6.04007877241)}),
                    ([1], {'median': 2.12107243998,
                           'hpdi_95': (0.04921140543,5.37594378422)}),
                    ([2], {'median': 8.57401519742,
                           'hpdi_95': (2.52257268432,13.5980140814)}),
                    ([3], {'median': 1.16517664273,
                           'hpdi_95': (0.03722408498,3.55286109134)}),
                    ([4], {'median': 4.67285220771,
                           'hpdi_95': (0.23966177613,9.11143589731)}),
                    ([5], {'median': 2.83958905335,
                           'hpdi_95': (0.53354021223,6.76893884696)}),
                    ([6], {'median': 2.9941635737,
                           'hpdi_95': (0.17159218528,7.07017845505)}),
                    ([7], {'median': 2.09198620615,
                           'hpdi_95': (0.04244849557,5.35559715967)}),
                    ([8], {'median': 7.87970675021,
                           'hpdi_95': (3.04965714237,14.1714879714)})]
        for i, (indices, d) in enumerate(m.iter_divergences()):
            indices_exp, d_exp = expected[i]
            self.assertEqual(indices, indices_exp)
            self.assertAlmostEqual(d['median'], d_exp['median'])
            self.assertAlmostEqual(d['hpdi_95'][0], d_exp['hpdi_95'][0])
            self.assertAlmostEqual(d['hpdi_95'][1], d_exp['hpdi_95'][1])
        m = dmr.models[2]
        expected = [([0,1,3,5,6,7], {'median': 1.74088590064,
                          'hpdi_95': (0.42087207209,3.06627558879)}),
                    ([2,4,8], {'median': 7.14681059644,
                           'hpdi_95': (3.96503221188,10.4840797416)})]
        for i, (indices, d) in enumerate(m.iter_divergences()):
            indices_exp, d_exp = expected[i]
            self.assertEqual(indices, indices_exp)
            self.assertAlmostEqual(d['median'], d_exp['median'])
            self.assertAlmostEqual(d['hpdi_95'][0], d_exp['hpdi_95'][0])
            self.assertAlmostEqual(d['hpdi_95'][1], d_exp['hpdi_95'][1])

    def test_ordered_by_num(self):
        dmr = sumresults.OrderedDivergenceModelResults(self.ordered_path,
                inclusion_threshold = 3)
        self.assertEqual(dmr.inclusion_threshold, 3)
        self.assertEqual(dmr.n, 3)
        self.assertEqual(len(dmr.models), 3)
        self.assertAlmostEqual(dmr.cumulative_prob, (0.0296 + 
                0.015 + 0.0118))
        m = dmr.models[0]
        self.assertAlmostEqual(m.prob, 0.0296)
        self.assertAlmostEqual(m.glm_prob, 0.00111786398348)
        self.assertEqual(m.partition, [0,0,0,0,0,0,0,0,0])
        expected = [(list(range(9)), {'median': 2.89806067731,
                         'hpdi_95': (1.62602022373,4.24572749138)})
                         ]
        for i, (indices, d) in enumerate(m.iter_divergences()):
            indices_exp, d_exp = expected[i]
            self.assertEqual(indices, indices_exp)
            self.assertAlmostEqual(d['median'], d_exp['median'])
            self.assertAlmostEqual(d['hpdi_95'][0], d_exp['hpdi_95'][0])
            self.assertAlmostEqual(d['hpdi_95'][1], d_exp['hpdi_95'][1])
        m = dmr.models[1]
        expected = [([0], {'median': 2.5934051444,
                          'hpdi_95': (0.00697993267,6.04007877241)}),
                    ([1], {'median': 2.12107243998,
                           'hpdi_95': (0.04921140543,5.37594378422)}),
                    ([2], {'median': 8.57401519742,
                           'hpdi_95': (2.52257268432,13.5980140814)}),
                    ([3], {'median': 1.16517664273,
                           'hpdi_95': (0.03722408498,3.55286109134)}),
                    ([4], {'median': 4.67285220771,
                           'hpdi_95': (0.23966177613,9.11143589731)}),
                    ([5], {'median': 2.83958905335,
                           'hpdi_95': (0.53354021223,6.76893884696)}),
                    ([6], {'median': 2.9941635737,
                           'hpdi_95': (0.17159218528,7.07017845505)}),
                    ([7], {'median': 2.09198620615,
                           'hpdi_95': (0.04244849557,5.35559715967)}),
                    ([8], {'median': 7.87970675021,
                           'hpdi_95': (3.04965714237,14.1714879714)})]
        for i, (indices, d) in enumerate(m.iter_divergences()):
            indices_exp, d_exp = expected[i]
            self.assertEqual(indices, indices_exp)
            self.assertAlmostEqual(d['median'], d_exp['median'])
            self.assertAlmostEqual(d['hpdi_95'][0], d_exp['hpdi_95'][0])
            self.assertAlmostEqual(d['hpdi_95'][1], d_exp['hpdi_95'][1])
        m = dmr.models[2]
        expected = [([0,1,3,5,6,7], {'median': 1.74088590064,
                          'hpdi_95': (0.42087207209,3.06627558879)}),
                    ([2,4,8], {'median': 7.14681059644,
                           'hpdi_95': (3.96503221188,10.4840797416)})]
        for i, (indices, d) in enumerate(m.iter_divergences()):
            indices_exp, d_exp = expected[i]
            self.assertEqual(indices, indices_exp)
            self.assertAlmostEqual(d['median'], d_exp['median'])
            self.assertAlmostEqual(d['hpdi_95'][0], d_exp['hpdi_95'][0])
            self.assertAlmostEqual(d['hpdi_95'][1], d_exp['hpdi_95'][1])

    def test_ordered_by_prob(self):
        dmr = sumresults.OrderedDivergenceModelResults(self.ordered_path,
                inclusion_threshold = 0.0565)
        self.assertEqual(dmr.inclusion_threshold, 0.0565)
        self.assertEqual(dmr.n, 4)
        self.assertEqual(len(dmr.models), 4)
        self.assertAlmostEqual(dmr.cumulative_prob, (0.0296 + 
                0.015 + 0.0118 + 0.0095))
        m = dmr.models[0]
        self.assertAlmostEqual(m.prob, 0.0296)
        self.assertAlmostEqual(m.glm_prob, 0.00111786398348)
        self.assertEqual(m.partition, [0,0,0,0,0,0,0,0,0])
        expected = [(list(range(9)), {'median': 2.89806067731,
                         'hpdi_95': (1.62602022373,4.24572749138)})
                         ]
        for i, (indices, d) in enumerate(m.iter_divergences()):
            indices_exp, d_exp = expected[i]
            self.assertEqual(indices, indices_exp)
            self.assertAlmostEqual(d['median'], d_exp['median'])
            self.assertAlmostEqual(d['hpdi_95'][0], d_exp['hpdi_95'][0])
            self.assertAlmostEqual(d['hpdi_95'][1], d_exp['hpdi_95'][1])
        m = dmr.models[1]
        expected = [([0], {'median': 2.5934051444,
                          'hpdi_95': (0.00697993267,6.04007877241)}),
                    ([1], {'median': 2.12107243998,
                           'hpdi_95': (0.04921140543,5.37594378422)}),
                    ([2], {'median': 8.57401519742,
                           'hpdi_95': (2.52257268432,13.5980140814)}),
                    ([3], {'median': 1.16517664273,
                           'hpdi_95': (0.03722408498,3.55286109134)}),
                    ([4], {'median': 4.67285220771,
                           'hpdi_95': (0.23966177613,9.11143589731)}),
                    ([5], {'median': 2.83958905335,
                           'hpdi_95': (0.53354021223,6.76893884696)}),
                    ([6], {'median': 2.9941635737,
                           'hpdi_95': (0.17159218528,7.07017845505)}),
                    ([7], {'median': 2.09198620615,
                           'hpdi_95': (0.04244849557,5.35559715967)}),
                    ([8], {'median': 7.87970675021,
                           'hpdi_95': (3.04965714237,14.1714879714)})]
        for i, (indices, d) in enumerate(m.iter_divergences()):
            indices_exp, d_exp = expected[i]
            self.assertEqual(indices, indices_exp)
            self.assertAlmostEqual(d['median'], d_exp['median'])
            self.assertAlmostEqual(d['hpdi_95'][0], d_exp['hpdi_95'][0])
            self.assertAlmostEqual(d['hpdi_95'][1], d_exp['hpdi_95'][1])
        m = dmr.models[2]
        expected = [([0,1,3,5,6,7], {'median': 1.74088590064,
                          'hpdi_95': (0.42087207209,3.06627558879)}),
                    ([2,4,8], {'median': 7.14681059644,
                           'hpdi_95': (3.96503221188,10.4840797416)})]
        for i, (indices, d) in enumerate(m.iter_divergences()):
            indices_exp, d_exp = expected[i]
            self.assertEqual(indices, indices_exp)
            self.assertAlmostEqual(d['median'], d_exp['median'])
            self.assertAlmostEqual(d['hpdi_95'][0], d_exp['hpdi_95'][0])
            self.assertAlmostEqual(d['hpdi_95'][1], d_exp['hpdi_95'][1])

    def test_unordered(self):
        self.assertRaises(Exception, sumresults.OrderedDivergenceModelResults,
                self.unordered_path)
        dmr = sumresults.OrderedDivergenceModelResults(self.post_ordered_path,
                inclusion_threshold = 10)
        self.assertEqual(dmr.inclusion_threshold, 10)
        self.assertEqual(dmr.n, 10)
        self.assertEqual(len(dmr.models), 10)
        self.assertTrue((dmr.cumulative_prob > 0.0) and (
                dmr.cumulative_prob < 1.0))
        prev_prob = 1.0
        for m in dmr.models:
            self.assertTrue(m.prob <= prev_prob)
            prev_prob = m.prob
            self.assertEqual(m.glm_prob, None)
            self.assertEqual(m.partition[0], 0)

if __name__ == '__main__':
    unittest.main()

