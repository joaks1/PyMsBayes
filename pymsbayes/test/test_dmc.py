#! /usr/bin/env python

import unittest
import os
import sys
import re
import random

from pymsbayes.test.support import package_paths
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.test import TestLevel, test_enabled
from pymsbayes.utils import GLOBAL_RNG
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class DmcTestCase(PyMsBayesTestCase):
    RESULT_DIR_PATTERN = re.compile(r'pymsbayes-results[-]*\d*')
    BASE_DIR_PATTERN = re.compile(r'pymsbayes-output[-]*\d*')
    DATA_DIR_PATTERN = re.compile(r'd\d+')
    MODEL_DIR_PATTERN = re.compile(r'm\d+')

    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')
        self.seed = GLOBAL_RNG.randint(1, 999999999)
        self.rng = random.Random()
        self.rng.seed(self.seed)
        self.output_dir = self.get_test_subdir(prefix='dmc-test-')
        self.output_prefix = self.temp_fs.token_id

    def tearDown(self):
        for base, dirs, files in os.walk(self.output_dir):
            for d in dirs:
                if self.BASE_DIR_PATTERN.match(d) or \
                        self.RESULT_DIR_PATTERN.match(d) or \
                        self.DATA_DIR_PATTERN.match(d) or \
                        self.MODEL_DIR_PATTERN.match(d) or \
                        d == 'prior-stats-summaries' or \
                        d == 'observed-summary-stats':
                    self.temp_fs._register_dir(os.path.join(base, d))
        self.tear_down()

    def get_result_paths(self, observed_idx, model_str, sim_idx,
            iter_count, compressed = False):
        prefix = '{0}d{1}-{2}-s{3}-{4}'.format(self.output_prefix,
                observed_idx, model_str, sim_idx, iter_count)
        p = os.path.join(self.output_dir, 'pymsbayes-results',
                'pymsbayes-output', 'd' + str(observed_idx), model_str,
                prefix)
        paths = {}
        paths['sample'] = p + '-posterior-sample.txt'
        paths['summary'] = p + '-posterior-summary.txt'
        paths['psi'] = p + '-psi-results.txt'
        paths['model'] = p + '-model-results.txt'
        paths['omega'] = p + '-omega-results.txt'
        paths['div'] = p + '-div-model-results.txt'
        paths['glm-summary'] = p + '-glm-posterior-summary.txt'
        paths['glm-density'] = p + '-glm-posterior-density-estimates.txt'
        if compressed:
            paths['sample'] += '.gz'
            paths['glm-summary'] += '.gz'
            paths['glm-density'] += '.gz'
        paths['prior-dir'] = os.path.join(self.output_dir, 'pymsbayes-results',
                'pymsbayes-output', 'prior-stats-summaries')
        return paths

    def _exe_dmc(self, args, stdout = None, stderr = None, return_code = 0):
        args += ['--output-prefix', self.output_prefix,
                 '--output-dir', self.output_dir]
        self._exe_script(script_name = 'dmc.py', args = args,
                stdout = stdout, stderr = stderr,
                return_code = return_code)

    def test_help(self):
        self._exe_dmc(['-h'], return_code=0)

    def test_no_args(self):
        self._exe_dmc([], return_code=2)

    def test_bogus_arg(self):
        self._exe_dmc(['--cheeseburger'], return_code=2)

    def test_basic(self):
        args = ['-o', self.cfg_path,
                '-p', self.cfg_path,
                '-r', 1,
                '-n', 400,
                '--num-posterior-samples', 200,
                '--num-standardizing-samples', 300,
                '-q', 100,
                '--np', 4,
                '--seed', self.seed,
                '--debug']
        self._exe_dmc(args, return_code=0)
        results = self.get_result_paths(1, 'm1', 1, 1)
        _LOG.debug('{0}\n'.format(results['prior-dir']))
        _LOG.debug('{0}\n'.format(results['sample']))
        self.assertTrue(os.path.exists(results['prior-dir']))
        self.assertTrue(os.path.exists(results['sample']))

    # def test_abc_team_multiple_models_multiple_obs_previous_priors(self):
    #     obs_worker1 = MsBayesWorker(
    #             temp_fs = self.temp_fs,
    #             sample_size = 2,
    #             config_path = self.cfg_path,
    #             schema = 'abctoolbox',
    #             write_stats_file = True)
    #     obs_worker1.start()
    #     obs_worker2 = MsBayesWorker(
    #             temp_fs = self.temp_fs,
    #             sample_size = 2,
    #             config_path = self.cfg_path,
    #             schema = 'abctoolbox',
    #             write_stats_file = True)
    #     obs_worker2.start()

    #     num_taxon_pairs = 4
    #     num_prior_samples = 2000
    #     num_processors = 4
    #     batch_size = 200
    #     num_standardizing_samples = 800
    #     num_posterior_samples = 100
    #     num_posterior_density_quantiles = 100
    #     num_batches = num_prior_samples / batch_size

    #     abct1 = ABCTeam(
    #             temp_fs = self.temp_fs,
    #             observed_stats_files = [obs_worker1.prior_stats_path,
    #                     obs_worker2.prior_stats_path],
    #             num_taxon_pairs = num_taxon_pairs,
    #             config_paths = [self.cfg_path, self.cfg_path2],
    #             num_prior_samples = num_prior_samples,
    #             num_processors = num_processors,
    #             num_standardizing_samples = num_standardizing_samples,
    #             num_posterior_samples = num_posterior_samples,
    #             num_posterior_density_quantiles = num_posterior_density_quantiles,
    #             batch_size = batch_size,
    #             output_dir = self.output_dir,
    #             output_prefix = self.test_id,
    #             rng = self.rng,
    #             abctoolbox_bandwidth = None,
    #             omega_threshold = 0.01,
    #             reporting_frequency = None,
    #             compress = False,
    #             keep_temps = False,
    #             global_estimate_only = False)
    #     self.assertFalse(abct1.finished)
    #     abct1.run()
    #     self.assertTrue(abct1.finished)
    #     self.assertEqual(abct1.num_samples_generated, 4000)
    #     self.assertEqual(abct1.num_samples_summarized, 1600)
    #     self.assertEqual(abct1.num_samples_processed[1], 8000)
    #     self.assertEqual(abct1.num_samples_processed[2], 8000)

    #     self.rng.seed(self.seed)

    #     abct2 = ABCTeam(
    #             temp_fs = self.temp_fs,
    #             observed_stats_files = [obs_worker1.prior_stats_path,
    #                     obs_worker2.prior_stats_path],
    #             num_taxon_pairs = num_taxon_pairs,
    #             config_paths = [self.cfg_path, self.cfg_path2],
    #             num_prior_samples = num_prior_samples,
    #             num_processors = num_processors,
    #             num_standardizing_samples = num_standardizing_samples,
    #             num_posterior_samples = num_posterior_samples,
    #             num_posterior_density_quantiles = num_posterior_density_quantiles,
    #             batch_size = batch_size,
    #             output_dir = self.output_dir,
    #             output_prefix = self.test_id,
    #             rng = self.rng,
    #             abctoolbox_bandwidth = None,
    #             omega_threshold = 0.01,
    #             reporting_frequency = None,
    #             compress = False,
    #             keep_temps = False,
    #             global_estimate_only = False,
    #             generate_prior_samples_only = True)
    #     self.assertFalse(abct2.finished)
    #     abct2.run()
    #     self.assertTrue(abct2.finished)
    #     self.assertEqual(abct2.num_samples_generated, 4000)
    #     self.assertEqual(abct2.num_samples_summarized, 1600)

    #     abct3 = ABCTeam(
    #             temp_fs = self.temp_fs,
    #             observed_stats_files = [obs_worker1.prior_stats_path,
    #                     obs_worker2.prior_stats_path],
    #             num_taxon_pairs = num_taxon_pairs,
    #             previous_prior_dir = abct2.summary_dir,
    #             num_processors = num_processors,
    #             num_posterior_samples = num_posterior_samples,
    #             num_posterior_density_quantiles = num_posterior_density_quantiles,
    #             output_dir = self.output_dir,
    #             output_prefix = self.test_id,
    #             rng = self.rng,
    #             abctoolbox_bandwidth = None,
    #             omega_threshold = 0.01,
    #             reporting_frequency = 1,
    #             compress = False,
    #             keep_temps = False,
    #             global_estimate_only = False,
    #             generate_prior_samples_only = False)
    #     self.assertFalse(abct3.finished)
    #     abct3.run()
    #     self.assertTrue(abct3.finished)
    #     self.assertEqual(abct3.num_samples_processed[1], 8000)
    #     self.assertEqual(abct3.num_samples_processed[2], 8000)

    #     for i in [1, 2]:
    #         for j in [1, 2, 'combined']:
    #             for k in [1, 2]:
    #                 res1 = self.get_result_paths(abct1, i, j, k)
    #                 res2 = self.get_result_paths(abct3, i, j, k)
    #                 self.assertSameFiles([res1['sample'], res2['sample']])

    #     observed_paths = {
    #             1: dict(zip(range(1,3),
    #                     [self.get_test_path(prefix='obs1') for i in range(2)])),
    #             2: dict(zip(range(1,3),
    #                     [self.get_test_path(prefix='obs2') for i in range(2)]))}
    #     with open(obs_worker1.prior_stats_path, 'rU') as obs_stream:
    #         head = obs_stream.next()
    #         for i, line in enumerate(obs_stream):
    #             with open(observed_paths[1][i+1], 'w') as out:
    #                 out.write(head)
    #                 out.write(line)
    #     with open(obs_worker2.prior_stats_path, 'rU') as obs_stream:
    #         head = obs_stream.next()
    #         for i, line in enumerate(obs_stream):
    #             with open(observed_paths[2][i+1], 'w') as out:
    #                 out.write(head)
    #                 out.write(line)

    #     prior_paths = {
    #             1: glob.glob(os.path.join(abct2.summary_dir,
    #                     '*m1-prior-sample.txt'))[0],
    #             2: glob.glob(os.path.join(abct2.summary_dir,
    #                     '*m2-prior-sample.txt'))[0]}
    #     summary_paths = {
    #             1: glob.glob(os.path.join(abct2.summary_dir,
    #                     '*m1-stat-means-and-std-devs.txt'))[0],
    #             2: glob.glob(os.path.join(abct2.summary_dir,
    #                     '*m2-stat-means-and-std-devs.txt'))[0],
    #             'combined': glob.glob(os.path.join(abct2.summary_dir,
    #                     '*m12-combined-stat-means-and-std-devs.txt'))[0]}

    #     tmp_post_paths = {}
    #     post_paths = {}
    #     for i in [1, 2]: # observed_file
    #         tmp_post_paths[i] = {}
    #         post_paths[i] = {}
    #         for j in [1, 2, 'combined']: # model
    #             tmp_post_paths[i][j] = {}
    #             post_paths[i][j] = {}
    #             for k in [1, 2]: # simulated dataset
    #                 tmp_post_paths[i][j][k] = self.get_test_path(
    #                         prefix='tmp-post-d{0}-m{1}-s{2}-'.format(i,j,k))
    #                 post_paths[i][j][k] = self.get_test_path(
    #                         prefix='post-d{0}-m{1}-s{2}-'.format(i,j,k))
    #     for i in [1, 2]: # observed_file
    #         for j in [1, 2]: # model
    #             for k in [1, 2]: # simulated dataset
    #                 rej_worker = EuRejectWorker(
    #                         temp_fs = self.temp_fs,
    #                         observed_path = observed_paths[i][k],
    #                         prior_paths = [prior_paths[j]],
    #                         posterior_path = tmp_post_paths[i][j][k],
    #                         num_posterior_samples = num_posterior_samples,
    #                         num_standardizing_samples = 0,
    #                         summary_in_path = summary_paths[j])
    #                 rej_worker.start()

    #                 samples = parse_parameters(tmp_post_paths[i][j][k])
    #                 ipc = IntegerPartitionCollection(samples['taus'])
    #                 div_models_to_indices = {}
    #                 for idx, key in enumerate(ipc.iterkeys()):
    #                     div_models_to_indices[key] = idx + 1
    #                 add_div_model_column(tmp_post_paths[i][j][k], post_paths[i][j][k],
    #                         div_models_to_indices,
    #                         compresslevel = None)

    #                 res1 = self.get_result_paths(abct1, i, j, k)
    #                 res2 = self.get_result_paths(abct3, i, j, k)
    #                 self.assertSameFiles([res1['sample'], res2['sample'],
    #                         post_paths[i][j][k]])

    #     for i in [1, 2]: # observed_file
    #         for j in ['combined']: # model
    #             for k in [1, 2]:
    #                 rej_worker = EuRejectWorker(
    #                         temp_fs = self.temp_fs,
    #                         observed_path = observed_paths[i][k],
    #                         prior_paths = [tmp_post_paths[i][1][k], tmp_post_paths[i][2][k]],
    #                         posterior_path = tmp_post_paths[i][j][k],
    #                         num_posterior_samples = num_posterior_samples,
    #                         num_standardizing_samples = 0,
    #                         summary_in_path = summary_paths['combined'])
    #                 rej_worker.start()

    #                 samples = parse_parameters(tmp_post_paths[i][j][k])
    #                 ipc = IntegerPartitionCollection(samples['taus'])
    #                 div_models_to_indices = {}
    #                 for idx, key in enumerate(ipc.iterkeys()):
    #                     div_models_to_indices[key] = idx + 1
    #                 add_div_model_column(tmp_post_paths[i][j][k], post_paths[i][j][k],
    #                         div_models_to_indices,
    #                         compresslevel = None)

    #                 res1 = self.get_result_paths(abct1, i, j, k)
    #                 res2 = self.get_result_paths(abct3, i, j, k)
    #                 self.assertSameFiles([res1['sample'], res2['sample'],
    #                         post_paths[i][j][k]])

if __name__ == '__main__':
    unittest.main()

