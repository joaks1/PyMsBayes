#! /usr/bin/env python

import sys
import os
import re
import glob

from configobj import ConfigObj

from pymsbayes import config
from pymsbayes.utils import probability
from pymsbayes.utils import functions
from pymsbayes.fileio import process_file_arg, expand_path
from pymsbayes.utils import stats
from pymsbayes.utils import parsing
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

def get_partitions_from_posterior_sample_file(file_obj,
        integer_partitions = False):
    post = parsing.parse_parameters(file_obj)
    if not 'taus' in post:
        name = getattr(file_obj, 'name', file_obj)
        raise Exception('posterior sample in {0} does not contain a '
                'divergence time vector'.format(name))
    if integer_partitions:
        return stats.IntegerPartitionCollection(post['taus'])
    return stats.PartitionCollection(post['taus'])

def parse_omega_results_file(file_obj):
    s_iter = parsing.spreadsheet_iter([file_obj], sep = '\t')
    i = -1
    for i, d in enumerate(s_iter):
        pass
    if i != 0:
        raise Exception('too many lines in omega results file {0!r}'.format(
                file_obj))
    try:
        threshold = float(d['omega_thresh'])
        prob_less = float(d['prob_less_than'])
        prob_less_glm = float(d['glm_prob_less_than'])
    except Exception:
        _LOG.error('bad format of omega results file {0!r}'.format(
                file_obj))
        raise
    return {'threshold': threshold,
            'prob_less': prob_less,
            'prob_less_glm': prob_less_glm}

def parse_psi_results_file(file_obj):
    s_iter = parsing.spreadsheet_iter([file_obj], sep = '\t')
    results = {}
    for d in s_iter:
        try:
            psi = int(d['num_of_div_events'])
            prob = float(d['estimated_prob'])
            prob_glm = float(d['glm_adjusted_prob'])
        except Exception:
            _LOG.error('bad format of psi results file {0!r}'.format(
                    file_obj))
            raise
        results[psi] = {'prob': prob, 'prob_glm': prob_glm}
    return results

def parse_model_results_file(file_obj):
    s_iter = parsing.spreadsheet_iter([file_obj], sep = '\t')
    results = {}
    for d in s_iter:
        try:
            model = int(d['model'])
            prob = float(d['estimated_prob'])
            prob_glm = float(d['glm_adjusted_prob'])
        except Exception:
            _LOG.error('bad format of model results file {0!r}'.format(
                    file_obj))
            raise
        results[model] = {'prob': prob, 'prob_glm': prob_glm}
    return results

def parse_posterior_summary_file(file_obj):
    return ConfigObj(file_obj)

def parse_data_key_file(path):
    wd = os.path.dirname(path)
    f, close = process_file_arg(path)
    observed_paths = {}
    for line in f:
        l = line.strip().split('=')
        if len(l) != 2:
            raise Exception('unexpected line {0!r} in data key file'.format(
                line))
        observed_index = l[0].strip().strip('d')
        p = os.path.abspath(os.path.join(wd, l[1].strip()))
        observed_paths[int(observed_index)] = p
    f.close()
    return observed_paths

def parse_model_key_file(path):
    wd = os.path.dirname(path)
    f, close = process_file_arg(path)
    model_paths = {}
    for line in f:
        l = line.strip().split('=')
        if len(l) != 2:
            raise Exception('unexpected line {0!r} in model key file'.format(
                line))
        model_index = l[0].strip().strip('m')
        p = os.path.abspath(os.path.join(wd, l[1].strip()))
        model_paths[int(model_index)] = p
    f.close()
    return model_paths


class DMCSimulationResults(object):
    result_file_name_pattern = re.compile(r'^d(?P<observed_index>\d+)-'
            'm(?P<prior_index>\d+)-s(?P<sim_index>\d+)-(?P<result_index>\d+)-'
            '(?P<contents>[a-zA-Z-]+).(?P<extension>[a-zA-Z.]+)$')

    def __init__(self, info_path):
        self.info_path = expand_path(info_path)
        self.results_dir = os.path.dirname(self.info_path)
        self.output_dir = os.path.join(self.results_dir, 'pymsbayes-output')
        self.observed_index_to_path = {}
        self.observed_path_to_index = {}
        self.observed_index_to_config = {}
        self.observed_config_to_index = {}
        self.observed_index_to_prior_index = {}
        self.prior_index_to_config = {}
        self.prior_config_to_index = {}
        self.num_taxon_pairs = None
        self.num_sim_reps = None
        self.sort_index = None
        self._parse_info_file()
        self.combined_prior_index = '{0}-combined'.format(''.join(
                [str(i) for i in sorted(self.prior_index_to_config.keys())]))
        self.combined_prior_results = True
        for i in self.observed_index_to_path.iterkeys():
            prefix = self.get_result_path_prefix(i, self.combined_prior_index,
                    1)
            pattern = prefix + '*-posterior-sample*'
            files = glob.glob(pattern)
            if not files:
                self.combined_prior_results = False
                break
        self.prior_configs = {}
        for k, v in self.prior_index_to_config.iteritems():
            self.prior_configs[k] = config.MsBayesConfig(v)
        self.compresslevel = 9
    
    def _parse_info_file(self):
        c = ConfigObj(self.info_path)
        settings = c.get('pymsbayes', None)
        if not settings:
            raise Exception('{0!r} is not a valid pymsbayes info file'.format(
                    self.info_path))
        if not set(['observed_configs', 'observed_paths',
                'prior_configs']).issubset(
                        set(settings.keys())):
            raise Exception('info file {0!r} is missing path information')
        for k in settings['observed_configs'].iterkeys():
            self.observed_index_to_config[int(k)] = os.path.abspath(
                    os.path.join(self.results_dir,
                            settings['observed_configs'][k]))
            self.observed_index_to_path[int(k)] = os.path.abspath(os.path.join(
                    self.results_dir,
                    settings['observed_paths'][k]))
        for k, v in settings['prior_configs'].iteritems():
            self.prior_index_to_config[int(k)] = os.path.abspath(os.path.join(
                    self.results_dir, v))
        omx = max(self.observed_index_to_path.iterkeys())
        pmx = max(self.prior_index_to_config.iterkeys())
        assert sorted(self.observed_index_to_path.keys()) == list(range(1,
                omx + 1))
        assert sorted(self.prior_index_to_config.keys()) == list(range(1,
                pmx + 1))
        self.observed_path_to_index = dict(zip([self.observed_index_to_path[
                k] for k in self.observed_index_to_path.iterkeys()],
                self.observed_index_to_path.iterkeys()))
        self.observed_config_to_index = dict(zip([self.observed_index_to_config[
                k] for k in self.observed_index_to_config.iterkeys()],
                self.observed_index_to_config.iterkeys()))
        self.prior_config_to_index = dict(zip([self.prior_index_to_config[
                k] for k in self.prior_index_to_config.iterkeys()],
                self.prior_index_to_config.iterkeys()))
        for i, obs_cfg in self.observed_index_to_config.iteritems():
            self.observed_index_to_prior_index[i] = \
                    self.prior_config_to_index.get(obs_cfg, -1)
        self.final_result_index = 1
        result_indices = self.get_result_indices(1, 1, 1)
        if result_indices:
            self.final_result_index = max(result_indices)
        self.num_taxon_pairs = int(settings['num_taxon_pairs'])
        self.num_sim_reps = int(settings['simulation_reps'])
        self.sort_index = int(settings['sort_index'])

    def get_result_dir(self, observed_index, prior_index):
        return os.path.join(self.output_dir, 'd' + str(observed_index),
                'm' + str(prior_index))

    def get_result_summary_path(self, observed_index, prior_index):
        fname = 'results.txt'
        if self.compresslevel:
            fname += '.gz'
        return os.path.join(self.get_result_dir(observed_index, prior_index),
                fname)

    def get_result_file_prefix(self, observed_index, prior_index, sim_index):
        return 'd{0}-m{1}-s{2}-'.format(observed_index, prior_index, sim_index)

    def get_result_path_prefix(self, observed_index, prior_index, sim_index):
        return os.path.join(self.get_result_dir(observed_index, prior_index),
                self.get_result_file_prefix(observed_index, prior_index,
                        sim_index))

    def get_result_indices(self, observed_index, prior_index, sim_index):
        prefix = self.get_result_path_prefix(observed_index, prior_index,
                sim_index)
        pattern = prefix + '*-posterior-*'
        result_files = [os.path.basename(x) for x in glob.glob(pattern)]
        result_indices = set()
        for f in result_files:
            m = self.result_file_name_pattern.match(f)
            result_indices.add(int(m.group('result_index')))
        return sorted(result_indices) 

    def result_iter(self, observed_index, prior_index):
        path_iter = self.result_path_iter(observed_index, prior_index)
        for i, (true_params, paths) in enumerate(path_iter):
            yield self.get_results_from_params_and_result_paths(
                    true_params, paths)

    def get_results_from_params_and_result_paths(self, true_params, paths):
        if not os.path.exists(paths['summary']):
            raise Exception('Posterior summary file {0!r} does not '
                    'exist'.format(paths['summary']))
        summary = parse_posterior_summary_file(paths['summary'])
        results = {}
        tau_true = float(true_params['PRI.E.t'])
        try:
            tau_mode_min = float(summary['PRI.E.t']['modes'][0].strip('()'))
            tau_mode_max = float(summary['PRI.E.t']['modes'][1].strip('()'))
            tau_median = float(summary['PRI.E.t']['median'])
            tau_mode_glm = float(summary['PRI.E.t']['mode_glm'])
        except:
            _LOG.error('Problem extracting "PRI.E.t" info from posterior '
                    'summary file {0!r'.format(paths['summary']))
            raise
        tau_mode = (tau_mode_min + tau_mode_max) / float(2)
        results['mean_tau'] = {'true': tau_true,
                'mode': tau_mode,
                'median': tau_median,
                'mode_glm': tau_mode_glm}
        omega_true = float(true_params['PRI.omega'])
        try:
            omega_mode_min = float(summary['PRI.omega']['modes'][0].strip('()'))
            omega_mode_max = float(summary['PRI.omega']['modes'][1].strip('()'))
            omega_median = float(summary['PRI.omega']['median'])
            omega_mode_glm = float(summary['PRI.omega']['mode_glm'])
        except:
            _LOG.error('Problem extracting "PRI.omega" info from posterior '
                    'summary file {0!r'.format(paths['summary']))
            raise
        omega_mode = (omega_mode_min + omega_mode_max) / float(2)
        omega_results = parse_omega_results_file(paths['omega'])
        results['omega'] = {'true': omega_true,
                'mode': omega_mode,
                'median': omega_median,
                'mode_glm': omega_mode_glm}
        results['omega'].update(omega_results)
        psi_true = int(true_params['PRI.Psi'])
        try:
            psi_mode = summary['PRI.Psi']['modes']
        except:
            _LOG.error('Problem extracting "PRI.Psi.modes" info from posterior '
                    'summary file {0!r'.format(paths['summary']))
            raise
        if not isinstance(psi_mode, str):
            psi_mode = psi_mode[0]
        psi_mode = int(psi_mode)
        try:
            psi_mode_glm = float(summary['PRI.Psi']['mode_glm'])
        except:
            _LOG.error('Problem extracting "PRI.Psi.mode_glm" info from '
                    'posterior summary file {0!r'.format(paths['summary']))
            raise
        psi_results = parse_psi_results_file(paths['psi'])
        results['psi'] = {'true': psi_true,
                'mode': psi_mode,
                'mode_glm': psi_mode_glm,
                'probs': psi_results}
        model_true = int(true_params['PRI.model'])
        try:
            model_mode = summary['PRI.model']['modes']
        except:
            _LOG.error('Problem extracting "PRI.model.modes" info from '
                    'posterior summary file {0!r'.format(paths['summary']))
            raise
        if not isinstance(model_mode, str):
            model_mode = model_mode[0]
        model_mode = int(model_mode)
        try:
            model_mode_glm = float(summary['PRI.model']['mode_glm'])
        except:
            _LOG.error('Problem extracting "PRI.model.mode_glm" info from '
                    'posterior summary file {0!r'.format(paths['summary']))
            raise
        model_results = parse_model_results_file(paths['model'])
        results['model'] = {'true': model_true,
                'mode': model_mode,
                'mode_glm': model_mode_glm,
                'probs': model_results}
        for k, v in true_params.iteritems():
            if (parsing.TAU_PATTERNS[0].match(k) or
                    parsing.D_THETA_PATTERNS[0].match(k) or
                    parsing.A_THETA_PATTERNS[0].match(k)):
                results[k] = float(v)
        return results

    def result_to_flat_dict(self, result):
        d = {'mean_tau_true': result['mean_tau']['true'],
             'mean_tau_mode': result['mean_tau']['mode'],
             'mean_tau_median': result['mean_tau']['median'],
             'mean_tau_mode_glm': result['mean_tau']['mode_glm'],
             'omega_true': result['omega']['true'],
             'omega_mode': result['omega']['mode'],
             'omega_median': result['omega']['median'],
             'omega_mode_glm': result['omega']['mode_glm'],
             'omega_threshold': result['omega']['threshold'],
             'omega_prob_less': result['omega']['prob_less'],
             'omega_prob_less_glm': result['omega']['prob_less_glm'],
             'psi_true': result['psi']['true'],
             'psi_mode': result['psi']['mode'],
             'psi_mode_glm': result['psi']['mode_glm'],
             'model_true': result['model']['true'],
             'model_mode': result['model']['mode'],
             'model_mode_glm': result['model']['mode_glm']}
        for i in result['psi']['probs'].iterkeys():
            d['psi_{0}_prob'.format(i)] = result['psi']['probs'][i]['prob']
            d['psi_{0}_prob_glm'.format(i)] = result['psi']['probs'][i][
                    'prob_glm']
        for i in result['model']['probs'].iterkeys():
            d['model_{0}_prob'.format(i)] = result['model']['probs'][i]['prob']
            d['model_{0}_prob_glm'.format(i)] = result['model']['probs'][i][
                    'prob_glm']
        for k, v in result.iteritems():
            if not isinstance(v, dict):
                d[k] = v
        return d

    def flat_result_iter(self, observed_index, prior_index,
            include_tau_exclusion_info = False):
        glm_failures = []
        for sim_idx, result in enumerate(self.result_iter(observed_index,
                prior_index)):
            r = self.result_to_flat_dict(result)
            if include_tau_exclusion_info:
                div_times = sorted([r['PRI.t.' + str(i)] for i in range(1,
                        self.num_taxon_pairs + 1)])
                model_index = r['model_mode']
                try:
                    model_index_glm = int(round(r['model_mode_glm']))
                except ValueError:
                    model_index_glm = model_index
                    glm_failures.append(sim_idx)
                tau_max = self.prior_configs[model_index].tau.maximum 
                tau_max_glm = self.prior_configs[model_index_glm].tau.maximum
                prob_of_exclusion = 0.0
                prob_of_exclusion_glm = 0.0
                prior_prob_of_exclusion = 0.0
                model_prior = float(1) / len(self.prior_configs.keys())
                bf_tau_max = []
                for i in self.prior_configs.iterkeys():
                    if max(div_times) > self.prior_configs[i].tau.maximum:
                        prob_of_exclusion += r['model_{0}_prob'.format(i)]
                        prob_of_exclusion_glm += r['model_{0}_prob_glm'.format(i)]
                        prior_prob_of_exclusion += model_prior
                        bf_tau_max.append(self.prior_configs[i].tau.maximum)
                if len(bf_tau_max) < 1:
                    bf_tau_max = float('inf')
                else:
                    bf_tau_max = max(bf_tau_max)
                ex = functions.get_sublist_greater_than(div_times, tau_max)
                ex_glm = functions.get_sublist_greater_than(div_times,
                        tau_max_glm)
                prior_odds = (prior_prob_of_exclusion /
                        (1 - prior_prob_of_exclusion))
                post_odds = (prob_of_exclusion / (1 - prob_of_exclusion))
                post_odds_glm = (prob_of_exclusion_glm /
                        (1 - prob_of_exclusion_glm))
                if probability.almost_equal(prior_odds, 0.0):
                    bf_of_exclusion = float('inf')
                    bf_of_exclusion_glm = float('inf')
                    if probability.almost_equal(bf_of_exclusion,
                            0.0):
                        bf_of_exclusion = 0.0
                    if probability.almost_equal(
                            bf_of_exclusion_glm, 0.0):
                        bf_of_exclusion_glm = 0.0
                else:
                    bf_of_exclusion = post_odds / prior_odds
                    bf_of_exclusion_glm = post_odds_glm / prior_odds
                bf_ex = []
                bf_ex_glm = []
                if bf_of_exclusion > 1.0:
                    bf_ex = functions.get_sublist_greater_than(div_times,
                            bf_tau_max)
                if bf_of_exclusion_glm > 1.0:
                    bf_ex_glm = functions.get_sublist_greater_than(div_times,
                            bf_tau_max)
                r['prob_of_exclusion'] = prob_of_exclusion
                r['prob_of_exclusion_glm'] = prob_of_exclusion_glm
                r['prior_prob_of_exclusion'] = prior_prob_of_exclusion
                r['bf_of_exclusion'] = bf_of_exclusion
                r['bf_of_exclusion_glm'] = bf_of_exclusion_glm
                r['tau_max'] = tau_max
                r['tau_max_glm'] = tau_max_glm
                r['num_excluded'] = len(ex)
                r['num_excluded_glm'] = len(ex_glm)
                r['bf_num_excluded'] = len(bf_ex)
                r['bf_num_excluded_glm'] = len(bf_ex_glm)
            yield r
        if len(glm_failures) > 0:
            _LOG.warning('WARNING: there were GLM-regression failures:\n'
                    'For observed index {0} prior index {1}, there were '
                    'failures at the following simulation indices:\n'
                    '{2}'.format(observed_index, prior_index, glm_failures))

    def result_path_iter(self, observed_index, prior_index):
        true_model = self.observed_index_to_prior_index[observed_index]
        out_dir = self.get_result_dir(observed_index, prior_index)
        if not os.path.isdir(out_dir):
            raise Exception('expected result direcory {0!r} does not '
                    'exist'.format(out_dir))
        observed_stream, close = process_file_arg(
                self.observed_index_to_path[observed_index])
        header = parsing.parse_header(observed_stream, sep = '\t', strict = True,
                seek = False)
        parameter_indices = functions.get_indices_of_patterns(header,
                parsing.PARAMETER_PATTERNS)
        for i, line in enumerate(observed_stream):
            l = line.strip().split()
            true_params = dict(zip([header[x] for x in parameter_indices],
                    [l[x] for x in parameter_indices]))
            true_params['PRI.model'] = str(true_model)
            result_prefix = '{0}{1}-'.format(self.get_result_path_prefix(
                    observed_index, prior_index, i + 1), 
                    self.final_result_index)
            summary_path = result_prefix + 'posterior-summary.txt'
            psi_path = result_prefix + 'psi-results.txt'
            omega_path = result_prefix + 'omega-results.txt'
            div_model_path = result_prefix + 'div-model-results.txt'
            model_path = result_prefix + 'model-results.txt'
            paths = {'summary': summary_path,
                     'psi': psi_path,
                     'omega': omega_path,
                     'div-model': div_model_path,
                     'model': model_path}
            yield true_params, paths
        observed_stream.close()

    def write_result_summaries(self, prior_indices = None, sep = '\t',
            include_tau_exclusion_info = False):
        if not prior_indices:
            prior_indices = self.prior_index_to_config.keys()
            if self.combined_prior_index:
                prior_indices.append(self.combined_prior_index)
        for prior_idx in prior_indices:
            for observed_idx in self.observed_index_to_path.iterkeys():
                out_path = self.get_result_summary_path(observed_idx, prior_idx)
                out_path = functions.get_new_path(out_path)
                out, close = process_file_arg(out_path, 'w',
                        compresslevel = self.compresslevel)
                keys = []
                for i, r in enumerate(self.flat_result_iter(observed_idx,
                        prior_idx, include_tau_exclusion_info)):
                    if i == 0:
                        keys = r.keys()
                        out.write('{0}\n'.format(sep.join(keys)))
                    out.write('{0}\n'.format(sep.join([str(r[k
                            ]) for k in keys])))
                out.close()
        

class UnorderedDivergenceModelResults(object):
    def __init__(self, div_model_results_path,
            inclusion_threshold = None):
        self.path = div_model_results_path
        self.inclusion_threshold = inclusion_threshold
        self.models = []
        self.n = 0
        self.cumulative_prob = 0.0
        header = parsing.parse_header(self.path)
        if 'div_model_with_conditional_age_estimates' in header:
            self._parse_results_file()
        else:
            self._parse_posterior_file()

    def _full(self):
        if not self.inclusion_threshold:
            return False
        if ((isinstance(self.inclusion_threshold, float)) and
                (self.cumulative_prob >= self.inclusion_threshold)):
            return True
        if ((isinstance(self.inclusion_threshold, int)) and
                (self.n >= self.inclusion_threshold)):
            return True
        return False

    def _parse_results_file(self):
        file_stream, close = process_file_arg(self.path)
        ss_iter = parsing.spreadsheet_iter([file_stream])
        for d in ss_iter:
            if self._full():
                if close:
                    file_stream.close()
                return
            try:
                dms = UnorderedDivergenceModelSummary(d)
            except:
                file_stream.close()
                raise
            self.n += 1
            self.cumulative_prob += dms.prob
            self.models.append(dms)
        if close:
            file_stream.close()

    def _parse_posterior_file(self):
        div_models = get_partitions_from_posterior_sample_file(
                self.path, integer_partitions = True)
        for k, m in div_models.iteritems():
            if self._full():
                return
            dms = UnorderedDivergenceModelSummary()
            dms.int_partition = m.integer_partition
            dms.prob = div_models.get_frequency(k)
            dms.age_info = [(i, s) for i, s in m.iter_item_summaries()]
            self.n += 1
            self.cumulative_prob += dms.prob
            self.models.append(dms)

class UnorderedDivergenceModelSummary(object):
    number_pattern_string = r'[\d\.Ee\-\+]+'
    age_info_pattern_string = (
            '(?P<index>\d+)\:{0}\[&age_median='
            '(?P<median>{0}),age_mean='
            '(?P<mean>{0}),age_n='
            '(?P<n>{0}),age_range={{'
            '(?P<range1>{0}),'
            '(?P<range2>{0})}},age_hpdi_95={{'
            '(?P<hpdi1>{0}),'
            '(?P<hpdi2>{0})}},age_qi_95={{'
            '(?P<qi1>{0}),'
            '(?P<qi2>{0})}}\]'.format(number_pattern_string))
    age_info_pattern = re.compile(age_info_pattern_string)

    def __init__(self, div_model_results_file_line_dict = None):
        self.int_partition = None
        self.age_info = None
        self.prob = None
        self.glm_prob = None
        if div_model_results_file_line_dict:
            self._parse_line_dict(div_model_results_file_line_dict)

    def _parse_line_dict(self, line_dict):
        self.prob = float(line_dict['estimated_prob'])
        self.glm_prob = float(line_dict['glm_adjusted_prob'])
        self._parse_model(model_key_string = line_dict['divergence_model'],
                model_summary_string = line_dict[
                        'div_model_with_conditional_age_estimates'])

    def _parse_model(self, model_key_string, model_summary_string):
        model_key = [int(x) for x in model_key_string.split(',')]
        matches = []
        if model_key[0] == 0:
            raise Exception('Divergence models appear ordered')
        self.int_partition = model_key
        for m in self.age_info_pattern.finditer(model_summary_string):
            d = {'index': int(m.group('index')),
                    'median': float(m.group('median')),
                    'mean': float(m.group('mean')),
                    'n': int(m.group('n')),
                    'range': (float(m.group('range1')),
                              float(m.group('range2'))),
                    'hpdi_95': (float(m.group('hpdi1')),
                             float(m.group('hpdi2'))),
                    'qi_95': (float(m.group('qi1')), float(m.group('qi2'))),
                    }
            matches.append(d)
        self.age_info = []
        for m in matches:
            index = m.pop('index')
            self.age_info.append((index, m))
        assert sorted(self.int_partition) == sorted(
                [k for k, d in self.age_info])

    def iter_divergences(self):
        for k, d in self.age_info:
            yield k, d

class OrderedDivergenceModelResults(object):
    def __init__(self, div_model_results_path,
            inclusion_threshold = None):
        self.path = div_model_results_path
        self.inclusion_threshold = inclusion_threshold
        self.models = []
        self.n = 0
        self.cumulative_prob = 0.0
        self.npairs = None
        header = parsing.parse_header(self.path)
        if 'div_model_with_conditional_age_estimates' in header:
            self._parse_results_file()
        else:
            self._parse_posterior_file()

    def _full(self):
        if not self.inclusion_threshold:
            return False
        if ((isinstance(self.inclusion_threshold, float)) and
                (self.cumulative_prob >= self.inclusion_threshold)):
            return True
        if ((isinstance(self.inclusion_threshold, int)) and
                (self.n >= self.inclusion_threshold)):
            return True
        return False

    def _parse_results_file(self):
        file_stream, close = process_file_arg(self.path)
        ss_iter = parsing.spreadsheet_iter([file_stream])
        for d in ss_iter:
            if self._full():
                if close:
                    file_stream.close()
                return
            try:
                dms = OrderedDivergenceModelSummary(d)
            except:
                file_stream.close()
                raise
            self.n += 1
            self.cumulative_prob += dms.prob
            self.models.append(dms)
            if not self.npairs:
                self.npairs = dms.npairs
            else:
                assert self.npairs == dms.npairs
        if close:
            file_stream.close()

    def _parse_posterior_file(self):
        div_models = get_partitions_from_posterior_sample_file(self.path)
        for k, m in div_models.iteritems():
            if self._full():
                return
            dms = OrderedDivergenceModelSummary()
            dms.prob = div_models.get_frequency(k)
            dms._parse_model(
                    model_key_string = k,
                    model_summary_string = m.value_summary_string())
            self.n += 1
            self.cumulative_prob += dms.prob
            self.models.append(dms)
            if not self.npairs:
                self.npairs = dms.npairs
            else:
                assert self.npairs == dms.npairs

class OrderedDivergenceModelCollection(OrderedDivergenceModelResults):
    def __init__(self, div_model_results_path):
        OrderedDivergenceModelResults.__init__(self, div_model_results_path)
        assert probability.almost_equal(self.cumulative_prob,
                1.0)

    def prob_of_shared_divergence(self, taxon_indices):
        taxon_indices = list(taxon_indices)
        for i in taxon_indices:
            if ((i < 0) or (i >= self.npairs)):
                raise ValueError('taxon index {0} is out of bounds'.format(i))
        prob_shared = 0.0
        for m in self.models:
            div_indices = [d for i, d in enumerate(
                    m.partition) if i in taxon_indices]
            if len(set(div_indices)) == 1:
                prob_shared += m.prob
        return prob_shared


class OrderedDivergenceModelSummary(object):
    number_pattern_string = r'[\d\.Ee\-\+]+'
    age_info_pattern_string = (
            '(?P<index>\d+)\:{0}\[&age_median='
            '(?P<median>{0}),age_mean='
            '(?P<mean>{0}),age_n='
            '(?P<n>{0}),age_range={{'
            '(?P<range1>{0}),'
            '(?P<range2>{0})}},age_hpdi_95={{'
            '(?P<hpdi1>{0}),'
            '(?P<hpdi2>{0})}},age_qi_95={{'
            '(?P<qi1>{0}),'
            '(?P<qi2>{0})}}\]'.format(number_pattern_string))
    age_info_pattern = re.compile(age_info_pattern_string)

    def __init__(self, div_model_results_file_line_dict = None):
        self.partition = None
        self.age_info = None
        self.prob = None
        self.glm_prob = None
        self.npairs = None
        if div_model_results_file_line_dict:
            self._parse_line_dict(div_model_results_file_line_dict)

    def _parse_line_dict(self, line_dict):
        self.prob = float(line_dict['estimated_prob'])
        self.glm_prob = float(line_dict['glm_adjusted_prob'])
        self._parse_model(model_key_string = line_dict['divergence_model'],
                model_summary_string = line_dict[
                        'div_model_with_conditional_age_estimates'])

    def _parse_model(self, model_key_string, model_summary_string):
        model_key = [int(x) for x in model_key_string.split(',')]
        matches = []
        if model_key[0] != 0:
            raise Exception('Divergence models appear unordered')
        self.partition = model_key
        self.npairs = len(self.partition)
        self.divergence_indices = sorted(set(self.partition))
        for m in self.age_info_pattern.finditer(model_summary_string):
            d = {'index': int(m.group('index')),
                    'median': float(m.group('median')),
                    'mean': float(m.group('mean')),
                    'n': int(m.group('n')),
                    'range': (float(m.group('range1')),
                              float(m.group('range2'))),
                    'hpdi_95': (float(m.group('hpdi1')),
                             float(m.group('hpdi2'))),
                    'qi_95': (float(m.group('qi1')), float(m.group('qi2'))),
                    }
            matches.append(d)
        self.age_info = {}
        for m in matches:
            index = m.pop('index')
            self.age_info[index] = m
        assert self.divergence_indices == sorted(self.age_info.keys())

    def iter_divergences(self):
        for k in self.divergence_indices:
            d = self.age_info[k]
            taxon_indices = [i for i, x in enumerate(self.partition) if x == k]
            yield taxon_indices, d

    def iter_per_element_divergences(self):
        for k in self.partition:
            d = self.age_info[k]
            yield k, d

