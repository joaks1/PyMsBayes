#! /usr/bin/env python

import sys
import os
import re
import glob
import math

from configobj import ConfigObj

import pymsbayes
from pymsbayes.fileio import process_file_arg, expand_path
from pymsbayes.config import MsBayesConfig
from pymsbayes.utils import MSBAYES_SORT_INDEX
from pymsbayes.utils.functions import get_sublist_greater_than, get_new_path
from pymsbayes.utils.probability import almost_equal
from pymsbayes.utils.errors import (SummaryFileParsingError,
        ParameterParsingError)
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

GZIP_FILE_PATTERN = re.compile(r'.*\.gz$', re.IGNORECASE)

def get_indices_of_patterns(target_list, regex_list, sort=True):
    indices = []
    for regex in regex_list:
        indices.extend([i for i, e in enumerate(target_list) if regex.match(e)])
    if sort:
        return sorted(indices)
    return indices

def get_indices_of_strings(target_list, string_list, sort=True):
    indices = []
    for s in string_list:
        indices.extend([i for i, e in enumerate(target_list) if s.strip() == e.strip()])
    if sort:
        return sorted(indices)
    return indices
    
def reduce_columns(in_file, out_file, column_indices, sep='\t',
        extra_tab=False):
    if not column_indices:
        raise Exception('no column indices to retain')
    in_stream, close_in = process_file_arg(in_file, 'rU')
    out_stream, close_out = process_file_arg(out_file, 'w')
    line_iter = iter(in_stream)
    for line_num, line in enumerate(line_iter):
        l = line.strip().split(sep)
        new_line = [l[i] for i in column_indices]
        if extra_tab:
            out_stream.write('%s\t\n' % sep.join(new_line))
        else:
            out_stream.write('%s\n' % sep.join(new_line))
    if close_in:
        in_stream.close()
    if close_out:
        out_stream.close()

##############################################################################
## patterns and functions for manipulating msBayes priors

HEADER_PATTERN = re.compile(r'^\s*\D.+')
PARAMETER_PATTERNS = [
        re.compile(r'^\s*PRI\.(?!numTauClass)\S+\s*$'),
        ]
DEFAULT_STAT_PATTERNS = [
        re.compile(r'^\s*pi\.\d+\s*$'),
        re.compile(r'^\s*wattTheta\.\d+\s*$'),
        re.compile(r'^\s*pi\.net\.\d+\s*$'),
        re.compile(r'^\s*tajD\.denom\.\d+\s*$'),
        ]
ALL_STAT_PATTERNS = [
        re.compile(r'^\s*(?!PRI)\S+\s*$'),
        ]
DUMMY_PATTERNS = [
        re.compile(r'^\s*PRI\.numTauClass\s*$')
        ]
MODEL_PATTERNS = [
        re.compile(r'^\s*PRI\.model\s*$'),
        ]
TAU_PATTERNS = [
        re.compile(r'^\s*PRI\.t\.\d+\s*$'),
        ]
D_THETA_PATTERNS = [
        re.compile(r'^\s*PRI\.d[12]Theta\.\d+\s*$'),
        ]
D1_THETA_PATTERNS = [
        re.compile(r'^\s*PRI\.d1Theta\.\d+\s*$'),
        ]
D2_THETA_PATTERNS = [
        re.compile(r'^\s*PRI\.d2Theta\.\d+\s*$'),
        ]
A_THETA_PATTERNS = [
        re.compile(r'^\s*PRI\.aTheta\.\d+\s*$'),
        ]
PSI_PATTERNS = [
        re.compile(r'^\s*PRI\.Psi\s*$'),
        ]
MEAN_TAU_PATTERNS = [
        re.compile(r'^\s*PRI\.E\.t\s*$'),
        ]
OMEGA_PATTERNS = [
        re.compile(r'^\s*PRI\.omega\s*$'),
        ]
DIV_MODEL_PATTERNS = [
        re.compile(r'^\s*PRI\.div\.model\s*$'),
        ]

def line_count(file_obj, ignore_headers=False):
    f, close = process_file_arg(file_obj)
    count = 0
    for line in f:
        if ignore_headers:
            if HEADER_PATTERN.match(line):
                continue
        count += 1
    if close:
        f.close()
    return count

def get_patterns_from_prefixes(prefixes, ignore_case=True):
    patterns = []
    for prefix in prefixes:
        pattern_str = r'^\s*{0}\d+\s*$'.format(prefix.replace('.', '\.'))
        if ignore_case:
            patterns.append(re.compile(pattern_str, re.IGNORECASE))
        else:
            patterns.append(re.compile(pattern_str))
    return patterns

def parse_header(file_obj, sep='\t', strict=True, seek=True):
    file_stream, close = process_file_arg(file_obj, 'rU')
    try:
        header_line = file_stream.next()
    except StopIteration:
        file_stream.close()
        if strict:
            raise Exception('did not find header in {0}'.format(file_stream.name))
        else:
            return None
    if not HEADER_PATTERN.match(header_line):
        file_stream.close()
        if strict:
            raise Exception('did not find header in {0}'.format(file_stream.name))
        else:
            return None
    header = header_line.strip().split(sep)
    if close:
        file_stream.close()
    else:
        if seek:
            file_stream.seek(0)
    return header

def spreadsheet_iter(spreadsheets, sep = '\t', header = None):
    head_line = False
    if not header:
        head_line = True
        header = parse_header(spreadsheets[0], sep = sep)
    d = dict(zip(header, [[] for i in range(len(header))]))
    for sheet_idx, ss in enumerate(spreadsheets):
        file_stream, close = process_file_arg(ss, 'rU')
        if head_line:
            h = file_stream.next().strip().split(sep)
            if header != h:
                raise Exception('headers do not match')
        for row_idx, row in enumerate(file_stream):
            if row.strip() == '':
                continue
            r = [el.strip() for el in row.strip().split(sep)]
            if len(r) != len(header):
                raise Exception('row {0} of spreadsheet {1} has {2} columns, '
                        'header has {3}'.format(row_idx + 1, sheet_idx + 1,
                                len(r), len(header)))
            yield dict(zip(header, r))
        if close:
            file_stream.close()

def dict_line_iter(d, sep = '\t', header = None):
    if not header:
        header = sorted(d.iterkeys())
    if sorted(header) != sorted(d.iterkeys()):
        raise ValueError('header does not match dict keys')
    yield '{0}\n'.format(sep.join(header))
    for i in range(len(d[header[0]])):
        yield '{0}\n'.format(sep.join([str(d[h][i]) for h in header]))

def get_dict_from_spreadsheets(spreadsheets, sep = '\t', header = None):
    ss_iter = spreadsheet_iter(spreadsheets, sep = sep, header = header)
    row_dict = ss_iter.next()
    d = dict(zip(row_dict.iterkeys(),
            [[row_dict[k]] for k in row_dict.iterkeys()]))
    for row_dict in ss_iter:
        for k in row_dict.iterkeys():
            d[k].append(row_dict[k])
    return d

def get_stats_by_time(spreadsheets, sep = '\t', header = None):
    tau_pattern = re.compile(r'^\s*(?P<prefix>PRI\.t)\.(?P<index>\d+)\s*$')
    stat_pattern = re.compile(r'^\s*(?!PRI)(?P<prefix>\S+)\.(?P<index>\d+)\s*$')
    tau_keys = []
    ss_iter = spreadsheet_iter(spreadsheets, sep = sep, header = header)
    row_dict = ss_iter.next()
    tau_keys = [k for k in row_dict.iterkeys() if tau_pattern.match(k)]
    stat_keys = [k for k in row_dict.iterkeys() if stat_pattern.match(k)]
    taxon_indices = sorted([int(k.split('.')[-1]) for k in tau_keys])
    ntaxa = max(taxon_indices)
    if taxon_indices != range(1, ntaxa + 1):
        raise Exception('unexpected taxon indices for tau parameters:\n\t'
                '{0}'.format(', '.join([str(i) for i in taxon_indices])))
    tau_prefixes = set()
    tau_prefixes.update([tau_pattern.match(k).group(
            'prefix') for k in tau_keys])
    if len(tau_prefixes) != 1:
        raise Exception('unexpected prefixes for tau parameters:\n\t'
                '{0}'.format(', '.join(tau_prefixes)))
    stat_prefixes = set()
    stat_prefixes.update([stat_pattern.match(k).group(
            'prefix') for k in stat_keys])
    for prefix in stat_prefixes:
        pattern = get_patterns_from_prefixes([prefix + '.'])[0]
        indices = sorted([int(
                k.split('.')[-1]) for k in stat_keys if pattern.match(k)])
        if indices != taxon_indices:
            raise Exception('unexpected taxon indices for stat {0!r}:\n\t'
                    '{1}'.format(prefix, ', '.join([str(i) for i in indices])))
    keys = list(tau_prefixes) + list(stat_prefixes)
    d = dict(zip(keys, [[] for i in range(len(keys))]))
    for k in d.iterkeys():
        for i in range(1, ntaxa + 1):
            d[k].append(float(row_dict['.'.join([k, str(i)])]))
    for row_dict in ss_iter:
        for k in d.iterkeys():
            for i in range(1, ntaxa + 1):
                d[k].append(float(row_dict['.'.join([k, str(i)])]))
    return d

def get_parameter_indices(header_list, parameter_patterns=PARAMETER_PATTERNS):
    return get_indices_of_patterns(header_list, parameter_patterns)

def get_stat_indices(header_list, stat_patterns=DEFAULT_STAT_PATTERNS):
    return get_indices_of_patterns(header_list, stat_patterns)

def get_dummy_indices(header_list, dummy_patterns=DUMMY_PATTERNS):
    return get_indices_of_patterns(header_list, dummy_patterns)
    
def observed_stats_for_abctoolbox(in_file, out_file,
        stat_patterns=DEFAULT_STAT_PATTERNS):
    header = parse_header(in_file)
    indices = get_stat_indices(header, stat_patterns=stat_patterns)
    reduce_columns(in_file, out_file, indices)
    return [header[i] for i in sorted(indices)]

def observed_parameters_for_abctoolbox(in_file, out_file,
        parameter_patterns=PARAMETER_PATTERNS):
    header = parse_header(in_file)
    indices = get_parameter_indices(header,
            parameter_patterns=parameter_patterns)
    reduce_columns(in_file, out_file, indices)
    return [header[i] for i in sorted(indices)]

def prior_for_abctoolbox(in_file, out_file,
        stat_patterns=DEFAULT_STAT_PATTERNS,
        parameter_patterns=PARAMETER_PATTERNS):
    header = parse_header(in_file)
    indices = get_parameter_indices(header,
            parameter_patterns=parameter_patterns)
    indices.extend(get_stat_indices(header, stat_patterns=stat_patterns))
    reduce_columns(in_file, out_file, sorted(indices), extra_tab=False)
    return [header[i] for i in sorted(indices)]

def prior_for_msreject(in_file, out_file,
        stat_patterns=DEFAULT_STAT_PATTERNS,
        parameter_patterns=PARAMETER_PATTERNS,
        dummy_patterns=DUMMY_PATTERNS,
        include_header=False):
    header = parse_header(in_file)
    in_file, close = process_file_arg(in_file)
    indices = get_parameter_indices(header,
            parameter_patterns=parameter_patterns)
    indices.extend(get_stat_indices(header, stat_patterns=stat_patterns))
    indices.extend(get_dummy_indices(header, dummy_patterns=DUMMY_PATTERNS))
    if not include_header:
        in_file.next()
    reduce_columns(in_file, out_file, sorted(indices), extra_tab=False)
    if close:
        in_file.close()
    return [header[i] for i in sorted(indices)]

def parameter_density_iter(parameter_density_file,
        parameter_patterns = DIV_MODEL_PATTERNS + MODEL_PATTERNS + \
                PSI_PATTERNS + MEAN_TAU_PATTERNS + OMEGA_PATTERNS):
    dens_file, close = process_file_arg(parameter_density_file)
    try:
        header = parse_header(dens_file, seek = False)
        parameter_indices = get_indices_of_patterns(header, parameter_patterns)
        indices_to_heads = dict(zip(parameter_indices,
                [header[i] for i in parameter_indices]))
        heads_to_dens_tups = dict(zip([header[i] for i in parameter_indices],
                [None for i in range(len(parameter_indices))]))
        if not len(parameter_indices) == len(set(indices_to_heads.itervalues())):
            dens_file.close()
            raise ParameterParsingError('some parameters were found in multiple '
                    'columns in density file {0!r}'.format(dens_file.name))
        for i, line in enumerate(dens_file):
            l = line.strip().split()
            if l:
                for idx in parameter_indices:
                    heads_to_dens_tups[indices_to_heads[idx]] = (float(l[idx]),
                            float(l[idx + 1]))
                yield heads_to_dens_tups
    except:
        raise
    finally:
        if close:
            dens_file.close()

def parse_parameter_density_file(parameter_density_file,
        parameter_patterns = DIV_MODEL_PATTERNS + MODEL_PATTERNS + \
                PSI_PATTERNS + MEAN_TAU_PATTERNS + OMEGA_PATTERNS):
    val_dens_tups = None
    for pd in parameter_density_iter(parameter_density_file):
        if not val_dens_tups:
            val_dens_tups = dict(zip(pd.keys(), [[] for i in range(len(pd))]))
        for k, vd_tup in pd.iteritems():
            val_dens_tups[k].append(vd_tup)
    return val_dens_tups

def parameter_iter(file_obj, include_line = False, include_thetas = False):
    indices = {}
    post_file, close = process_file_arg(file_obj)
    header = parse_header(post_file, seek = False)
    mean_t_indices = get_indices_of_patterns(header, MEAN_TAU_PATTERNS)
    if len(mean_t_indices) > 1:
        post_file.close()
        raise ParameterParsingError('posterior file {0} has {1} mean '
                'tau columns'.format(post_file.name, len(mean_t_indices)))
    if mean_t_indices:
        indices['mean_tau'] = mean_t_indices
    omega_indices = get_indices_of_patterns(header, OMEGA_PATTERNS)
    if len(omega_indices) > 1:
        post_file.close()
        raise ParameterParsingError('posterior file {0} has {1} omega '
                'columns'.format(post_file.name, len(omega_indices)))
    if omega_indices:
        indices['omega'] = omega_indices
    t_indices = get_indices_of_patterns(header, TAU_PATTERNS)
    if t_indices:
        indices['taus'] = t_indices
    if include_thetas:
        a_theta_indices = get_indices_of_patterns(header, A_THETA_PATTERNS)
        d1_theta_indices = get_indices_of_patterns(header, D1_THETA_PATTERNS)
        d2_theta_indices = get_indices_of_patterns(header, D2_THETA_PATTERNS)
        if a_theta_indices:
            indices['a_thetas'] = a_theta_indices
        if d1_theta_indices:
            indices['d1_thetas'] = d1_theta_indices
        if d2_theta_indices:
            indices['d2_thetas'] = d2_theta_indices
    psi_indices = get_indices_of_patterns(header, PSI_PATTERNS)
    if len(psi_indices) > 1:
        post_file.close()
        raise ParameterParsingError('posterior file {0} has {1} psi '
                'columns'.format(post_file.name, len(psi_indices)))
    if psi_indices:
        indices['psi'] = psi_indices
    model_indices = get_indices_of_patterns(header, MODEL_PATTERNS)
    if len(model_indices) > 1:
        post_file.close()
        raise ParameterParsingError('posterior file {0} has {1} model '
                'columns'.format(post_file.name, len(model_indices)))
    if model_indices:
        indices['model'] = model_indices
    div_model_indices = get_indices_of_patterns(header, DIV_MODEL_PATTERNS)
    if len(div_model_indices) > 1:
        post_file.close()
        raise ParameterParsingError('posterior file {0} has {1} div model '
                'columns'.format(post_file.name, len(div_model_indices)))
    if div_model_indices:
        indices['div_model'] = div_model_indices
    samples = dict(zip(indices.keys(), [None for i in range(len(indices))]))
    for i, line in enumerate(post_file):
        l = line.strip().split()
        if l:
            if len(l) != len(header):
                post_file.close()
                raise ParameterParsingError('posterior file {0} has '
                        '{1} columns at line {2}; expecting {3}'.format(
                                post_file.name, len(l), i + 2, len(header)))
            for k, idx_list in indices.iteritems():
                if k in ['mean_tau', 'omega']:
                    samples[k] = [float(l[i]) for i in idx_list]
                elif k in ['psi', 'model', 'div_model']:
                    samples[k] = [int(l[i]) for i in idx_list]
                elif k in ['taus', 'a_thetas', 'd1_thetas', 'd2_thetas']:
                    samples[k] = [[float(l[i]) for i in idx_list]]
                else:
                    post_file.close()
                    raise ParameterParsingError('unexpected key {0!r}; '
                            'posterior file {1}, line {2}'.format(
                                k, post_file.name, i+2))
            if include_line:
                yield samples, l
            else:
                yield samples
    if close:
        post_file.close()

def parse_parameters(file_obj, include_thetas = False):
    samples = None
    for s in parameter_iter(file_obj, include_thetas = include_thetas):
        if not samples:
            samples = dict(zip(s.keys(), [[] for i in range(len(s))]))
        for k, v in s.iteritems():
            samples[k].extend(v)
    return samples

def add_div_model_column(in_file, out_file, div_models_to_indices,
        compresslevel = None):
    header = parse_header(in_file)
    if get_indices_of_patterns(header, DIV_MODEL_PATTERNS) != []:
        raise ParameterParsingError('posterior file {0} already has a '
                'divergence model column'.format(
                getattr(in_file, 'name', in_file)))
    header.insert(0, 'PRI.div.model')
    out, close = process_file_arg(out_file, 'w', compresslevel=compresslevel)
    out.write('{0}\n'.format('\t'.join(header)))
    other_index = max(div_models_to_indices.itervalues()) + 1
    for parameters, line in parameter_iter(in_file, include_line = True):
        if not parameters.has_key('taus'):
            out.close()
            raise ParameterParsingError('posterior file {0} does not contain '
                    'divergence time vector'.format(
                    getattr(file_obj, 'name', file_obj)))
        if MSBAYES_SORT_INDEX.current_value() == 0:
            ip = pymsbayes.utils.stats.Partition(parameters['taus'][0])
        else:
            ip = pymsbayes.utils.stats.IntegerPartition(parameters['taus'][0])
        idx = div_models_to_indices.get(ip.key, other_index)
        line.insert(0, str(idx))
        out.write('{0}\n'.format('\t'.join(line)))
    if close:
        out.close()

def strip_div_model_column(in_file, out_file, compresslevel = None):
    header = parse_header(in_file)
    div_indices = set(get_indices_of_patterns(header, DIV_MODEL_PATTERNS))
    indices = list(set(range(len(header))) - div_indices)
    in_stream, close_in = process_file_arg(in_file, 'rU')
    out_stream, close_out = process_file_arg(out_file, 'w',
            compresslevel=compresslevel)
    for line_num, line in enumerate(in_stream):
        l = line.strip().split()
        out_stream.write('{0}\n'.format('\t'.join([l[i] for i in indices])))
    if close_in:
        in_stream.close()
    if close_out:
        out_stream.close()

def parse_abctoolbox_summary_file(file_obj):
    sum_file, close = process_file_arg(file_obj, 'rU')
    header = sum_file.next().strip().split()
    param_names = header[1:]
    params_to_indices = dict(zip(param_names,
            [i for i in range(len(param_names))]))
    summaries = dict(zip(param_names, [{} for i in range(len(param_names))]))
    for line in sum_file:
        l = line.strip().split()
        stat_name = l.pop(0)
        for k, d in summaries.iteritems():
            d[stat_name] = float(l[params_to_indices[k]])
    if close:
        sum_file.close()
    return summaries

##############################################################################
## ABACUS output parsers

def parse_summary_file(file_obj):
    f, close = process_file_arg(file_obj)
    lines = []
    for l in f:
        l = l.strip()
        if l:
            lines.append(l)
    if close:
        f.close()
    if len(lines) != 4:
        raise SummaryFileParsingError('summary file {0} has {1} lines'.format(
                f.name, len(lines)))
    header = lines[0].split()
    means = [float(x) for x in lines[1].split()]
    std_devs = [float(x) for x in lines[2].split()]
    sample_sizes = [int(x) for x in lines[3].split()]
    if not len(header) == len(means) == len(std_devs) == len(sample_sizes):
        raise SummaryFileParsingError('lines of summary file {0} have unequal '
                'numbers of columns'.format(f.name))
    d = {}
    for i in range(len(header)):
        d[header[i]] = {'mean': means[i],
                        'std_deviation': std_devs[i],
                        'n': sample_sizes[i]}
    return d, header

##############################################################################
## results parsing

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
            self.prior_configs[k] = MsBayesConfig(v)
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
        pattern = prefix + '*-posterior-sample*'
        result_files = [os.path.basename(x) for x in glob.glob(pattern)]
        result_indices = []
        for f in result_files:
            m = self.result_file_name_pattern.match(f)
            result_indices.append(int(m.group('result_index')))
        result_indices.sort()
        return result_indices

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
            if TAU_PATTERNS[0].match(k) or D_THETA_PATTERNS[0].match(k) or \
                    A_THETA_PATTERNS[0].match(k):
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
                ex = get_sublist_greater_than(div_times, tau_max)
                ex_glm = get_sublist_greater_than(div_times, tau_max_glm)
                prior_odds = (prior_prob_of_exclusion /
                        (1 - prior_prob_of_exclusion))
                post_odds = (prob_of_exclusion / (1 - prob_of_exclusion))
                post_odds_glm = (prob_of_exclusion_glm /
                        (1 - prob_of_exclusion_glm))
                if almost_equal(prior_odds, 0.0):
                    bf_of_exclusion = float('inf')
                    bf_of_exclusion_glm = float('inf')
                    if almost_equal(bf_of_exclusion, 0.0):
                        bf_of_exclusion = 0.0
                    if almost_equal(bf_of_exclusion_glm, 0.0):
                        bf_of_exclusion_glm = 0.0
                else:
                    bf_of_exclusion = post_odds / prior_odds
                    bf_of_exclusion_glm = post_odds_glm / prior_odds
                bf_ex = []
                bf_ex_glm = []
                if bf_of_exclusion > 1.0:
                    bf_ex = get_sublist_greater_than(div_times, bf_tau_max)
                if bf_of_exclusion_glm > 1.0:
                    bf_ex_glm = get_sublist_greater_than(div_times, bf_tau_max)
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
        header = parse_header(observed_stream, sep = '\t', strict = True,
                seek = False)
        parameter_indices = get_indices_of_patterns(header, PARAMETER_PATTERNS)
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
                out_path = get_new_path(out_path)
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
        
def parse_omega_results_file(file_obj):
    s_iter = spreadsheet_iter([file_obj], sep = '\t')
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
    s_iter = spreadsheet_iter([file_obj], sep = '\t')
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
    s_iter = spreadsheet_iter([file_obj], sep = '\t')
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


class UnorderedDivergenceModelResults(object):
    def __init__(self, div_model_results_path,
            inclusion_threshold = None):
        self.path = div_model_results_path
        self.inclusion_threshold = inclusion_threshold
        self.models = []
        self.n = 0
        self.cumulative_prob = 0.0
        header = parse_header(self.path)
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
        ss_iter = spreadsheet_iter([file_stream])
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
        post = parse_parameters(self.path)
        if not post.has_key('taus'):
            raise Exception('posterior sample in {0} does not contain a '
                    'divergence time vector'.format(self.posterior_path))
        div_models = pymsbayes.utils.stats.IntegerPartitionCollection(
                post['taus'])
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

class NumberOfDivergencesSummary(object):
    def __init__(self,
            config_path,
            psi_results_path,
            posterior_summary_path = None,
            num_prior_samples = 100000,
            num_processors = 4):
        self.config_path = config_path
        self.config = MsBayesConfig(config_path)
        self.psi_results_path = psi_results_path
        self.posterior_summary_path = posterior_summary_path
        self.num_prior_samples = num_prior_samples
        self.num_processors = num_processors
        self.psi_posterior_probs = {}
        self.psi_prior_probs = None
        self.unordered_model_prior_probs = None
        self.ordered_model_prior_probs = None
        self.psi_bayes_factors = {}
        self.omega = None
        self.omega_hpd = None
        self.plot = None
        self._parse_psi_results_path()
        self._parse_posterior_summary_file()
        self._simulate_psi_prior_probs()
        self._get_unordered_model_prior_probs()
        self._get_ordered_model_prior_probs()
        self._get_psi_bayes_factors()

    def _simulate_psi_prior_probs(self):
        if self.config.div_model_prior == 'dpp':
            prob_team = pymsbayes.teams.ModelProbabilityEstimatorTeam(
                    config_paths = [self.config_path],
                    num_samples = self.num_prior_samples,
                    num_processors = self.num_processors)
            prob_team.start()
            self.psi_prior_probs = prob_team.psi_probs[self.config_path]
        elif config.div_model_prior == 'uniform':
            ips = stats.IntegerPartition.number_of_int_partitions_by_k(
                    num_elements = npairs)
            n = sum(ips)
            for i in range(1, self.config.npairs + 1):
                self.psi_prior_probs[i] = ips[i-1] / float(n)
        elif config.div_model_prior == 'psi':
            self.psi_prior_probs = {}
            for i in range(1, self.config.npairs + 1):
                self.psi_prior_probs[i] = 1.0 / npairs

    def _get_unordered_model_prior_probs(self):
        if not self.psi_prior_probs:
            return
        ips = pymsbayes.stats.IntegerPartition.number_of_int_partitions_by_k(
                num_elements = len(self.psi_prior_probs))
        self.unordered_model_prior_probs = {}
        for k, p in self.psi_prior_probs.iteritems():
            self.unordered_model_prior_probs[k] = p / float(ips[k-1])

    def _get_ordered_model_prior_probs(self):
        if not self.psi_prior_probs:
            return
        part = pymsbayes.stats.Partition([0] * len(self.psi_prior_probs))
        self.ordered_model_prior_probs = {}
        for k, p in self.psi_prior_probs.iteritems():
            self.ordered_model_prior_probs[k] = p / float(
                    part.number_of_partitions_into_k_subsets(k))

    def _get_psi_bayes_factors(self):
        for k, p in self.psi_posterior_probs.iteritems():
            p_not = 1.0 - p
            prior = self.psi_prior_probs[k]
            prior_not = 1.0 - p
            self.psi_bayes_factors[k] = 2 * math.log(
                    (p / p_not) / (prior / prior_not))

    def _parse_psi_results_path(self):
        for d in spreadsheet_iter([self.psi_results_path]):
            self.psi_posterior_probs[int(d['num_of_div_events'])] = float(
                    d['estimated_prob'])

    def _parse_posterior_summary_file(self):
        results = parse_posterior_summary_file(
                self.posterior_summary_path)
        self.omega = float(results['PRI.omega']['median'])
        self.omega_hpd = [float(results['PRI.omega']['HPD_95_interval'][0]),
                float(results['PRI.omega']['HPD_95_interval'][1])]
        if self.omega_hpd[0] < 0.0:
            self.omega_hpd[0] = 0.0

    def create_plot(self,
            plot_label_size = 10.0,
            right_text_size = 10.0,
            x_label_size = 10.0,
            y_label_size = 10.0,
            xtick_label_size = 10.0,
            ytick_label_size = 8.0,
            height = 6.0,
            width = 3.0,
            margin_bottom = 0.0,
            margin_left = 0.0,
            margin_top = 0.98,
            margin_right = 1.0,
            padding_between_vertical = 1.0):
        if not pymsbayes.plotting.MATPLOTLIB_AVAILABLE:
            _LOG.warning('matplotlib is not available; cannot create plot')
            return
        right_text = ''
        if not self.omega is None:
            right_text = r'$D_T = {0:.2f} ({1:.2f}-{2:.2f})$'.format(self.omega,
                    self.omega_hpd[0],
                    self.omega_hpd[1])
        keys = sorted(self.psi_posterior_probs.keys())
        psi_posterior_bar_data = pymsbayes.plotting.BarData(
                    values = [self.psi_posterior_probs[k] for k in keys],
                    labels = keys,
                    width = 1.0,
                    orientation = 'vertical',
                    color = '0.5',
                    edgecolor = '0.5',
                    label_size = xtick_label_size,
                    measure_tick_label_size = ytick_label_size,
                    zorder = 0)
        psi_prior_bar_data = pymsbayes.plotting.BarData(
                    values = [self.psi_prior_probs[k] for k in keys],
                    labels = keys,
                    width = 1.0,
                    orientation = 'vertical',
                    color = '0.5',
                    edgecolor = '0.5',
                    label_size = xtick_label_size,
                    measure_tick_label_size = ytick_label_size,
                    zorder = 0)
        psi_bayes_factor_bar_data = pymsbayes.plotting.BarData(
                    values = [self.psi_bayes_factors[k] for k in keys],
                    labels = keys,
                    width = 1.0,
                    orientation = 'vertical',
                    color = '0.5',
                    edgecolor = '0.5',
                    label_size = xtick_label_size,
                    measure_tick_label_size = ytick_label_size,
                    zorder = 0)
        plots = [
                pymsbayes.plotting.ScatterPlot(
                        bar_data_list = [psi_posterior_bar_data],
                        right_text = right_text,
                        y_label = 'Posterior probability',
                        y_label_size = y_label_size,
                        ),
                pymsbayes.plotting.ScatterPlot(
                        bar_data_list = [psi_prior_bar_data],
                        y_label = 'Prior probability',
                        y_label_size = y_label_size,
                        ),
                pymsbayes.plotting.ScatterPlot(
                        bar_data_list = [psi_bayes_factor_bar_data],
                        y_label = '2ln(Bayes factor)',
                        y_label_size = y_label_size,
                        x_label = r'Number of divergence events, $|\tau|$',
                        x_label_size = x_label_size,
                        ),
                ]
        for p in plots:
            p.right_text_size = right_text_size
            p.plot_label_size = plot_label_size
            yticks = [i for i in p.ax.get_yticks()]
            ytick_labels = [i for i in yticks]
            if len(ytick_labels) > 5:
                for i in range(1, len(ytick_labels), 2):
                    ytick_labels[i] = ''
            yticks_obj = pymsbayes.plotting.Ticks(ticks = yticks,
                    labels = ytick_labels,
                    size = ytick_label_size)
            p.yticks_obj = yticks_obj

        self.plot = pymsbayes.plotting.PlotGrid(subplots = plots,
                num_columns = 1,
                share_x = True,
                share_y = False,
                height = height,
                width = width,
                auto_height = False)
        self.plot.auto_adjust_margins = False
        self.plot.margin_top = margin_top
        self.plot.margin_bottom = margin_bottom
        self.plot.margin_right = margin_right
        self.plot.margin_left = margin_left
        self.plot.padding_between_vertical = padding_between_vertical
        self.plot.reset_figure()
        self.plot.reset_figure()

    def save_plot(self, path):
        if self.plot:
            self.plot.savefig(path)

