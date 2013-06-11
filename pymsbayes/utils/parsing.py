#! /usr/bin/env python

import sys
import os
import re

from pymsbayes.fileio import process_file_arg
from pymsbayes.utils.errors import (SummaryFileParsingError,
        ParameterParsingError)
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

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
        re.compile(r'\s*PRI\.(?!numTauClass)\S+\s*$'),
        ]
DEFAULT_STAT_PATTERNS = [
        re.compile(r'\s*pi\.\d+\s*'),
        re.compile(r'\s*wattTheta\.\d+\s*'),
        re.compile(r'\s*pi\.net\.\d+\s*'),
        re.compile(r'\s*tajD\.denom\.\d+\s*'),
        ]
ALL_STAT_PATTERNS = [
        re.compile(r'\s*(?!PRI)\S+\s*$'),
        ]
DUMMY_PATTERNS = [
        re.compile(r'\s*PRI\.numTauClass\s*')
        ]
MODEL_PATTERNS = [
        re.compile(r'\s*PRI\.model\s*'),
        ]
TAU_PATTERNS = [
        re.compile(r'\s*PRI\.t\.\d+\s*'),
        ]
D_THETA_PATTERNS = [
        re.compile(r'\s*PRI\.d[12]Theta\.\d+\s*'),
        ]
A_THETA_PATTERNS = [
        re.compile(r'\s*PRI\.aTheta\.\d+\s*'),
        ]
PSI_PATTERNS = [
        re.compile(r'\s*PRI\.Psi\s*'),
        ]
MEAN_TAU_PATTERNS = [
        re.compile(r'\s*PRI\.E\.t\s*'),
        ]
OMEGA_PATTERNS = [
        re.compile(r'\s*PRI\.omega\s*'),
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
        pattern_str = r'\s*{0}.*\s*'.format(prefix.replace('.', '\.'))
        if ignore_case:
            patterns.append(re.compile(pattern_str, re.IGNORECASE))
        else:
            patterns.append(re.compile(pattern_str))
    return patterns

def parse_header(file_obj, sep='\t'):
    file_stream, close = process_file_arg(file_obj, 'rU')
    header_line = file_stream.next()
    if not HEADER_PATTERN.match(header_line):
        raise Exception('did not find header in {0}'.format(file_stream.name))
    header = header_line.strip().split(sep)
    if close:
        file_stream.close()
    else:
        file_stream.seek(0)
    return header

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

def parse_parameters(file_obj):
    samples = {}
    indices = {}
    post_file, close = process_file_arg(file_obj)
    header = parse_header(post_file)
    mean_t_indices = get_indices_of_patterns(header, MEAN_TAU_PATTERNS)
    if len(mean_t_indices) > 1:
        raise ParameterParsingError('posterior file {0} has {1} mean '
                'tau columns'.format(post_file.name, len(mean_t_indices)))
    if mean_t_indices:
        samples['mean_tau'] = []
        indices['mean_tau'] = mean_t_indices
    omega_indices = get_indices_of_patterns(header, OMEGA_PATTERNS)
    if len(omega_indices) > 1:
        raise ParameterParsingError('posterior file {0} has {1} omega '
                'columns'.format(post_file.name, len(omega_indices)))
    if omega_indices:
        samples['omega'] = []
        indices['omega'] = omega_indices
    t_indices = get_indices_of_patterns(header, TAU_PATTERNS)
    if t_indices:
        samples['taus'] = []
        indices['taus'] = t_indices
    psi_indices = get_indices_of_patterns(header, PSI_PATTERNS)
    if len(psi_indices) > 1:
        raise ParameterParsingError('posterior file {0} has {1} psi '
                'columns'.format(post_file.name, len(psi_indices)))
    if psi_indices:
        samples['psi'] = []
        indices['psi'] = psi_indices
    model_indices = get_indices_of_patterns(header, MODEL_PATTERNS)
    if len(model_indices) > 1:
        raise ParameterParsingError('posterior file {0} has {1} model '
                'columns'.format(post_file.name, len(model_indices)))
    if model_indices:
        samples['model'] = []
        indices['model'] = model_indices
    post_file.next() # header
    for i, line in enumerate(post_file):
        l = line.strip().split()
        if l:
            if len(l) != len(header):
                raise ParameterParsingError('posterior file {0} has '
                        '{1} columns at line {2}; expecting {3}'.format(
                                post_file.name, len(l), i + 2, len(header)))
            for k, idx_list in indices.iteritems():
                if k in ['mean_tau', 'omega']:
                    samples[k].extend([float(l[i]) for i in idx_list])
                elif k in ['psi', 'model']:
                    samples[k].extend([int(l[i]) for i in idx_list])
                elif k == 'taus':
                    samples[k].append([float(l[i]) for i in idx_list])
                else:
                    raise ParameterParsingError('unexpected key {0!r}; '
                            'posterior file {1}, line {2}'.format(
                                k, post_file.name, i+2))
    if close:
        post_file.close()
    return samples

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

