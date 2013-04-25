#! /usr/bin/env python

import os
import sys
import re
import subprocess
import time
import shutil
import traceback
from cStringIO import StringIO

from pymsbayes.fileio import expand_path, process_file_arg, FileStream, open
from pymsbayes.utils.tempfs import TempFileSystem
from pymsbayes.utils import BIN_DIR
from pymsbayes.utils.functions import (get_random_int, get_indices_of_patterns,
        reduce_columns, is_dir)
from pymsbayes.utils.errors import WorkerExecutionError, PriorMergeError
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)


##############################################################################
## msBayes prior header patterns

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

##############################################################################
## other globals

VALID_REJECTION_TOOLS = ['msreject', 'abctoolbox']
VALID_REGRESSION_METHODS = ['llr', 'glm']

##############################################################################
## functions for manipulating prior files

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


##############################################################################
## Base class for all workers

class Worker(object):
    total = 0
    def __init__(self,
            stdout_path = None,
            stderr_path = None,
            tag = ''):
        self.__class__.total += 1
        self.stdout_path = stdout_path
        self.stderr_path = stderr_path
        self.cmd = []
        self.process = None
        self.finished = False
        self.exit_code = None
        self.stdout = None
        self.stderr = None
        self.tag = tag

    def get_stderr(self):
        if not self.stderr_path:
            return self.stderr
        try:
            return open(self.stderr_path, 'rU').read()
        except IOError, e:
            _LOG.error('Could not open stderr file')
            raise e

    def start(self):
        try:
            self._pre_process()
        except:
            e = StringIO()
            traceback.print_exc(file=e)
            _LOG.error('Error during pre-processing:\n{0}'.format(
                    e.getvalue()))
            raise
        _LOG.info('Starting process with following command:\n\t'
                '{0}'.format(' '.join(self.cmd)))
        if self.stdout_path:
            sout = open(self.stdout_path, 'a')
        else:
            sout = subprocess.PIPE
        if self.stderr_path:
            serr = open(self.stderr_path, 'w')
        else:
            serr = subprocess.PIPE
        self.process = subprocess.Popen(self.cmd,
                stdout = sout,
                stderr = serr,
                shell = False)
        self.stdout, self.stderr = self.process.communicate()
        self.exit_code = self.process.wait()
        if self.stdout_path:
            sout.close()
        if self.stderr_path:
            serr.close()
        if self.exit_code != 0:
            if self.stderr_path:
                self.stderr = open(self.stderr_path, 'rU').read()
            _LOG.error('execution failed')
            raise WorkerExecutionError('{0} failed. stderr:\n{1}'.format(
                self.name, self.stderr))
        try:
            self._post_process()
        except:
            e = StringIO()
            traceback.print_exc(file=e)
            _LOG.error('Error during post-processing:\n{0}'.format(
                    e.getvalue()))
            raise
        self.finished = True

    def _post_process(self):
        pass

    def _pre_process(self):
        pass


##############################################################################
## functions for managing msbayes workers

def merge_priors(workers, prior_path, header_path=None, include_header=False):
    out, close = process_file_arg(prior_path, 'w')
    h = None
    std = None
    for w in workers:
        if not h:
            h = w.header
            std = w
            if include_header:
                out.write('{0}\n'.format('\t'.join(h)))
        if w.header != h:
            raise PriorMergeError('Workers {0} and {1} have different prior '
                    'headers. Cannot merge!'.format(std.name, w.name))
        prior_file = open(w.prior_path, 'rU')
        for line in prior_file:
            if not HEADER_PATTERN.match(line.strip()):
                out.write(line)
        prior_file.close()
    if close:
        out.close()
    if not header_path:
        header_path = prior_path + '.header'
    out, close = process_file_arg(header_path, 'w')
    out.write('{0}\n'.format('\t'.join(h)))
    if close:
        out.close()
    return prior_path, header_path

def merge_prior_files(paths, dest_path):
    out, close = process_file_arg(dest_path, 'w')
    h = None
    ncols = None
    for p in paths:
        header = None
        f, f_close = process_file_arg(p, 'rU')
        for line_num, line in enumerate(f):
            if line_num == 0:
                if HEADER_PATTERN.match(line.strip()):
                    header = line.strip().split()
                    if h and h != header:
                        raise PriorMergeError('prior files {0} and {1} have '
                                'different headers'.format(f.name, std.name))
                if not ncols:
                    ncols = len(line.strip().split())
                    std = f
                    if header:
                        h = header
                        out.write(line)
                        continue
                if len(line.strip().split()) != ncols:
                    raise PriorMergeError('prior files {0} and {1} do not '
                            'have the same number of columns. Cannot '
                            'merge!'.format(f.name, std.name))
                if not header:
                    out.write(line)
            else:
                out.write(line)
        if f_close:
            f.close()
    if close:
        out.close()

##############################################################################
## msBayes class for generating prior files

class MsBayesWorker(Worker):
    count = 0
    valid_schemas = ['msreject', 'abctoolbox']

    def __init__(self,
            temp_fs,
            sample_size,
            config_path,
            exe_path = None,
            model_index = None,
            sort_index = None,
            report_parameters = True,
            seed = None,
            schema = 'msreject',
            include_header = False,
            stat_patterns=DEFAULT_STAT_PATTERNS,
            parameter_patterns=PARAMETER_PATTERNS,
            stdout_path = None,
            stderr_path = None,
            staging_dir = None,
            write_stats_file = False,
            tag = ''):
        Worker.__init__(self,
                stdout_path = stdout_path,
                stderr_path = stderr_path,
                tag = tag)
        self.__class__.count += 1
        self.name = 'MsBayesWorker-' + str(self.count)
        self.temp_fs = temp_fs
        self.sample_size = int(sample_size)
        self.config_path = expand_path(config_path)
        self.output_dir = self.temp_fs.create_subdir(prefix = self.name + '-')
        if not exe_path:
            exe_path = os.path.join(BIN_DIR, 'msbayes.pl')
        self.exe_path = expand_path(exe_path)
        self.model_index = None
        if model_index != None:
            self.model_index = int(model_index)
        self.sort_index = None
        if sort_index != None:
            self.sort_index = int(sort_index)
        self.report_parameters = report_parameters
        if seed is None:
            self.seed = get_random_int()
        else:
            self.seed = int(seed)
        self.prior_path = self.temp_fs.get_file_path(
                parent = self.output_dir,
                prefix = 'prior-{0}-{1}-'.format(
                        self.sample_size,
                        self.seed),
                create = False)
        self.header_path = self.temp_fs.get_file_path(
                parent = self.output_dir,
                prefix = 'prior-{0}-{1}-header-'.format(
                        self.sample_size,
                        self.seed),
                create = False)
        if not schema.lower() in self.valid_schemas:
            raise ValueError(
                    'schema {0} is not valid. Options are: {1}'.format(
                        schema, ','.join(self.valid_schemas)))
        self.schema = schema.lower()
        self.include_header = include_header
        self.prior_stats_path = None
        if write_stats_file:
            self.prior_stats_path = self.prior_path + '.stats.txt'
            self.temp_fs._register_file(self.prior_stats_path)
        self.stat_patterns = stat_patterns
        self.parameter_patterns = parameter_patterns
        self.header = None
        self.parameter_indices = None
        self.stat_indices = None
        self.staging_dir = None
        self.staging_prior_path = None
        if is_dir(staging_dir):
            self.staging_dir = staging_dir
            self.staging_prior_path = os.path.join(staging_dir,
                    os.path.basename(self.prior_path))
            self.staging_prior_stats_path = self.staging_prior_path
            if write_stats_file:
                self.staging_prior_stats_path = self.staging_prior_path + \
                        '.stats.txt'
        self._update_cmd()

    def _update_cmd(self):
        cmd = [self.exe_path,
               '-r', str(self.sample_size),
               '-S', str(self.seed),
               '-c', self.config_path,]
        if self.staging_dir:
            cmd.extend(['-o', self.staging_prior_path])
        else:
            cmd.extend(['-o', self.prior_path])
        if self.sort_index != None:
            cmd.extend(['-s', str(self.sort_index)])
        if self.model_index != None:
            cmd.extend(['-m', str(self.model_index)])
        if self.report_parameters:
            cmd.append('-p')
        self.cmd = cmd

    def _post_process(self):
        prior_path = self.prior_path
        prior_stats_path = self.prior_stats_path
        if self.staging_dir:
            prior_path = self.staging_prior_path
            prior_stats_path = self.staging_prior_stats_path
        raw_prior_path = prior_path + '.raw'
        shutil.move(prior_path, raw_prior_path)
        if self.prior_stats_path:
            h = observed_stats_for_abctoolbox(
                    in_file = raw_prior_path,
                    out_file = prior_stats_path,
                    stat_patterns = self.stat_patterns)
        if self.schema == 'msreject':
            header = prior_for_msreject(
                    in_file = raw_prior_path,
                    out_file = prior_path,
                    stat_patterns = self.stat_patterns,
                    parameter_patterns = self.parameter_patterns,
                    dummy_patterns = DUMMY_PATTERNS,
                    include_header = self.include_header)
        elif self.schema == 'abctoolbox':
            header = prior_for_abctoolbox(
                    in_file = raw_prior_path,
                    out_file = prior_path,
                    stat_patterns = self.stat_patterns,
                    parameter_patterns = self.parameter_patterns)
        else:
            raise ValueError(
                    'schema {0} is not valid. Options are: {1}'.format(
                        self.schema, ','.join(self.valid_schemas)))
        os.remove(raw_prior_path)
        if header:
            out = open(self.header_path, 'w')
            out.write('{0}\n'.format('\t'.join(header)))
            out.close()
            self._set_header(header)
        if self.staging_dir:
            shutil.move(self.staging_prior_path, self.prior_path)
            if self.prior_stats_path:
                shutil.move(self.staging_prior_stats_path,
                        self.prior_stats_path)

    def _set_header(self, header):
        self.header = header
        self.parameter_indices = get_parameter_indices(
                header_list = self.header,
                parameter_patterns = PARAMETER_PATTERNS)
        self.stat_indices = get_stat_indices(
                header_list = self.header,
                stat_patterns = ALL_STAT_PATTERNS)

        
##############################################################################
## functions for managing rejection workers

def assemble_rejection_workers(temp_fs,
        observed_sims_file,
        prior_path,
        results_dir,
        tolerance = None,
        num_posterior_samples = None,
        num_prior_samples = None,
        posterior_prefix='unadjusted-posterior',
        stat_indices = None,
        continuous_parameter_indices = None,
        discrete_parameter_indices = None,
        reject_path = None,
        regress_path = None,
        regress = True,
        rejection_tool = 'abctoolbox',
        regression_method = 'glm',
        keep_temps = False,
        bandwidth = None,
        num_posterior_quantiles = None):
    if tolerance is None and num_posterior_samples is None:
        raise ValueError('must specify either `tolerance` or '
                '`num_posterior_samples`')
    if rejection_tool.lower() not in VALID_REJECTION_TOOLS:
        raise ValueError('invalid rejection tool {0}; valid options are: '
                '{1}'.format(rejection_tool, ', '.join(VALID_REJECTION_TOOLS)))
    if regression_method.lower() not in VALID_REGRESSION_METHODS:
        raise ValueError('invalid regression method {0}; valid options are: '
                '{1}'.format(regression_method, ', '.join(
                        VALID_REGRESSION_METHODS)))
    if (tolerance is None or num_posterior_samples is None) and \
            num_prior_samples is None:
        num_prior_samples = line_count(prior_path, ignore_headers=True)
    if tolerance is None:
        tolerance = num_posterior_samples / float(num_prior_samples)
    if num_posterior_samples is None:
        num_posterior_samples = int(tolerance * num_prior_samples)
    if num_prior_samples is None:
        num_prior_samples = int(num_posterior_samples / tolerance)
    obs_file, close = process_file_arg(observed_sims_file)
    header = parse_header(obs_file)
    if rejection_tool == 'abctoolbox':
        header = parse_header(prior_path)
    all_stat_indices = get_stat_indices(header,
            stat_patterns=ALL_STAT_PATTERNS)
    parameter_indices = []
    if continuous_parameter_indices:
        parameter_indices += continuous_parameter_indices
    if discrete_parameter_indices:
        parameter_indices += discrete_parameter_indices
    obs_temp_dir = temp_fs.create_subdir(prefix = 'observed-files-')
    header_line = obs_file.next()
    reject_workers = []
    for i, line in enumerate(obs_file):
        rej_obs_path = temp_fs.get_file_path(parent = obs_temp_dir,
                prefix = 'observed-{0}-reject-'.format(i+1),
                create = False)
        reg_obs_path = temp_fs.get_file_path(parent = obs_temp_dir,
                prefix = 'observed-{0}-regress-'.format(i+1),
                create = False)
        l = line.strip().split()
        if rejection_tool.lower() == 'msreject':
            with open(rej_obs_path, 'w') as out:
                out.write(line)
        elif rejection_tool.lower() == 'abctoolbox':
            with open(rej_obs_path, 'w') as out:
                out.write('{0}\n{1}\n'.format(
                        '\t'.join([header[idx] for idx in all_stat_indices]),
                        '\t'.join([str(l[idx]) for idx in all_stat_indices])))
        else:
            raise ValueError('Unexpected rejection tool {0}'.format(
                    rejection_tool))
        if regression_method.lower() == 'llr':
            with open(reg_obs_path, 'w') as out:
                out.write(header_line)
                out.write(line)
        elif regression_method.lower() == 'glm':
            if rejection_tool.lower() == 'abctoolbox':
                reg_obs_path = rej_obs_path
            else:
                with open(reg_obs_path, 'w') as out:
                    out.write('{0}\n{1}\n'.format(
                            '\t'.join([header[idx] for idx in all_stat_indices]),
                            '\t'.join([str(l[idx]) for idx in all_stat_indices])))
        else:
            raise ValueError('Unexpected regression method {0}'.format(
                    regression_method))
        posterior_file_name = '{0}.{1}.txt'.format(posterior_prefix, i+1)
        posterior_path = os.path.join(results_dir, posterior_file_name)
        regression_worker = None
        if regress and regression_method.lower() == 'llr':
            regression_worker = RegressionWorker(
                    observed_path = reg_obs_path,
                    posterior_path = posterior_path,
                    tolerance = 1.0,
                    stat_indices = stat_indices,
                    continuous_parameter_indices = continuous_parameter_indices,
                    discrete_parameter_indices = discrete_parameter_indices,
                    exe_path = regress_path)
        if regress and regression_method.lower() == 'glm':
            regression_worker = ABCToolBoxRegressWorker(
                    temp_fs = temp_fs,
                    observed_path = reg_obs_path,
                    posterior_path = posterior_path,
                    parameter_indices = sorted(parameter_indices),
                    exe_path = regress_path,
                    keep_temps = keep_temps,
                    num_posterior_samples = num_posterior_samples,
                    bandwidth = bandwidth,
                    num_posterior_quantiles = num_posterior_quantiles)
        if rejection_tool.lower() == 'msreject':
            reject_workers.append(MsRejectWorker(
                    header = header,
                    observed_path = rej_obs_path,
                    prior_path = prior_path,
                    tolerance = tolerance,
                    posterior_path = posterior_path,
                    stat_indices = stat_indices,
                    regression_worker = regression_worker,
                    exe_path = reject_path))
        if rejection_tool.lower() == 'abctoolbox':
            reject_workers.append(ABCToolBoxRejectWorker(
                    temp_fs = temp_fs,
                    observed_path = rej_obs_path,
                    prior_path = prior_path,
                    num_posterior_samples = num_posterior_samples,
                    posterior_path = posterior_path,
                    regression_worker = regression_worker,
                    exe_path = reject_path,
                    keep_temps = keep_temps,
                    max_read_sims = int(num_prior_samples + \
                            (0.01 * num_prior_samples))))
    if close:
        obs_file.close()
    return reject_workers


##############################################################################
## worker class for rejection via msreject

class MsRejectWorker(Worker):
    count = 0
    def __init__(self,
            header,
            observed_path,
            prior_path,
            tolerance,
            posterior_path = None,
            stat_indices = None,
            regression_worker = None,
            exe_path = None,
            stderr_path = None,
            tag = ''):
        Worker.__init__(self,
                stdout_path = None,
                stderr_path = stderr_path,
                tag = tag)
        self.__class__.count += 1
        self.name = 'MsRejectWorker-' + str(self.count)
        if not exe_path:
            exe_path = os.path.join(BIN_DIR, 'msReject')
        self.exe_path = expand_path(exe_path)
        self.observed_path = expand_path(observed_path)
        self.prior_path = expand_path(prior_path)
        if not posterior_path:
            posterior_path = observed_path + '.posterior'
        self.posterior_path = expand_path(posterior_path)
        self.header = header
        potential_stat_indices = get_stat_indices(
                self.header,
                stat_patterns = ALL_STAT_PATTERNS)
        if not stat_indices:
            self.stat_indices = potential_stat_indices
        else:
            diff = set(stat_indices) - set(potential_stat_indices)
            if len(diff) > 0:
                raise ValueError('stat indices are not valid')
            self.stat_indices = stat_indices
        self.stdout_path = self.posterior_path
        self.tolerance = float(tolerance)
        self.regression_worker = regression_worker
        self._update_cmd()

    def _update_cmd(self):
        cmd = [self.exe_path,
               self.observed_path,
               self.prior_path,
               str(self.tolerance)]
        cmd.extend([str(i+1) for i in self.stat_indices])
        self.cmd = cmd
    
    def _post_process(self):
        if self.regression_worker:
            self.regression_worker.start()

    def _pre_process(self):
        stdout = open(self.stdout_path, 'w')
        stdout.write('{0}\n'.format('\t'.join(self.header)))
        stdout.close()
        

##############################################################################
## worker class for regression via regress_cli.r

class RegressionWorker(Worker):
    count = 0
    def __init__(self,
            observed_path,
            posterior_path,
            summary_path = None,
            adjusted_path = None,
            tolerance = 1.0,
            stat_indices = None,
            continuous_parameter_indices = None,
            discrete_parameter_indices = None,
            exe_path = None,
            stdout_path = None,
            stderr_path = None,
            tag = ''):
        Worker.__init__(self,
                stdout_path = stdout_path,
                stderr_path = stderr_path,
                tag = tag)
        self.__class__.count += 1
        self.name = 'RegressionWorker-' + str(self.count)
        if not exe_path:
            exe_path = os.path.join(BIN_DIR, 'regress_cli.r')
        self.exe_path = expand_path(exe_path)
        self.observed_path = expand_path(observed_path)
        self.posterior_path = expand_path(posterior_path)
        self.tolerance = float(tolerance)
        self.header = None
        self.stat_indices = stat_indices
        self.continuous_parameter_indices = continuous_parameter_indices
        self.discrete_parameter_indices = discrete_parameter_indices
        if not summary_path:
            summary_path = self.posterior_path + '.regression-summary.txt'
        self.summary_path = summary_path
        if not adjusted_path:
            adjusted_path = self.posterior_path + '.regression-adjusted.txt'
        self.adjusted_path = adjusted_path

    def _pre_process(self):
        self.header = parse_header(self.posterior_path)
        potential_stat_indices = get_stat_indices(
                self.header,
                stat_patterns = ALL_STAT_PATTERNS)
        if not self.stat_indices:
            self.stat_indices = potential_stat_indices
        else:
            diff = set(self.stat_indices) - set(potential_stat_indices)
            if len(diff) > 0:
                raise ValueError('stat indices are not valid')
        if not self.continuous_parameter_indices:
            self.continuous_parameter_indices = get_parameter_indices(
                    self.header,
                    parameter_patterns=(MEAN_TAU_PATTERNS+OMEGA_PATTERNS))
        if not self.discrete_parameter_indices:
            self.discrete_parameter_indices = get_parameter_indices(
                    self.header,
                    parameter_patterns=(MODEL_PATTERNS+PSI_PATTERNS))
        self._update_cmd()

    def _update_cmd(self):
        cmd = [self.exe_path,
               '-t', str(self.tolerance),
               '--observed-path={0}'.format(self.observed_path),
               '--posterior-path={0}'.format(self.posterior_path),
               '--summary-path={0}'.format(self.summary_path),
               '--adjusted-path={0}'.format(self.adjusted_path),
               '-c', ','.join(
                    [str(i+1) for i in self.continuous_parameter_indices]),
               '-d', ','.join(
                    [str(i+1) for i in self.discrete_parameter_indices]),
               '-s', ','.join(
                    [str(i+1) for i in self.stat_indices]),]
        self.cmd = cmd

class ABCToolBoxRejectWorker(Worker):
    count = 0
    def __init__(self,
            temp_fs,
            observed_path,
            prior_path,
            num_posterior_samples,
            posterior_path = None,
            regression_worker = None,
            exe_path = None,
            stdout_path = None,
            stderr_path = None,
            keep_temps = False,
            max_read_sims = 10000000,
            tag = ''):
        Worker.__init__(self,
                stdout_path = stdout_path,
                stderr_path = stderr_path,
                tag = tag)
        self.__class__.count += 1
        self.name = 'ABCToolBoxRejectWorker-' + str(self.count)
        self.temp_fs = temp_fs
        if not exe_path:
            exe_path = os.path.join(BIN_DIR, 'ABCestimator')
        self.exe_path = expand_path(exe_path)
        self.output_dir = self.temp_fs.create_subdir(prefix = self.name + '-')
        self.output_prefix = os.path.join(self.output_dir, 
                self.temp_fs.token_id + '_ABC_GLM_')
        self.cfg_path = self.temp_fs.get_file_path(parent = self.output_dir,
                prefix = 'cfg-',
                create = False)
        self.observed_path = expand_path(observed_path)
        self.prior_path = expand_path(prior_path)
        if not posterior_path:
            posterior_path = observed_path + '.posterior'
        self.posterior_path = expand_path(posterior_path)
        self.header = None
        self.stats_header = None
        self.parameter_indices = None
        self.num_posterior_samples = int(num_posterior_samples)
        self.regression_worker = regression_worker
        self.max_read_sims = int(max_read_sims)

    def _pre_process(self):
        self.header = parse_header(self.prior_path)
        self.parameter_indices = sorted(
                get_parameter_indices(self.header,
                        parameter_patterns=(PARAMETER_PATTERNS)) + \
                get_dummy_indices(self.header))
        self.stats_header = parse_header(self.observed_path)
        cfg = self._compose_cfg_string()
        out = open(self.cfg_path, 'w')
        out.write(cfg.getvalue())
        out.close()
        self._update_cmd()

    def _post_process(self):
        post_path = self.output_prefix + 'BestSimsParamStats_Obs0.txt'
        if not os.path.exists(post_path):
            raise Exception('{0} did not produce a posterior file.\n'
                    'Here is the stderr from the subprocess:\n{1}'.format(
                            self.name, self.get_stderr()))
        # remove extra two columns added by abctoolbox
        with open(self.posterior_path, 'w') as o:
            with open(post_path, 'rU') as i:
                for l in i:
                    o.write('{0}\n'.format('\t'.join(l.strip().split()[2:])))
        self.temp_fs.remove_dir(self.output_dir)
        if self.regression_worker:
            self.regression_worker.start()

    def _update_cmd(self):
        cmd = [self.exe_path, self.cfg_path]
        self.cmd = cmd
        
    def _compose_cfg_string(self):
        cfg = StringIO()
        cfg.write('estimationType standard\n')
        cfg.write('simName {0}\n'.format(self.prior_path))
        cfg.write('obsName {0}\n'.format(self.observed_path))
        cfg.write('params {0}\n'.format(','.join(
                [str(i+1) for i in self.parameter_indices])))
        cfg.write('numRetained {0}\n'.format(self.num_posterior_samples))
        cfg.write('diracPeakWidth 0\n') # No GLM regression
        cfg.write('posteriorDensityPoints 0\n')
        cfg.write('stadardizeStats 1\n')
        cfg.write('writeRetained 1\n')
        cfg.write('maxReadSims {0}\n'.format(self.max_read_sims))
        cfg.write('outputPrefix {0}\n'.format(self.output_prefix))
        return cfg

class ABCToolBoxRegressWorker(Worker):
    count = 0
    def __init__(self,
            temp_fs,
            observed_path,
            posterior_path,
            parameter_indices = None,
            summary_path = None,
            adjusted_path = None,
            exe_path = None,
            stdout_path = None,
            stderr_path = None,
            keep_temps = False,
            bandwidth = None,
            num_posterior_samples = None,
            num_posterior_quantiles = None,
            tag = '',
            ):
        Worker.__init__(self,
                stdout_path = stdout_path,
                stderr_path = stderr_path,
                tag = tag)
        self.__class__.count += 1
        self.name = 'ABCToolBoxRegressWorker-' + str(self.count)
        self.temp_fs = temp_fs
        if not exe_path:
            exe_path = os.path.join(BIN_DIR, 'ABCestimator')
        self.exe_path = expand_path(exe_path)
        self.output_dir = self.temp_fs.create_subdir(prefix = self.name + '-')
        self.output_prefix = os.path.join(self.output_dir, 
                self.temp_fs.token_id + '_ABC_GLM_')
        self.cfg_path = self.temp_fs.get_file_path(prefix = 'cfg-',
                create = False)
        self.observed_path = expand_path(observed_path)
        self.posterior_path = expand_path(posterior_path)
        self.parameter_indices = parameter_indices
        self.header = None
        self.stats_header = None
        if not num_posterior_quantiles:
            num_posterior_quantiles = 1000
        self.num_posterior_quantiles = int(num_posterior_quantiles)
        self.num_posterior_samples = num_posterior_samples
        self.bandwidth = bandwidth
        if not summary_path:
            summary_path = self.posterior_path + '.regression-summary.txt'
        self.summary_path = summary_path
        if not adjusted_path:
            adjusted_path = self.posterior_path + '.regression-adjusted.txt'
        self.adjusted_path = adjusted_path

    def _pre_process(self):
        self.header = parse_header(self.posterior_path)
        if not self.num_posterior_samples:
            self.num_posterior_samples = line_count(self.posterior_path,
                    ignore_headers=True)
        if not self.bandwidth:
            self.bandwidth = 2 / float(self.num_posterior_samples)
        potential_stat_indices = get_stat_indices(
                self.header,
                stat_patterns = ALL_STAT_PATTERNS)
        if not self.parameter_indices:
            self.parameter_indices = sorted(get_parameter_indices(
                    self.header,
                    parameter_patterns=(MEAN_TAU_PATTERNS + \
                            OMEGA_PATTERNS + \
                            MODEL_PATTERNS + \
                            PSI_PATTERNS)))
        if len(set.intersection(set(potential_stat_indices),
                set(self.parameter_indices))) > 0:
            raise ValueError('parameter indices are not valid. '
                    'they contain stat indices!')
        self.stats_header = parse_header(self.observed_path)
        cfg = self._compose_cfg_string()
        out = open(self.cfg_path, 'w')
        out.write(cfg.getvalue())
        out.close()
        self._update_cmd()

    def _post_process(self):
        summary_path = self.output_prefix + 'PosteriorCharacteristics_Obs0.txt'
        adjusted_path = self.output_prefix + 'PosteriorEstimates_Obs0.txt'
        shutil.move(summary_path, self.summary_path)
        shutil.move(adjusted_path, self.adjusted_path)
        self.temp_fs.remove_dir(self.output_dir)

    def _update_cmd(self):
        cmd = [self.exe_path, self.cfg_path]
        self.cmd = cmd
        
    def _compose_cfg_string(self):
        cfg = StringIO()
        cfg.write('estimationType standard\n')
        cfg.write('simName {0}\n'.format(self.posterior_path))
        cfg.write('obsName {0}\n'.format(self.observed_path))
        cfg.write('params {0}\n'.format(','.join(
                [str(i+1) for i in self.parameter_indices])))
        cfg.write('diracPeakWidth {0}\n'.format(self.bandwidth))
        cfg.write('posteriorDensityPoints {0}\n'.format(
                self.num_posterior_quantiles))
        cfg.write('stadardizeStats 1\n')
        cfg.write('writeRetained 0\n')
        cfg.write('maxReadSims {0}\n'.format(self.num_posterior_samples + 100))
        cfg.write('outputPrefix {0}\n'.format(self.output_prefix))
        return cfg

