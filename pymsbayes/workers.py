#! /usr/bin/env python

import os
import sys
import re
import subprocess
import time
import shutil
import traceback
from cStringIO import StringIO
from configobj import ConfigObj

from pymsbayes.fileio import (expand_path, process_file_arg, FileStream, open,
        GzipFileStream)
from pymsbayes.utils.tempfs import TempFileSystem
from pymsbayes.utils import get_tool_path
from pymsbayes.utils.functions import (get_random_int, get_indices_of_patterns,
        reduce_columns, is_dir)
from pymsbayes.utils.errors import WorkerExecutionError, PriorMergeError
from pymsbayes.utils.stats import *
from pymsbayes.utils.parsing import *
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)


##############################################################################
## globals

VALID_REJECTION_TOOLS = ['msreject', 'abctoolbox']
VALID_REGRESSION_METHODS = ['llr', 'glm']

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
        self.tag = str(tag)

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
    if not workers:
        raise Exception('no prior workers to merge.')
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
            summary_worker = None,
            rejection_worker = None,
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
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.temp_fs = temp_fs
        self.sample_size = int(sample_size)
        self.config_path = expand_path(config_path)
        self.output_dir = self.temp_fs.create_subdir(prefix = self.name + '-')
        if not exe_path:
            exe_path = get_tool_path('msbayes-old')
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
        self.summary_worker = summary_worker
        self.rejection_worker = rejection_worker
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
        if self.summary_worker:
            self.summary_worker.start()
        if self.rejection_worker:
            self.rejection_worker.start()

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
        self.name = self.__class__.__name__ + '-' + str(self.count)
        if not exe_path:
            exe_path = get_tool_path('msreject')
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
            regress_summary_path = None,
            regress_posterior_path = None,
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
        self.name = self.__class__.__name__ + '-' + str(self.count)
        if not exe_path:
            exe_path = get_tool_path('regress_cli')
        self.exe_path = expand_path(exe_path)
        self.observed_path = expand_path(observed_path)
        self.posterior_path = expand_path(posterior_path)
        self.tolerance = float(tolerance)
        self.header = None
        self.stat_indices = stat_indices
        self.continuous_parameter_indices = continuous_parameter_indices
        self.discrete_parameter_indices = discrete_parameter_indices
        if not regress_summary_path:
            regress_summary_path = self.posterior_path + \
                    '.regression-summary.txt'
        self.regress_summary_path = regress_summary_path
        if not regress_posterior_path:
            regress_posterior_path = self.posterior_path + \
                    '.regression-adjusted.txt'
        self.regress_posterior_path = regress_posterior_path

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
                    parameter_patterns=(MEAN_TAU_PATTERNS + OMEGA_PATTERNS))
        if not self.discrete_parameter_indices:
            self.discrete_parameter_indices = get_parameter_indices(
                    self.header,
                    parameter_patterns=(MODEL_PATTERNS + PSI_PATTERNS + \
                            DIV_MODEL_PATTERNS))
        self._update_cmd()

    def _update_cmd(self):
        cmd = [self.exe_path,
               '-t', str(self.tolerance),
               '--observed-path={0}'.format(self.observed_path),
               '--posterior-path={0}'.format(self.posterior_path),
               '--summary-path={0}'.format(self.regress_summary_path),
               '--adjusted-path={0}'.format(self.regress_posterior_path),
               '-c', ','.join(
                    [str(i+1) for i in self.continuous_parameter_indices]),
               '-d', ','.join(
                    [str(i+1) for i in self.discrete_parameter_indices]),
               '-s', ','.join(
                    [str(i+1) for i in self.stat_indices]),]
        self.cmd = cmd


class EuRejectSummaryMerger(object):
    count = 0
    def __init__(self, eureject_workers, tag = ''):
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.eureject_workers = eureject_workers
        self.sample_sum_collection = None
        self.finished = False
        self.tag = str(tag)

    def _check_worker(self, eureject_worker):
        if not eureject_worker.finished:
            raise Exception('{0} has not finished running'.format(
                    eureject_worker.name))
        if eureject_worker.num_summarized < 1:
            raise Exception('{0} did not summarize any samples'.format(
                    eureject_worker.name))
        if eureject_worker.num_summarized != \
                eureject_worker.num_standardizing_samples:
            raise Exception('{0} summarized {1} samples, but was expected '
                    'to summarize {2}'.format(eureject_worker.name,
                            eureject_worker.num_summarized,
                            eureject_worker.num_standardizing_samples))

    def start(self):
        for ew in self.eureject_workers:
            self._check_worker(ew)
            ssc = SampleSummaryCollection.get_from_summary_file(
                    ew.summary_out_path)
            assert(ssc.sample_sums[ssc.keys[0]].n == ew.num_summarized)
            if not self.sample_sum_collection:
                self.sample_sum_collection = ssc
            else:
                self.sample_sum_collection.update(ssc)
        self.finished = True
    
    def write_summary(self, path):
        self.sample_sum_collection.write(path)


class EuRejectWorker(Worker):
    count = 0
    stderr_pattern_string = (
            "Files\s+used\s+for\s+calculating[^:]*:\s*"
            "(?P<standardizing_files>[^\n]+)\n"
            "Number\s+of\s+samples\s+used\s+for[^:]*:\s*"
            "(?P<num_summarized>\d+)\s*\n"
            "Files\s+processed[^:]*:\s*"
            "(?P<rejection_files>[^\n]+)\n"
            "Total\s+number\s+of\s+samples\s+processed[^:]*:\s*"
            "(?P<num_processed>\d+)\s*\n"
            "Number\s+of\s+samples\s+retained[^:]*:\s*"
            "(?P<num_retained>\d+)\s*\n")
    stderr_pattern = re.compile(stderr_pattern_string, re.IGNORECASE)
    def __init__(self,
            temp_fs,
            observed_path,
            prior_paths,
            num_posterior_samples,
            num_standardizing_samples = 10000,
            summary_in_path = None,
            summary_out_path = None,
            posterior_path = None,
            regression_worker = None,
            exe_path = None,
            stderr_path = None,
            keep_temps = False,
            tag = ''):
        Worker.__init__(self,
                stdout_path = None,
                stderr_path = None,
                tag = tag)
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.temp_fs = temp_fs
        if not exe_path:
            exe_path = get_tool_path('eureject')
        self.exe_path = expand_path(exe_path)
        self.output_dir = self.temp_fs.create_subdir(prefix = self.name + '-')
        if not stderr_path:
            stderr_path = self.temp_fs.get_file_path(
                    parent = self.output_dir,
                    prefix = self.name + '-stderr-',
                    create = False)
        self.stderr_path = expand_path(stderr_path)
        if not summary_in_path:
            self.summary_provided = False
        else:
            self.summary_provided = True
            self.summary_in_path = expand_path(summary_in_path)
        if not summary_out_path:
            summary_out_path = self.temp_fs.get_file_path(
                    parent = self.output_dir,
                    prefix = self.name + '-summary-',
                    create = False)
        self.summary_out_path = expand_path(summary_out_path)
        self.observed_path = expand_path(observed_path)
        self.prior_paths = [expand_path(p) for p in prior_paths]
        if not posterior_path:
            posterior_path = observed_path + '.posterior'
        self.posterior_path = expand_path(posterior_path)
        self.stdout_path = posterior_path
        self.num_posterior_samples = int(num_posterior_samples)
        self.num_standardizing_samples = int(num_standardizing_samples)
        self.regression_worker = regression_worker
        self.keep_temps = keep_temps
        self.num_retained = None
        self.num_summarized = None
        self.num_processed = None
        self.standardizing_files = None
        self.rejection_files = None
        self._update_cmd()

    def _pre_process(self):
        self._update_cmd()

    def _post_process(self):
        self._parse_stderr()
        if self.regression_worker:
            self.regression_worker.start()
        if not self.keep_temps:
            self.temp_fs.remove_dir(self.output_dir)

    def _parse_stderr(self):
        se = self.get_stderr()
        if not se:
            raise WorkerExecutionError('unexpected std error from '
                    '{0}: {1}'.format(self.name, se))
        m = self.stderr_pattern.search(se)
        if m is None:
            raise WorkerExecutionError('unexpected std error from '
                    '{0}: {1}'.format(self.name, se))
        self.num_retained = int(m.group('num_retained'))
        self.num_summarized = int(m.group('num_summarized'))
        self.num_processed = int(m.group('num_processed'))
        self.rejection_files = [x.strip() for x in m.group(
                'rejection_files').split(',')]
        self.standardizing_files = [x.strip() for x in m.group(
                'standardizing_files').split(',')]
        if hasattr(se, 'close'):
            se.close()

    def _update_cmd(self):
        cmd = [self.exe_path,
               '-f', self.observed_path,
               '-k', str(self.num_posterior_samples),
               '-n', str(self.num_standardizing_samples),
               '-o', self.summary_out_path]
        if self.summary_provided:
            cmd.extend(['-s', self.summary_in_path])
        cmd.extend(self.prior_paths)
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
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.temp_fs = temp_fs
        if not exe_path:
            exe_path = get_tool_path('abcestimator')
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
            regress_summary_path = None,
            regress_posterior_path = None,
            exe_path = None,
            stdout_path = None,
            stderr_path = None,
            keep_temps = False,
            bandwidth = None,
            num_posterior_samples = None,
            num_posterior_quantiles = 1000,
            compress = False,
            tag = '',
            ):
        Worker.__init__(self,
                stdout_path = stdout_path,
                stderr_path = stderr_path,
                tag = tag)
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.temp_fs = temp_fs
        if not exe_path:
            exe_path = get_tool_path('abcestimator')
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
        self.num_posterior_quantiles = int(num_posterior_quantiles)
        self.num_posterior_samples = num_posterior_samples
        self.bandwidth = bandwidth
        if not regress_summary_path:
            regress_summary_path = self.posterior_path + \
                    '.regression-summary.txt'
        self.regress_summary_path = regress_summary_path
        if not regress_posterior_path:
            regress_posterior_path = self.posterior_path + \
                    '.regression-adjusted.txt'
        self.regress_posterior_path = regress_posterior_path
        self.compress = bool(compress)

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
                            PSI_PATTERNS + \
                            DIV_MODEL_PATTERNS)))
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
        regress_summary_path = self.output_prefix + \
                'PosteriorCharacteristics_Obs0.txt'
        regress_posterior_path = self.output_prefix + \
                'PosteriorEstimates_Obs0.txt'
        if self.compress:
            summary_stream = open(regress_summary_path, 'rU')
            out = GzipFileStream(self.regress_summary_path, 'w',
                    compresslevel = 9)
            for line in summary_stream:
                out.write(line)
            summary_stream.close()
            out.close()
            adjusted_stream = open(regress_posterior_path, 'rU')
            out = GzipFileStream(self.regress_posterior_path, 'w',
                    compresslevel = 9)
            for line in adjusted_stream:
                out.write(line)
            adjusted_stream.close()
            out.close()
        else:
            shutil.move(regress_summary_path, self.regress_summary_path)
            shutil.move(regress_posterior_path, self.regress_posterior_path)
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

class PosteriorWorker(object):
    count = 0
    def __init__(self,
            temp_fs,
            observed_path,
            posterior_path,
            num_taxon_pairs,
            posterior_out_path = None,
            output_prefix = None,
            model_indices = None,
            regress_posterior_path = None,
            abctoolbox_exe_path = None,
            abctoolbox_stdout_path = None,
            abctoolbox_stderr_path = None,
            abctoolbox_bandwidth = None,
            abctoolbox_num_posterior_quantiles = None,
            omega_threshold = 0.01,
            compress = False,
            keep_temps = False,
            tag = ''):
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.temp_fs = temp_fs
        self.temp_output_dir = self.temp_fs.create_subdir(
                prefix = self.name + '-')
        self.temp_posterior_path = self.temp_fs.get_file_path(
                parent = self.temp_output_dir,
                prefix = 'temp-posterior-')
        self.posterior_path = expand_path(posterior_path)
        self.observed_path = expand_path(observed_path)
        self.num_taxon_pairs = int(num_taxon_pairs)
        if not posterior_out_path:
            posterior_out_path = self.posterior_path
        self.posterior_out_path = expand_path(posterior_out_path)
        if not output_prefix:
            output_prefix = self.posterior_out_path
        self.output_prefix = output_prefix
        self.div_model_results_path = self.output_prefix + \
                '-div-model-results.txt'
        self.psi_results_path = self.output_prefix + '-psi-results.txt'
        self.model_results_path = self.output_prefix + '-model-results.txt'
        self.omega_results_path = self.output_prefix + '-omega-results.txt'
        self.posterior_summary_path = self.output_prefix + \
                '-posterior-summary.txt'
        self.regress_summary_path = self.output_prefix + \
                '-glm-posterior-summary.txt'
        if not regress_posterior_path:
            regress_posterior_path = self.output_prefix + \
                    '-glm-posterior-density-estimates.txt'
        self.regress_posterior_path = regress_posterior_path
        self.compress = bool(compress)
        if self.compress:
            self.regress_summary_path += '.gz'
            self.regress_posterior_path += '.gz'
            self.posterior_out_path += '.gz'
        if model_indices is not None:
            model_indices = set(model_indices)
        self.model_indices = model_indices
        self.finished = False
        self.tag = str(tag)
        self.top_div_models_to_indices = {}
        self.num_posterior_samples = None
        self.regression_worker = None
        self.keep_temps = bool(keep_temps)
        self.abctoolbox_exe_path = abctoolbox_exe_path
        self.abctoolbox_stdout_path = abctoolbox_stdout_path
        self.abctoolbox_stderr_path = abctoolbox_stderr_path
        self.abctoolbox_bandwidth = abctoolbox_bandwidth
        self.abctoolbox_num_posterior_quantiles = \
                abctoolbox_num_posterior_quantiles
        self.omega_threshold = float(omega_threshold)
        self.psi_probs = dict(zip(
                [i+1 for i in range(self.num_taxon_pairs)],
                [None for i in range(self.num_taxon_pairs)]))
        self.adjusted_psi_probs = dict(zip(
                [i+1 for i in range(self.num_taxon_pairs)],
                [None for i in range(self.num_taxon_pairs)]))
        self.model_probs = None
        self.adjusted_model_probs = None
        self.prob_omega_zero = None
        self.adjusted_prob_omega_zero = None
        self.div_model_summary = None
        self.unadjusted_summaries = {}

    def _process_posterior_sample(self):
        post = parse_parameters(self.posterior_path)
        if not post.has_key('taus'):
            raise Exception('posterior sample in {0} does not contain a '
                    'divergence time vector'.format(self.posterior_path))
        div_models = IntegerPartitionCollection(post['taus'])
        self.num_posterior_samples = div_models.n
        self.div_model_summary = div_models.get_summary()
        self._map_top_div_models(div_models)
        dmodels = []
        for k, dm in div_models.iteritems():
            dmodels.extend([self.top_div_models_to_indices[k]] * dm.n)
        post['div_model'] = dmodels
        psi_freqs = get_freqs(post['psi'])
        if max(psi_freqs.iterkeys()) > self.num_taxon_pairs:
            raise ValueError('number of taxon pairs is {0}, but found '
                    'psi estimates of {1} in posterior {2}'.format(
                        self.num_taxon_pairs,
                        max(psi_freqs.iterkeys()),
                        self.posterior_path))
        for i in range(self.num_taxon_pairs):
            self.psi_probs[i+1] = psi_freqs.get(i+1, 0.0)
        if post.has_key('model'):
            model_freqs = get_freqs(post['model'])
            if self.model_indices:
                self.model_probs = {}
                if not set(model_freqs.keys()).issubset(self.model_indices):
                    raise ValueError('model indices in posterior ({0}) are not '
                            'a subset of indices provided ({1})'.format(
                                ','.join([str(i) for i in model_freqs.keys()]),
                                ','.join([str(i) for i in self.model_indices])))
                for i in self.model_indices:
                    self.model_probs[i] = model_freqs.get(i, 0.0)
            else:
                self.model_probs = model_freqs
        self.prob_omega_zero = freq_less_than(post['omega'],
                self.omega_threshold)
        self.unadjusted_summaries['PRI.omega'] = get_summary(post['omega'])
        self.unadjusted_summaries['PRI.E.t'] = get_summary(post['mean_tau'])
        self.unadjusted_summaries['PRI.Psi'] = get_summary(post['psi'])
        if post.has_key('model'):
            self.unadjusted_summaries['PRI.model'] = get_summary(post['model'])
        self.unadjusted_summaries['PRI.div.model'] = get_summary(
                post['div_model'])

    def _map_top_div_models(self, div_models):
        for i, k in enumerate(div_models.iterkeys()):
            self.top_div_models_to_indices[k] = i + 1

    def _add_div_model_column_to_posterior(self):
        add_div_model_column(self.posterior_path, self.temp_posterior_path,
                self.top_div_models_to_indices,
                compresslevel = None)

    def _return_posterior_sample(self):
        if self.compress:
            post_stream = open(self.temp_posterior_path, 'rU')
            out = GzipFileStream(self.posterior_out_path, 'w',
                    compresslevel = 9)
            for line in post_stream:
                out.write(line)
            out.close()
            post_stream.close()
            if self.posterior_out_path == self.posterior_path + '.gz':
                os.remove(self.posterior_path)
        else:
            shutil.move(self.temp_posterior_path, self.posterior_out_path)

    def _prep_regression_worker(self):
        self.regression_worker = ABCToolBoxRegressWorker(
                temp_fs = self.temp_fs,
                observed_path = self.observed_path,
                posterior_path = self.temp_posterior_path,
                parameter_indices = None,
                regress_summary_path = self.regress_summary_path,
                regress_posterior_path = self.regress_posterior_path,
                exe_path = self.abctoolbox_exe_path,
                stdout_path = self.abctoolbox_stdout_path,
                stderr_path = self.abctoolbox_stderr_path,
                keep_temps = self.keep_temps,
                bandwidth = self.abctoolbox_bandwidth,
                num_posterior_samples = self.num_posterior_samples,
                num_posterior_quantiles = \
                        self.abctoolbox_num_posterior_quantiles,
                compress = self.compress)

    def _process_regression_results(self):
        discrete_probs = summarize_discrete_parameters_from_densities(
                self.regress_posterior_path,
                discrete_parameter_patterns = PSI_PATTERNS + MODEL_PATTERNS + \
                        DIV_MODEL_PATTERNS,
                include_omega_summary = True,
                omega_threshold = self.omega_threshold)
        if max(discrete_probs['PRI.Psi'].iterkeys()) > self.num_taxon_pairs:
            raise ValueError('number of taxon pairs is {0}, but found '
                    'psi estimates of {1} in posterior {2}'.format(
                        self.num_taxon_pairs,
                        max(discrete_probs['PRI.Psi'].iterkeys()),
                        self.regress_posterior_path))
        for i in range(self.num_taxon_pairs):
            self.adjusted_psi_probs[i+1] = discrete_probs['PRI.Psi'].get(i+1, 0.0)
        for k, s in self.div_model_summary:
            s['adjusted_frequency'] = discrete_probs['PRI.div.model'].get(
                    self.top_div_models_to_indices[k], 0.0)
        self.adjusted_prob_omega_zero = discrete_probs['PRI.omega'].get(0, 0.0)
        if self.model_probs:
            self.adjusted_model_probs = {}
            for i in self.model_probs.iterkeys():
                self.adjusted_model_probs[i] = discrete_probs['PRI.model'].get(
                        i, 0.0)

    def _write_div_model_results(self):
        out, close = process_file_arg(self.div_model_results_path, 'w')
        out.write('divergence_model\testimated_prob\tglm_adjusted_prob\t'
                'div_model_with_conditional_age_estimates\n')
        for div_model, summary in self.div_model_summary:
            out.write('{0}\t{1}\t{2}\t{3}\n'.format(
                    div_model,
                    summary['frequency'],
                    summary['adjusted_frequency'],
                    summary['string']))
        if close:
            out.close()

    def _write_psi_results(self):
        out, close = process_file_arg(self.psi_results_path, 'w')
        out.write('num_of_div_events\testimated_prob\t'
                'glm_adjusted_prob\n')
        assert sorted(self.psi_probs.keys()) == sorted(
                self.adjusted_psi_probs.keys())
        for k in self.psi_probs.iterkeys():
            out.write('{0}\t{1}\t{2}\n'.format(k, self.psi_probs[k],
                    self.adjusted_psi_probs[k]))
        if close:
            out.close()

    def _write_model_results(self):
        if not self.model_probs:
            return
        out, close = process_file_arg(self.model_results_path, 'w')
        out.write('model\testimated_prob\t'
                'glm_adjusted_prob\n')
        for k in self.model_probs.iterkeys():
            out.write('{0}\t{1}\t{2}\n'.format(k, self.model_probs[k],
                    self.adjusted_model_probs[k]))
        if close:
            out.close()

    def _write_omega_results(self):
        out, close = process_file_arg(self.omega_results_path, 'w')
        out.write('omega_thresh\tprob_less_than\t'
                'glm_prob_less_than\n')
        out.write('{0}\t{1}\t{2}\n'.format(
                self.omega_threshold,
                self.prob_omega_zero,
                self.adjusted_prob_omega_zero))
        if close:
            out.close()

    def _write_summary(self):
        glm = parse_abctoolbox_summary_file(self.regress_summary_path)
        out, close = process_file_arg(self.posterior_summary_path, 'w')
        for param, summary in self.unadjusted_summaries.iteritems():
            out.write('[{0}]\n'.format(param))
            if isinstance(summary['modes'][0], tuple):
                out.write('    modes = {0}\n'.format(', '.join([
                    '({0}, {1})'.format(x, y) for x, y in summary['modes']])))
            else:
                out.write('    modes = {0}\n'.format(', '.join(
                    [str(x) for x in summary['modes']])))
            out.write('    mode_glm = {0}\n'.format(glm[param]['mode']))
            out.write('    median = {0}\n'.format(summary['median']))
            out.write('    median_glm = {0}\n'.format(glm[param]['median']))
            out.write('    mean = {0}\n'.format(summary['mean']))
            out.write('    mean_glm = {0}\n'.format(glm[param]['mean']))
            out.write('    n = {0}\n'.format(summary['n']))
            out.write('    range = {0}, {1}\n'.format(summary['range'][0],
                    summary['range'][1]))
            out.write('    HPD_95_interval = {0}, {1}\n'.format(
                    summary['hpdi_95'][0],
                    summary['hpdi_95'][1]))
            out.write('    HPD_95_interval_glm = {0}, {1}\n'.format(
                    glm[param]['HPD_95_lower_bound'],
                    glm[param]['HPD_95_upper_bound']))
            out.write('    quantile_95_interval = {0}, {1}\n'.format(
                    summary['qi_95'][0],
                    summary['qi_95'][1]))
            out.write('    quantile_95_interval_glm = {0}, {1}\n'.format(
                    glm[param]['quantile_95_lower_bound'],
                    glm[param]['quantile_95_upper_bound']))
        if close:
            out.close()

    def _post_process(self):
        if not self.keep_temps:
            self.temp_fs.remove_dir(self.temp_output_dir)

    def start(self):
        self._process_posterior_sample()
        self._add_div_model_column_to_posterior()
        self._prep_regression_worker()
        self.regression_worker.start()
        self._return_posterior_sample()
        self._process_regression_results()
        self._write_div_model_results()
        self._write_psi_results()
        self._write_model_results()
        self._write_omega_results()
        self._write_summary()
        self._post_process()
        self.finished = True

