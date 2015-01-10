#! /usr/bin/env python

import os
import sys
import re
import subprocess
import time
import shutil
import traceback
import tempfile
from cStringIO import StringIO

from pymsbayes.fileio import (expand_path, process_file_arg, FileStream, open,
        GzipFileStream)
from pymsbayes.config import MsBayesConfig
from pymsbayes.utils.tempfs import TempFileSystem
from pymsbayes.utils import (ToolPathManager, MSBAYES_SORT_INDEX,
        GLOBAL_RNG)
from pymsbayes.utils import probability
from pymsbayes.utils import functions
from pymsbayes.utils import errors
from pymsbayes.utils import stats
from pymsbayes.utils import parsing
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
            append_stdout = True,
            subprocess_kwargs = {},
            tag = None):
        self.__class__.total += 1
        self.stdout_path = stdout_path
        self.stderr_path = stderr_path
        self.cmd = []
        self.process = None
        self.finished = False
        self.exit_code = None
        self.stdout = None
        self.stderr = None
        self.append_stdout = append_stdout
        self.error = None
        self.trace_back = None
        self.subprocess_kwargs = subprocess_kwargs
        self.tag = tag

    def get_stderr(self):
        if not self.stderr_path:
            return self.stderr
        try:
            se = open(self.stderr_path, 'rU')
        except IOError, e:
            _LOG.error('Could not open stderr file')
            raise e
        msg = se.read()
        se.close()
        return msg

    def start(self):
        try:
            self._start()
        except Exception as e:
            self.error = e
            f = StringIO()
            traceback.print_exc(file=f)
            self.trace_back = f.getvalue()

    def _start(self):
        try:
            self._pre_process()
        except:
            e = StringIO()
            traceback.print_exc(file=e)
            _LOG.error('Error during pre-processing:\n{0}'.format(
                    e.getvalue()))
            raise
        _LOG.debug('Starting process with following command:\n\t'
                '{0}'.format(' '.join([str(x) for x in self.cmd])))
        if self.stdout_path:
            if self.append_stdout:
                sout = open(self.stdout_path, 'a')
            else:
                sout = open(self.stdout_path, 'w')
        else:
            sout = subprocess.PIPE
        if self.stderr_path:
            serr = open(self.stderr_path, 'w')
        else:
            serr = subprocess.PIPE
        self.process = subprocess.Popen(self.cmd,
                stdout = sout,
                stderr = serr,
                shell = False,
                **self.subprocess_kwargs)
        self.stdout, self.stderr = self.process.communicate()
        self.exit_code = self.process.wait()
        if hasattr(sout, 'close'):
            sout.close()
        if hasattr(serr, 'close'):
            serr.close()
        if self.exit_code != 0:
            _LOG.error('execution failed')
            raise errors.WorkerExecutionError('{0} failed.\ninvocation:\n{1}\n'
                    'stderr:\n{2}'.format(
                    self.name, ' '.join([str(x) for x in self.cmd]),
                    self.get_stderr()))
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


class GenericWorker(Worker):
    count = 0
    def __init__(self,
            exe,
            args,
            stdout_path = None,
            stderr_path = None,
            append_stdout = True,
            subprocess_kwargs = {},
            tag = None):
        Worker.__init__(self,
                stdout_path = stdout_path,
                stderr_path = stderr_path,
                append_stdout = append_stdout,
                subprocess_kwargs = subprocess_kwargs,
                tag = tag)
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.exe = ToolPathManager.get_tool_path(exe)
        if isinstance(args, str):
            args = args.strip().split()
        self.args = [str(x) for x in args]
        self.cmd = [self.exe] + self.args

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
            raise errors.PriorMergeError('Workers {0} and {1} have different '
                    'prior headers. Cannot merge!'.format(std.name, w.name))
        prior_file = open(w.prior_path, 'rU')
        for line in prior_file:
            if not parsing.HEADER_PATTERN.match(line.strip()):
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

def merge_prior_files(paths, dest_path, append = True, compresslevel = None):
    h = None
    header_exists = False
    ncols = None
    mode = 'a'
    if not append:
        mode = 'w'
    if append and os.path.exists(dest_path):
        h = parsing.parse_header(dest_path, strict = False)
        if h:
            header_exists = True
            ncols = len(h)
        else:
            dest_stream, d_close = process_file_arg(dest_path, mode = 'rU',
                    compresslevel = compresslevel)
            try:
                ncols = len(dest_stream.next().strip().split())
            except StopIteration:
                ncols = None
            finally:
                if d_close:
                    dest_stream.close()
    out, out_close = process_file_arg(dest_path, mode = mode,
            compresslevel = compresslevel)
    std = None
    if ncols:
        std = out
    for p in paths:
        header = None
        f, f_close = process_file_arg(p, 'rU')
        for line_num, line in enumerate(f):
            if line_num == 0:
                if parsing.HEADER_PATTERN.match(line.strip()):
                    header = line.strip().split()
                    if h and h != header:
                        if out_close:
                            out.close()
                        if f_close:
                            f.close()
                        raise errors.PriorMergeError('prior files {0} and {1} '
                                'have different headers'.format(f.name,
                                        std.name))
                if not ncols:
                    ncols = len(line.strip().split())
                    std = f
                    if header:
                        if not h:
                            h = header
                        if not header_exists:
                            out.write(line)
                        continue
                if len(line.strip().split()) != ncols:
                    if out_close:
                        out.close()
                    if f_close:
                        f.close()
                    raise errors.PriorMergeError('prior files {0} and {1} do '
                            'not have the same number of columns. Cannot '
                            'merge!'.format(f.name, std.name))
                if not header:
                    out.write(line)
            else:
                out.write(line)
        if f_close:
            f.close()
    if out_close:
        out.close()

##############################################################################
## msBayes class for generating prior files

class ObsSumStatsWorker(Worker):
    valid_schemas = ['abctoolbox']
    count = 0
    def __init__(self,
            temp_fs,
            config_path,
            output_path = None,
            exe_path = None,
            schema = 'abctoolbox',
            stat_patterns = parsing.DEFAULT_STAT_PATTERNS,
            stderr_path = None,
            keep_temps = False,
            tag = None):
        Worker.__init__(self,
                stdout_path = None,
                stderr_path = stderr_path,
                tag = tag)
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.temp_fs = temp_fs
        self.config_path = expand_path(config_path)
        if not output_path:
            output_path = self.config_path + '-summary-statistics.txt'
        self.output_path = expand_path(output_path)
        if not exe_path:
            exe_path = 'obsSumStats.pl'
        self.exe_path = ToolPathManager.get_tool_path(exe_path)
        self.sort_index = MSBAYES_SORT_INDEX.current_value()
        self.schema = schema.lower()
        if not self.schema in self.valid_schemas:
            raise ValueError('{0} is not a valid schema, options are: '
                    '{1}'.format(self.schema, ', '.join(self.valid_schemas)))
        self.stat_patterns = stat_patterns
        self.keep_temps = bool(keep_temps)
        self.header = None
        self.output_dir = None
        self.temp_output_path = None

    def _pre_process(self):
        self.output_dir = self.temp_fs.create_subdir(prefix = self.name + '-')
        self.temp_output_path = self.temp_fs.get_file_path(
                parent = self.output_dir,
                prefix = self.name + '-temp-summary-stats-')
        self.stdout_path = self.temp_output_path
        self._update_cmd()

    def _update_cmd(self):
        cmd = [self.exe_path,
               '-s', str(self.sort_index)]
        cmd.append(self.config_path)
        self.cmd = cmd

    def _post_process(self):
        _LOG.debug('{0}\n'.format(self.get_stderr()))
        if not os.path.isfile(self.temp_output_path):
            raise Exception('{0} did not produce the output file {1}; '
                    'here is the std error:\n{2}\n'.format(
                        self.name, self.temp_output_path, self.get_stderr()))
        self.header = parsing.observed_stats_for_abctoolbox(
                    in_file = self.temp_output_path,
                    out_file = self.output_path,
                    stat_patterns = self.stat_patterns)
        if not self.keep_temps:
            self.temp_fs.remove_dir(self.output_dir)


class MsBayesWorker(Worker):
    count = 0
    valid_schemas = ['msreject', 'abctoolbox']

    def __init__(self,
            temp_fs,
            sample_size,
            config_path,
            prior_path = None,
            header_path = None,
            exe_path = None,
            model_index = None,
            report_parameters = True,
            seed = None,
            schema = 'msreject',
            include_header = False,
            stat_patterns=parsing.DEFAULT_STAT_PATTERNS,
            parameter_patterns=parsing.PARAMETER_PATTERNS,
            summary_worker = None,
            rejection_worker = None,
            stdout_path = None,
            stderr_path = None,
            staging_dir = None,
            write_stats_file = False,
            write_header_file = True,
            stats_file_path = None,
            tag = None):
        Worker.__init__(self,
                stdout_path = stdout_path,
                stderr_path = stderr_path,
                tag = tag)
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.temp_fs = temp_fs
        self.output_dir = self.temp_fs.create_subdir(prefix = self.name + '-')
        self.sample_size = int(sample_size)
        if self.sample_size < 1:
            raise ValueError('sample size {0} is not a positive '
                    'integer'.format(sample_size))
        self.config_path = expand_path(config_path)
        if not exe_path:
            cfg = MsBayesConfig(self.config_path)
            if cfg.implementation.lower() == 'new':
                exe_path = 'dpp-msbayes.pl'
            else:
                exe_path = 'msbayes.pl'
        self.exe_path = ToolPathManager.get_tool_path(exe_path)
        self.model_index = None
        if model_index != None:
            self.model_index = int(model_index)
        self.sort_index = MSBAYES_SORT_INDEX.current_value()
        self.report_parameters = report_parameters
        if seed is None:
            self.seed = functions.get_random_int()
        else:
            self.seed = int(seed)
        if prior_path:
            self.prior_path = expand_path(prior_path)
        else:
            self.prior_path = self.temp_fs.get_file_path(
                    parent = self.output_dir,
                    prefix = 'prior-{0}-{1}-'.format(
                            self.sample_size,
                            self.seed))
        self.header_path = None
        if header_path:
            self.header_path = expand_path(header_path)
        if not schema.lower() in self.valid_schemas:
            raise ValueError(
                    'schema {0} is not valid. Options are: {1}'.format(
                        schema, ','.join(self.valid_schemas)))
        self.schema = schema.lower()
        self.include_header = include_header
        self.prior_stats_path = None
        self.write_stats_file = write_stats_file
        if self.write_stats_file:
            if stats_file_path:
                self.prior_stats_path = expand_path(stats_file_path)
            else:
                self.prior_stats_path = self.prior_path + '.stats.txt'
                self.temp_fs._register_file(self.prior_stats_path)
        self.write_header_file = write_header_file
        self.stat_patterns = stat_patterns
        self.parameter_patterns = parameter_patterns
        self.header = None
        self.parameter_indices = None
        self.stat_indices = None
        self.staging_dir = None
        self.staging_prior_path = None
        if functions.is_dir(staging_dir):
            self.staging_dir = staging_dir
            self.staging_prior_path = os.path.join(self.staging_dir,
                    os.path.basename(self.prior_path))
            self.staging_prior_stats_path = self.staging_prior_path
            if self.write_stats_file:
                self.staging_prior_stats_path = self.staging_prior_path + \
                        '.stats.txt'
        self.summary_worker = summary_worker
        self.rejection_worker = rejection_worker
        self._update_cmd()

    def _pre_process(self):
        self._update_cmd()

    def _update_cmd(self):
        cmd = [self.exe_path,
               '-r', str(self.sample_size),
               '-S', str(self.seed),
               '-s', str(self.sort_index),
               '-c', self.config_path,]
        if self.staging_dir:
            cmd.extend(['-o', self.staging_prior_path])
        else:
            cmd.extend(['-o', self.prior_path])
        if self.model_index != None:
            cmd.extend(['-m', str(self.model_index)])
        if self.report_parameters:
            cmd.append('-p')
        self.cmd = cmd

    def _post_process(self):
        _LOG.debug('{0}\n'.format(self.get_stderr()))
        prior_path = self.prior_path
        prior_stats_path = self.prior_stats_path
        if self.staging_dir:
            prior_path = self.staging_prior_path
            prior_stats_path = self.staging_prior_stats_path
        raw_prior_path = prior_path + '.raw'
        shutil.move(prior_path, raw_prior_path)
        if self.prior_stats_path:
            h = parsing.observed_stats_for_abctoolbox(
                    in_file = raw_prior_path,
                    out_file = prior_stats_path,
                    stat_patterns = self.stat_patterns)
        if self.schema == 'msreject':
            header = parsing.prior_for_msreject(
                    in_file = raw_prior_path,
                    out_file = prior_path,
                    stat_patterns = self.stat_patterns,
                    parameter_patterns = self.parameter_patterns,
                    dummy_patterns = parsing.DUMMY_PATTERNS,
                    include_header = self.include_header)
        elif self.schema == 'abctoolbox':
            header = parsing.prior_for_abctoolbox(
                    in_file = raw_prior_path,
                    out_file = prior_path,
                    stat_patterns = self.stat_patterns,
                    parameter_patterns = self.parameter_patterns)
        else:
            raise ValueError(
                    'schema {0} is not valid. Options are: {1}'.format(
                        self.schema, ','.join(self.valid_schemas)))
        os.remove(raw_prior_path)
        if header and self.write_header_file:
            self.header_path = self.temp_fs.get_file_path(
                    parent = self.output_dir,
                    prefix = 'prior-{0}-{1}-header-'.format(
                            self.sample_size,
                            self.seed))
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
        self.parameter_indices = parsing.get_parameter_indices(
                header_list = self.header,
                parameter_patterns = parsing.PARAMETER_PATTERNS)
        self.stat_indices = parsing.get_stat_indices(
                header_list = self.header,
                stat_patterns = parsing.ALL_STAT_PATTERNS)

    def purge(self):
        self.temp_fs.remove_dir(self.output_dir)

        
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
        num_prior_samples = parsing.line_count(prior_path, ignore_headers=True)
    if tolerance is None:
        tolerance = num_posterior_samples / float(num_prior_samples)
    if num_posterior_samples is None:
        num_posterior_samples = int(tolerance * num_prior_samples)
    if num_prior_samples is None:
        num_prior_samples = int(num_posterior_samples / tolerance)
    obs_file, close = process_file_arg(observed_sims_file)
    header = parsing.parse_header(obs_file)
    if rejection_tool == 'abctoolbox':
        header = parsing.parse_header(prior_path)
    all_stat_indices = parsing.get_stat_indices(header,
            stat_patterns=parsing.ALL_STAT_PATTERNS)
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
                prefix = 'observed-{0}-reject-'.format(i+1))
        reg_obs_path = temp_fs.get_file_path(parent = obs_temp_dir,
                prefix = 'observed-{0}-regress-'.format(i+1))
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
            tag = None):
        Worker.__init__(self,
                stdout_path = None,
                stderr_path = stderr_path,
                tag = tag)
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        if not exe_path:
            exe_path = 'msReject'
        self.exe_path = ToolPathManager.get_tool_path(exe_path)
        self.observed_path = expand_path(observed_path)
        self.prior_path = expand_path(prior_path)
        if not posterior_path:
            posterior_path = observed_path + '.posterior'
        self.posterior_path = expand_path(posterior_path)
        self.header = header
        potential_stat_indices = parsing.get_stat_indices(
                self.header,
                stat_patterns = parsing.ALL_STAT_PATTERNS)
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
            tag = None):
        Worker.__init__(self,
                stdout_path = stdout_path,
                stderr_path = stderr_path,
                tag = tag)
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        if not exe_path:
            exe_path = 'regress_cli.r'
        self.exe_path = ToolPathManager.get_tool_path(exe_path)
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
        self.header = parsing.parse_header(self.posterior_path)
        potential_stat_indices = parsing.get_stat_indices(
                self.header,
                stat_patterns = parsing.ALL_STAT_PATTERNS)
        if not self.stat_indices:
            self.stat_indices = potential_stat_indices
        else:
            diff = set(self.stat_indices) - set(potential_stat_indices)
            if len(diff) > 0:
                raise ValueError('stat indices are not valid')
        if not self.continuous_parameter_indices:
            self.continuous_parameter_indices = parsing.get_parameter_indices(
                    self.header,
                    parameter_patterns=(parsing.MEAN_TAU_PATTERNS +
                            parsing.OMEGA_PATTERNS +
                            parsing.CV_PATTERNS))
        if not self.discrete_parameter_indices:
            self.discrete_parameter_indices = parsing.get_parameter_indices(
                    self.header,
                    parameter_patterns=(parsing.MODEL_PATTERNS + \
                            parsing.PSI_PATTERNS + \
                            parsing.DIV_MODEL_PATTERNS))
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
    def __init__(self, eureject_workers, tag = None):
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.eureject_workers = eureject_workers
        self.sample_sum_collection = None
        self.finished = False
        self.tag = tag

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
            ssc = stats.SampleSummaryCollection.get_from_summary_file(
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
            tag = None):
        Worker.__init__(self,
                stdout_path = None,
                stderr_path = None,
                append_stdout = False,
                tag = tag)
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.temp_fs = temp_fs
        if not exe_path:
            exe_path = 'eureject'
        self.exe_path = ToolPathManager.get_tool_path(exe_path)
        self.stderr_path = None
        if stderr_path:
            self.stderr_path = expand_path(stderr_path)
        if not summary_in_path:
            self.summary_provided = False
        else:
            self.summary_provided = True
            self.summary_in_path = expand_path(summary_in_path)
        self.summary_out_path = None
        if summary_out_path:
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
        self.output_dir = None

    def _pre_process(self):
        self.output_dir = self.temp_fs.create_subdir(prefix = self.name + '-')
        if not self.stderr_path:
            self.stderr_path = self.temp_fs.get_file_path(
                    parent = self.output_dir,
                    prefix = self.name + '-stderr-')
        if not self.summary_out_path:
            self.summary_out_path = self.temp_fs.get_file_path(
                    parent = self.output_dir,
                    prefix = self.name + '-summary-')
        self._update_cmd()

    def _post_process(self):
        self._parse_stderr()
        if self.regression_worker:
            self.regression_worker.start()
        if not self.keep_temps:
            self.temp_fs.remove_dir(self.output_dir)

    def _parse_stderr(self):
        se = self.get_stderr()
        if se is None:
            raise errors.WorkerExecutionError('no std error from '
                    '{0}: {1}'.format(self.name, se))
        m = self.stderr_pattern.search(se)
        if m is None:
            raise errors.WorkerExecutionError('unexpected std error from '
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
            tag = None):
        Worker.__init__(self,
                stdout_path = stdout_path,
                stderr_path = stderr_path,
                tag = tag)
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.temp_fs = temp_fs
        if not exe_path:
            exe_path = 'ABCestimator'
        self.exe_path = ToolPathManager.get_tool_path(exe_path)
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
        self.keep_temps = keep_temps
        self.output_dir = None
        self.output_prefix = None
        self.cfg_path = None

    def _pre_process(self):
        self.output_dir = self.temp_fs.create_subdir(prefix = self.name + '-')
        self.output_prefix = os.path.join(self.output_dir, 
                self.temp_fs.token_id + '_ABC_GLM_')
        self.cfg_path = self.temp_fs.get_file_path(parent = self.output_dir,
                prefix = 'cfg-')
        self.header = parsing.parse_header(self.prior_path)
        self.parameter_indices = sorted(
                parsing.get_parameter_indices(self.header,
                        parameter_patterns=(parsing.PARAMETER_PATTERNS)) + \
                        parsing.get_dummy_indices(self.header))
        self.stats_header = parsing.parse_header(self.observed_path)
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
        if self.regression_worker:
            self.regression_worker.start()
        if not self.keep_temps:
            self.temp_fs.remove_dir(self.output_dir)

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
            tag = None,
            ):
        Worker.__init__(self,
                stdout_path = stdout_path,
                stderr_path = stderr_path,
                tag = tag)
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.temp_fs = temp_fs
        if not exe_path:
            exe_path = 'ABCestimator'
        self.exe_path = ToolPathManager.get_tool_path(exe_path)
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
        self.compress = compress
        self.keep_temps = keep_temps
        self.failed = False
        self.output_dir = None
        self.output_prefix = None
        self.cfg_path = None

    def _pre_process(self):
        self.output_dir = self.temp_fs.create_subdir(prefix = self.name + '-')
        self.output_prefix = os.path.join(self.output_dir, 
                self.temp_fs.token_id + '_ABC_GLM_')
        self.cfg_path = self.temp_fs.get_file_path(parent = self.output_dir,
                prefix = 'cfg-')
        self.header = parsing.parse_header(self.posterior_path)
        if not self.num_posterior_samples:
            self.num_posterior_samples = parsing.line_count(self.posterior_path,
                    ignore_headers=True)
        if not self.bandwidth:
            self.bandwidth = 2 / float(self.num_posterior_samples)
        potential_stat_indices = parsing.get_stat_indices(
                self.header,
                stat_patterns = parsing.ALL_STAT_PATTERNS)
        if not self.parameter_indices:
            self.parameter_indices = sorted(parsing.get_parameter_indices(
                    self.header,
                    parameter_patterns=(parsing.MEAN_TAU_PATTERNS + \
                            parsing.OMEGA_PATTERNS + \
                            parsing.CV_PATTERNS + \
                            parsing.MODEL_PATTERNS + \
                            parsing.PSI_PATTERNS + \
                            parsing.DIV_MODEL_PATTERNS)))
        if len(set.intersection(set(potential_stat_indices),
                set(self.parameter_indices))) > 0:
            raise ValueError('parameter indices are not valid. '
                    'they contain stat indices!')
        self.stats_header = parsing.parse_header(self.observed_path)
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
            if os.path.exists(regress_summary_path):
                try:
                    summary_stream = open(regress_summary_path, 'rU')
                    out = GzipFileStream(self.regress_summary_path, 'w',
                            compresslevel = 9)
                    for line in summary_stream:
                        out.write(line)
                    summary_stream.close()
                    out.close()
                except:
                    out.close()
                    _LOG.warning('{0}: problem compressing {1} to {2}... '
                            'trying without compression'.format(
                                self.name, regress_summary_path,
                                self.regress_summary_path))
                    base, ext = os.path.splitext(self.regress_summary_path)
                    if ext.lower() == '.gz':
                        self.regress_summary_path = base
                    shutil.move(regress_summary_path, self.regress_summary_path)
            else:
                _LOG.warning('{0}: did not produce summary file'.format(
                        self.name))
                self.regress_summary_path = None
                self.failed = True
            if os.path.exists(regress_posterior_path):
                try:
                    adjusted_stream = open(regress_posterior_path, 'rU')
                    out = GzipFileStream(self.regress_posterior_path, 'w',
                            compresslevel = 9)
                    for line in adjusted_stream:
                        out.write(line)
                    adjusted_stream.close()
                    out.close()
                except:
                    out.close()
                    _LOG.warning('{0}: problem compressing {1} to {2}... '
                            'trying without compression'.format(
                                self.name, regress_posterior_path,
                                self.regress_posterior_path))
                    base, ext = os.path.splitext(self.regress_posterior_path)
                    if ext.lower() == '.gz':
                        self.regress_posterior_path = base
                    shutil.move(regress_posterior_path,
                            self.regress_posterior_path)
            else:
                _LOG.warning('{0}: did not produce posterior density '
                        'file'.format(self.name))
                self.regress_summary_path = None
                self.failed = True
        else:
            if os.path.exists(regress_summary_path):
                shutil.move(regress_summary_path, self.regress_summary_path)
            else:
                _LOG.warning('{0}: did not produce summary file'.format(
                        self.name))
                self.regress_summary_path = None
                self.failed = True
            if os.path.exists(regress_posterior_path):
                shutil.move(regress_posterior_path, self.regress_posterior_path)
            else:
                _LOG.warning('{0}: did not produce posterior density '
                        'file'.format(self.name))
                self.regress_posterior_path = None
                self.failed = True
        if not self.keep_temps:
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
            cv_threshold = 0.01,
            compress = False,
            keep_temps = False,
            tag = None):
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.temp_fs = temp_fs
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
        self.cv_results_path = self.output_prefix + '-cv-results.txt'
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
        self.tag = tag
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
        self.cv_threshold = float(cv_threshold)
        self.cv_included = False
        self.psi_probs = dict(zip(
                [i+1 for i in range(self.num_taxon_pairs)],
                [None for i in range(self.num_taxon_pairs)]))
        self.adjusted_psi_probs = dict(zip(
                [i+1 for i in range(self.num_taxon_pairs)],
                [None for i in range(self.num_taxon_pairs)]))
        self.model_probs = None
        self.adjusted_model_probs = None
        self.prob_omega_zero = None
        self.prob_cv_zero = None
        self.adjusted_prob_omega_zero = None
        self.adjusted_prob_cv_zero = None
        self.div_model_summary = None
        self.unadjusted_summaries = {}
        self.parameter_patterns_to_remove = set()
        self.regression_failed = False
        self.temp_output_dir = None
        self.temp_posterior_path = None

    def _pre_process(self):
        self.temp_output_dir = self.temp_fs.create_subdir(
                prefix = self.name + '-')
        self.temp_posterior_path = self.temp_fs.get_file_path(
                parent = self.temp_output_dir,
                prefix = 'temp-posterior-')

    def _process_posterior_sample(self):
        post = parsing.parse_parameters(self.posterior_path)
        if not post.has_key('taus'):
            raise Exception('posterior sample in {0} does not contain a '
                    'divergence time vector'.format(self.posterior_path))
        if MSBAYES_SORT_INDEX.current_value() == 0:
            div_models = stats.PartitionCollection(post['taus'])
        else:
            div_models = stats.IntegerPartitionCollection(post['taus'])
        self.num_posterior_samples = div_models.n
        self.div_model_summary = div_models.get_summary()
        self._map_top_div_models(div_models)
        dmodels = []
        for k, dm in div_models.iteritems():
            dmodels.extend([self.top_div_models_to_indices[k]] * dm.n)
        post['div_model'] = dmodels
        psi_freqs = stats.get_freqs(post['psi'])
        if max(psi_freqs.iterkeys()) > self.num_taxon_pairs:
            raise ValueError('number of taxon pairs is {0}, but found '
                    'psi estimates of {1} in posterior {2}'.format(
                        self.num_taxon_pairs,
                        max(psi_freqs.iterkeys()),
                        self.posterior_path))
        for i in range(self.num_taxon_pairs):
            self.psi_probs[i+1] = psi_freqs.get(i+1, 0.0)
        if post.has_key('model'):
            if len(set(post['model'])) < 2:
                self.parameter_patterns_to_remove.update(parsing.MODEL_PATTERNS)
            model_freqs = stats.get_freqs(post['model'])
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
        self.prob_omega_zero = stats.freq_less_than(post['omega'],
                self.omega_threshold)
        if len(set(post['div_model'])) < 2:
            self.parameter_patterns_to_remove.update(parsing.DIV_MODEL_PATTERNS)
        if len(set(post['psi'])) < 2:
            self.parameter_patterns_to_remove.update(parsing.PSI_PATTERNS)
        if len(set(post['omega'])) < 2:
            self.parameter_patterns_to_remove.update(parsing.OMEGA_PATTERNS)
        if len(set(post['cv'])) < 2:
            self.parameter_patterns_to_remove.update(parsing.CV_PATTERNS)
        if len(set(post['mean_tau'])) < 2:
            self.cv_included = True
            self.parameter_patterns_to_remove.update(parsing.MEAN_TAU_PATTERNS)
        self.unadjusted_summaries['PRI.omega'] = stats.get_summary(
                post['omega'])
        if post.has_key('cv'):
            self.cv_included = True
            self.prob_cv_zero = stats.freq_less_than(post['cv'],
                    self.cv_threshold)
            self.unadjusted_summaries['PRI.cv'] = stats.get_summary(
                    post['cv'])
        self.unadjusted_summaries['PRI.E.t'] = stats.get_summary(
                post['mean_tau'])
        self.unadjusted_summaries['PRI.Psi'] = stats.get_summary(post['psi'])
        if post.has_key('model'):
            self.unadjusted_summaries['PRI.model'] = stats.get_summary(
                    post['model'])
        self.unadjusted_summaries['PRI.div.model'] = stats.get_summary(
                post['div_model'])
        if MSBAYES_SORT_INDEX.current_value() == 0:
            for i in range(len(post['taus'][0])):
                self.unadjusted_summaries['PRI.t.' + str(i+1)] = stats.get_summary(
                        [v[i] for v in post['taus']])

    def _map_top_div_models(self, div_models):
        for i, k in enumerate(div_models.iterkeys()):
            self.top_div_models_to_indices[k] = i + 1

    def _add_div_model_column_to_posterior(self):
        self.add_div_model_column(self.posterior_path,
                self.temp_posterior_path,
                self.top_div_models_to_indices,
                compresslevel = None)

    @classmethod
    def add_div_model_column(cls, in_file, out_file, div_models_to_indices,
            compresslevel = None):
        header = parsing.parse_header(in_file)
        if functions.get_indices_of_patterns(header,
                parsing.DIV_MODEL_PATTERNS) != []:
            raise errors.ParameterParsingError('posterior file {0} already has '
                    'a divergence model column'.format(
                    getattr(in_file, 'name', in_file)))
        header.insert(0, 'PRI.div.model')
        out, close = process_file_arg(out_file, 'w',
                compresslevel=compresslevel)
        out.write('{0}\n'.format('\t'.join(header)))
        other_index = max(div_models_to_indices.itervalues()) + 1
        for parameters, line in parsing.parameter_iter(in_file,
                include_line = True):
            if not parameters.has_key('taus'):
                out.close()
                raise errors.ParameterParsingError('posterior file {0} does not '
                        'contain divergence time vector'.format(
                        getattr(file_obj, 'name', file_obj)))
            if parsing.MSBAYES_SORT_INDEX.current_value() == 0:
                ip = stats.Partition(parameters['taus'][0])
            else:
                ip = stats.IntegerPartition(parameters['taus'][0])
            idx = div_models_to_indices.get(ip.key, other_index)
            line.insert(0, str(idx))
            out.write('{0}\n'.format('\t'.join(line)))
        if close:
            out.close()

    @classmethod
    def strip_div_model_column(cls, in_file, out_file, compresslevel = None):
        header = parsing.parse_header(in_file)
        div_indices = set(functions.get_indices_of_patterns(header,
                parsing.DIV_MODEL_PATTERNS))
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

    def _return_posterior_sample(self):
        if self.compress:
            try:
                post_stream = open(self.temp_posterior_path, 'rU')
                out = GzipFileStream(self.posterior_out_path, 'w',
                        compresslevel = 9)
                for line in post_stream:
                    out.write(line)
                out.close()
                post_stream.close()
                if self.posterior_out_path == self.posterior_path + '.gz':
                    os.remove(self.posterior_path)
            except:
                _LOG.warning('{0}: problem compressing posterior sample... '
                        'returning uncompressed'.format(self.name))
                base, ext = os.path.splitext(self.posterior_out_path)
                if ext.lower() == '.gz':
                    self.posterior_out_path = base
                self.posterior_out_path = self.posterior_out_path
                shutil.move(self.temp_posterior_path, self.posterior_out_path)
                self.temp_posterior_path = self.posterior_out_path
        else:
            shutil.move(self.temp_posterior_path, self.posterior_out_path)
            self.temp_posterior_path = self.posterior_out_path

    def _prep_regression_worker(self):
        header = parsing.parse_header(self.temp_posterior_path)
        all_patterns = set(parsing.MEAN_TAU_PATTERNS + parsing.OMEGA_PATTERNS +
                parsing.CV_PATTERNS + parsing.MODEL_PATTERNS +
                parsing.PSI_PATTERNS + parsing.DIV_MODEL_PATTERNS)
        # Adding tau parameters causes ABCEstimator to fail often
        # if MSBAYES_SORT_INDEX.current_value() == 0:
        #     all_patterns.update(parsing.TAU_PATTERNS)
        parameter_patterns = sorted(list(all_patterns - \
                self.parameter_patterns_to_remove))
        parameter_indices = sorted(parsing.get_parameter_indices(header,
                parameter_patterns = parameter_patterns))
        self.regression_worker = ABCToolBoxRegressWorker(
                temp_fs = self.temp_fs,
                observed_path = self.observed_path,
                posterior_path = self.temp_posterior_path,
                parameter_indices = parameter_indices,
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
        if self.regression_worker.failed:
            self.regression_failed = True
        discrete_parameter_patterns = list(set(parsing.PSI_PATTERNS +
                parsing.MODEL_PATTERNS +
                parsing.DIV_MODEL_PATTERNS) - self.parameter_patterns_to_remove)
        try:
            discrete_probs = stats.summarize_discrete_parameters_from_densities(
                    self.regress_posterior_path,
                    discrete_parameter_patterns = discrete_parameter_patterns,
                    include_omega_summary = True,
                    omega_threshold = self.omega_threshold,
                    include_cv_summary = True,
                    cv_threshold = self.cv_threshold)
        except:
            _LOG.warning('{0}: Problem parsing posterior density file {1!r}'
                    '... Returning unadjusted results'.format(
                            self.name, self.regress_posterior_path))
            discrete_probs = {}
        if max(discrete_probs.get('PRI.Psi', {-1: None}).iterkeys()) > \
                self.num_taxon_pairs:
            raise ValueError('number of taxon pairs is {0}, but found '
                    'psi estimates of {1} in posterior {2}'.format(
                        self.num_taxon_pairs,
                        max(discrete_probs['PRI.Psi'].iterkeys()),
                        self.regress_posterior_path))
        for i in range(self.num_taxon_pairs):
            self.adjusted_psi_probs[i+1] = discrete_probs.get('PRI.Psi', 
                    {i+1: float('nan')}).get(i+1, 0.0)
        for k, s in self.div_model_summary:
            s['adjusted_frequency'] = discrete_probs.get('PRI.div.model', 
                    {self.top_div_models_to_indices[k]: float('nan')}).get(
                            self.top_div_models_to_indices[k], 0.0)
        self.adjusted_prob_omega_zero = discrete_probs.get(
                'PRI.omega', {0: float('nan')}).get(0, 0.0)
        self.adjusted_prob_cv_zero = discrete_probs.get(
                'PRI.cv', {0: float('nan')}).get(0, 0.0)
        if self.model_probs:
            self.adjusted_model_probs = {}
            for i in self.model_probs.iterkeys():
                self.adjusted_model_probs[i] = discrete_probs.get(
                        'PRI.model', {i: float('nan')}).get(i, 0.0)

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

    def _write_cv_results(self):
        if self.cv_included:
            out, close = process_file_arg(self.cv_results_path, 'w')
            out.write('cv_thresh\tprob_less_than\t'
                    'glm_prob_less_than\n')
            out.write('{0}\t{1}\t{2}\n'.format(
                    self.cv_threshold,
                    self.prob_cv_zero,
                    self.adjusted_prob_cv_zero))
            if close:
                out.close()

    def _write_summary(self):
        try:
            glm = parsing.parse_abctoolbox_summary_file(
                    self.regress_summary_path)
        except:
            _LOG.warning('{0}: Parsing of regression summary file {1!r} '
                    'failed... Returning unadjusted summary'.format(
                        self.name, self.regress_summary_path))
            glm = {}
        out, close = process_file_arg(self.posterior_summary_path, 'w')
        for param, summary in self.unadjusted_summaries.iteritems():
            out.write('[{0}]\n'.format(param))
            if isinstance(summary['modes'][0], tuple):
                out.write('    modes = {0}\n'.format(', '.join([
                    '({0}, {1})'.format(x, y) for x, y in summary['modes']])))
            else:
                out.write('    modes = {0}\n'.format(', '.join(
                    [str(x) for x in summary['modes']])))
            out.write('    mode_glm = {0}\n'.format(
                    glm.get(param, {'mode': float('nan')})['mode']))
            out.write('    median = {0}\n'.format(summary['median']))
            out.write('    median_glm = {0}\n'.format(
                    glm.get(param, {'median': float('nan')})['median']))
            out.write('    mean = {0}\n'.format(summary['mean']))
            out.write('    mean_glm = {0}\n'.format(
                    glm.get(param, {'mean': float('nan')})['mean']))
            out.write('    n = {0}\n'.format(summary['n']))
            out.write('    range = {0}, {1}\n'.format(summary['range'][0],
                    summary['range'][1]))
            out.write('    HPD_95_interval = {0}, {1}\n'.format(
                    summary['hpdi_95'][0],
                    summary['hpdi_95'][1]))
            out.write('    HPD_95_interval_glm = {0}, {1}\n'.format(
                    glm.get(param, {'HPD_95_lower_bound': float('nan')})[
                            'HPD_95_lower_bound'],
                    glm.get(param, {'HPD_95_upper_bound': float('nan')})[
                            'HPD_95_upper_bound']))
            out.write('    quantile_95_interval = {0}, {1}\n'.format(
                    summary['qi_95'][0],
                    summary['qi_95'][1]))
            out.write('    quantile_95_interval_glm = {0}, {1}\n'.format(
                    glm.get(param, {'quantile_95_lower_bound': float('nan')})[
                            'quantile_95_lower_bound'],
                    glm.get(param, {'quantile_95_upper_bound': float('nan')})[
                            'quantile_95_upper_bound']))
        if close:
            out.close()

    def _post_process(self):
        if not self.keep_temps:
            self.temp_fs.remove_dir(self.temp_output_dir)

    def start(self):
        self._pre_process()
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
        self._write_cv_results()
        self._write_summary()
        self._post_process()
        self.finished = True

class ModelProbabilityEstimator(object):
    count = 0
    def __init__(self, config,
            num_samples = 1000,
            omega_threshold = 0.01,
            cv_threshold = 0.01,
            rng = None,
            tag = None):
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.config = config
        if not isinstance(config, MsBayesConfig):
            self.config = MsBayesConfig(config)
        self.npairs = self.config.npairs
        self.num_samples = num_samples
        self.omega_threshold = omega_threshold
        self.cv_threshold = cv_threshold
        self.shared_div_summary = dict(zip(
                [i + 1 for i in range(1, self.npairs)],
                [stats.SampleSummary() for i in range(1, self.npairs)]))
        self.psi_summary = dict(zip([i + 1 for i in range(self.npairs)],
                [stats.SampleSummary() for i in range(self.npairs)]))
        self.omega_summary = stats.SampleSummary()
        self.cv_summary = stats.SampleSummary()
        self.rng = rng
        self.tag = tag
        self.finished = False

    def get_prior_sample_iter(self):
        p = stats.Partition()
        if self.config.div_model_prior == 'psi':
            return p.psi_multinomial_draw_iter(
                    num_samples = self.num_samples,
                    num_elements = self.config.npairs,
                    base_distribution = self.config.tau,
                    rng = self.rng)
        elif self.config.div_model_prior == 'uniform':
            return p.uniform_integer_partition_draw_iter(
                    num_samples = self.num_samples,
                    num_elements = self.config.npairs,
                    base_distribution = self.config.tau,
                    rng = self.rng)
        elif self.config.div_model_prior == 'dpp':
            return p.dirichlet_process_draw_iter(
                    alpha = self.config.dpp_concentration,
                    num_samples = self.num_samples,
                    num_elements = self.config.npairs,
                    base_distribution = self.config.tau,
                    rng = self.rng)
        else:
            raise Exception('divergence model prior {0!r} is not '
                    'supported'.format(self.config.div_model_prior))

    def start(self):
        codiv_summarizer = dict(zip([i + 1 for i in range(1, self.npairs)],
                [stats.SampleSummarizer() for i in range(1, self.npairs)]))
        psi_summarizer = dict(zip([i + 1 for i in range(self.npairs)],
                [stats.SampleSummarizer() for i in range(self.npairs)]))
        omega_summarizer = stats.SampleSummarizer()
        cv_summarizer = stats.SampleSummarizer()
        prior_sample_iter = self.get_prior_sample_iter()
        rng = self.rng
        if not rng:
            rng = GLOBAL_RNG
        for partition in prior_sample_iter:
            for k in psi_summarizer.iterkeys():
                if len(partition.values) == k:
                    psi_summarizer[k].add_sample(1)
                else:
                    psi_summarizer[k].add_sample(0)
            for k in codiv_summarizer.iterkeys():
                s = rng.sample(partition.partition, k)
                if len(set(s)) == 1:
                    codiv_summarizer[k].add_sample(1)
                else:
                    codiv_summarizer[k].add_sample(0)
            ss = stats.SampleSummarizer(
                    samples = partition.get_element_vector())
            d = ss.variance / ss.mean
            if d < self.omega_threshold:
                omega_summarizer.add_sample(1)
            else:
                omega_summarizer.add_sample(0)
            cv = ss.std_deviation / ss.mean
            if cv < self.cv_threshold:
                cv_summarizer.add_sample(1)
            else:
                cv_summarizer.add_sample(0)
        for k in self.psi_summary.iterkeys():
            self.psi_summary[k].update(stats.SampleSummary(
                sample_size = psi_summarizer[k].n,
                mean = psi_summarizer[k].mean,
                variance = psi_summarizer[k].variance))
        for k in self.shared_div_summary.iterkeys():
            self.shared_div_summary[k].update(stats.SampleSummary(
                sample_size = codiv_summarizer[k].n,
                mean = codiv_summarizer[k].mean,
                variance = codiv_summarizer[k].variance))
        self.omega_summary.update(stats.SampleSummary(
                sample_size = omega_summarizer.n,
                mean = omega_summarizer.mean,
                variance = omega_summarizer.variance))
        self.cv_summary.update(stats.SampleSummary(
                sample_size = cv_summarizer.n,
                mean = cv_summarizer.mean,
                variance = cv_summarizer.variance))
        self.finished = True

class DivModelSimulator(object):
    count = 0
    def __init__(self, config,
            num_samples = 1000,
            rng = None,
            tag = None):
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.config = config
        if not isinstance(config, MsBayesConfig):
            self.config = MsBayesConfig(config)
        self.npairs = self.config.npairs
        self.num_samples = num_samples
        self.div_models = stats.PartitionCollection()
        self.rng = rng
        self.tag = tag
        self.finished = False

    def get_prior_sample_iter(self):
        p = stats.Partition()
        if self.config.div_model_prior == 'psi':
            return p.psi_multinomial_draw_iter(
                    num_samples = self.num_samples,
                    num_elements = self.config.npairs,
                    base_distribution = self.config.tau,
                    rng = self.rng)
        elif self.config.div_model_prior == 'uniform':
            return p.uniform_integer_partition_draw_iter(
                    num_samples = self.num_samples,
                    num_elements = self.config.npairs,
                    base_distribution = self.config.tau,
                    rng = self.rng)
        elif self.config.div_model_prior == 'dpp':
            return p.dirichlet_process_draw_iter(
                    alpha = self.config.dpp_concentration,
                    num_samples = self.num_samples,
                    num_elements = self.config.npairs,
                    base_distribution = self.config.tau,
                    rng = self.rng)
        else:
            raise Exception('divergence model prior {0!r} is not '
                    'supported'.format(self.config.div_model_prior))

    def start(self):
        self.div_models.add_iter(self.get_prior_sample_iter())
        self.finished = True

class DppSimWorker(object):
    count = 0
    def __init__(self,
            alpha,
            num_elements,
            num_samples = 1000,
            base_distribution = None,
            rng = None,
            tag = None):
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.alpha = alpha
        self.num_elements = num_elements
        self.base_distribution = base_distribution
        self.num_samples = num_samples
        self.div_models = stats.PartitionCollection()
        self.rng = rng
        self.tag = tag
        self.psi_summary = dict(zip([i + 1 for i in range(self.num_elements)],
                [stats.SampleSummary() for i in range(self.num_elements)]))
        self.finished = False

    def get_prior_sample_iter(self):
        p = stats.Partition()
        return p.dirichlet_process_draw_iter(
                alpha = self.alpha,
                num_samples = self.num_samples,
                num_elements = self.num_elements,
                base_distribution = self.base_distribution,
                rng = self.rng)

    def start(self):
        psi_summarizer = dict(zip([i + 1 for i in range(self.num_elements)],
                [stats.SampleSummarizer() for i in range(self.num_elements)]))
        rng = self.rng
        if not rng:
            rng = GLOBAL_RNG
        for partition in self.get_prior_sample_iter():
            self.div_models.add(partition)
            for k in psi_summarizer.iterkeys():
                if len(partition.values) == k:
                    psi_summarizer[k].add_sample(1)
                else:
                    psi_summarizer[k].add_sample(0)
        for k in self.psi_summary.iterkeys():
            self.psi_summary[k].update(stats.SampleSummary(
                sample_size = psi_summarizer[k].n,
                mean = psi_summarizer[k].mean,
                variance = psi_summarizer[k].variance))
        self.finished = True


class GeneTreeSimWorker(GenericWorker):
    count = 0
    def __init__(self,
            args,
            stdout_path = None,
            stderr_path = None,
            tag = None):
        GenericWorker.__init__(self,
                exe = 'genetree_in_sptree_sim',
                args = args,
                stdout_path = stdout_path,
                stderr_path = stderr_path,
                tag = tag)
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)


class SequenceSimWorker(GenericWorker):
    count = 0
    def __init__(self,
            args,
            output_path,
            stdout_path = None,
            stderr_path = None,
            tag = None):
        GenericWorker.__init__(self,
                exe = 'simulate_sequences',
                args = args,
                stdout_path = stdout_path,
                stderr_path = stderr_path,
                tag = tag)
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.output_path = functions.get_new_path(output_path)
        self.tmp_dir = None

    def _pre_process(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.subprocess_kwargs['cwd'] = self.tmp_dir

    def _post_process(self):
        f = os.listdir(self.tmp_dir)
        if not ((len(f) == 1) and (os.path.splitext(f[0])[-1] == '.nex')):
            raise Exception('{0} produced unexpected output'.format(self.name))
        p = os.path.join(self.tmp_dir, f[0])
        shutil.move(p, self.output_path)
        os.rmdir(self.tmp_dir)
