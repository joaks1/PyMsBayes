#! /usr/bin/env python

import os
import sys
import re
import multiprocessing
import subprocess
import time

from pymsbayes.utils import BIN_DIR
from pymsbayes.utils.functions import (get_random_int, expand_path,
        get_indices_of_patterns, reduce_columns, process_file_arg)
from pymsbayes.utils.errors import WorkerExecutionError
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)
_LOCK = multiprocessing.Lock()


##############################################################################
## msBayes prior header patterns

PARAMETER_PATTERNS = [
        re.compile(r'\s*PRI\.(?!numTauClass)\S+\s*$'),
        ]
DEFAULT_STAT_PATTERNS = [
        re.compile(r'\s*pi\.\d+\s*'),
        re.compile(r'\s*wattTheta\.\d+\s*'),
        re.compile(r'\s*pi\.net\.\d+\s*'),
        re.compile(r'\s*tajD\.denom\.\d+\s*'),
        ]
DUMMY_PATTERNS = [
        re.compile(r'\s*PRI\.numTauClass\s*')
        ]

##############################################################################
## functions for manipulating prior files

def parse_header(file_obj, sep='\t'):
    file_stream, close = process_file_arg(in_file, 'rU')
    header = file_stream.next().strip().split(sep)
    if close:
        file_stream.close()
    else:
        file_stream.seek(0)
    return header

def get_parameters_indices(header_list, parameter_patterns=PARAMETER_PATTERNS):
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

def observed_parameters_for_abctoolbox(in_file, out_file,
        parameter_patterns=PARAMETER_PATTERNS):
    header = parse_header(in_file)
    indices = get_parameter_indices(header,
            parameter_patterns=parameter_patterns)
    reduce_columns(in_file, out_file, indices)

def observed_for_acceptrej(in_file, out_file,
        stat_patterns=DEFAULT_STAT_PATTERNS,
        parameter_patterns=PARAMETER_PATTERNS):
    # in_file, close_in = process_file_arg(in_file, 'rU')
    # out_file, close_out = process_file_arg(out_file, 'w')
    header = parse_header(in_file)
    parameter_indices = get_parameter_indices(header,
            parameter_patterns=parameter_patterns)
    stat_indices = get_stat_indices(header,
            stat_patterns=stat_patterns)
    indices = sorted(parameter_indices + stat_indices)
    reduce_columns(in_file, out_file, indices)
    # new_head = [header[i] for i in (parameter_indices + stat_indices)]
    # out_file.write('%s\t\n' % '\t'.join(new_head))
    # line_iter = iter(in_file)
    # line_iter.next()
    # for line_num, line in enumerate(line_iter):
    #     l = line.strip().split()
    #     new_line = ['1', '0', '1', '0'] + [l[i] for i in stat_indices]
    #     out_file.write('%s\t\n' % '\t'.join(new_line))
    # if close_in:
    #     in_file.close()
    # if close_out:
    #     out_file.close()

# def observed_for_msreject(in_file, out_file,
#         stat_patterns=DEFAULT_STAT_PATTERNS,
#         parameter_patterns=PARAMETER_PATTERNS,
#         dummy_patterns=DUMMY_PATTERNS):
#     # in_file, close_in = process_file_arg(in_file, 'rU')
#     # out_file, close_out = process_file_arg(out_file, 'w')
#     header = parse_header(in_file)
#     parameter_indices = get_parameter_indices(header,
#             parameter_patterns=parameter_patterns)
#     stat_indices = get_stat_indices(header,
#             stat_patterns=stat_patterns)
#     dummy_indices = get_dummy_indices(header,
#             dummy_patterns=DUMMY_PATTERNS)
#     indices = sorted(dummy_indices + parameter_indices + stat_indices)
#     reduce_columns(in_file, out_file, indices)
#     # new_head = [header[i] for i in (dummy_indices + parameter_indices + stat_indices)]
#     # out_file.write('%s\t\n' % '\t'.join(new_head))
#     # line_iter = iter(in_file)
#     # line_iter.next()
#     # for line_num, line in enumerate(line_iter):
#     #     l = line.strip().split()
#     #     new_line = ['0', '1', '0', '1', '0'] + [l[i] for i in stat_indices]
#     #     out_file.write('%s\n' % '\t'.join(new_line))

def prior_for_abctoolbox(in_file, out_file,
        stat_patterns=DEFAULT_STAT_PATTERNS,
        parameter_patterns=PARAMETER_PATTERNS):
    header = parse_header(in_file)
    indices = get_parameter_indices(header,
            parameter_patterns=parameter_patterns)
    indices.extend(get_stat_indices(header, stat_patterns=stat_patterns))
    reduce_columns(in_file, out_file, sorted(indices))
    
def prior_for_msreject(in_file, out_file,
        stat_patterns=DEFAULT_STAT_PATTERNS,
        parameter_patterns=PARAMETER_PATTERNS,
        dummy_patterns=DUMMY_PATTERNS,
        header=False):
    in_file, close = process_file_arg(in_file)
    header = parse_header(in_file)
    indices = get_parameter_indices(header,
            parameter_patterns=parameter_patterns)
    indices.extend(get_stat_indices(header, stat_patterns=stat_patterns))
    indices.extend(get_dummy_indices(header, dummy_patterns=DUMMY_PATTERNS))
    if not header:
        in_file.next()
    reduce_columns(in_file, out_file, sorted(indices), extra_tab=True)
    if close:
        in_file.close()


##############################################################################
## Base class for all workers

class Worker(multiprocessing.Process):
    total = 0
    def __init__(self, **kwargs):
        self.__class__.total += 1
        multiprocessing.Process.__init__(self)
        self.log = kwargs.get('log', _LOG)
        self.lock = kwargs.get('lock', _LOCK)
        self.queue = kwargs.get('queue', multiprocessing.Queue())
        self.stdout_path = kwargs.get('stdout_path', None)
        self.stderr_path = kwargs.get('stderr_path', None)
        self.cmd = []
        self.finished = False
        self.subprocess_exit_code = None

    def send_msg(self, msg, method_str='info'):
        msg = '{0} ({1}): {2}'.format(self.name, self.pid, msg)
        self.lock.acquire()
        try:
            getattr(self.log, method_str)(msg)
        finally:
            self.lock.release()

    def send_debug(self, msg):
        self.send_msg(msg, method_str='debug')

    def send_info(self, msg):
        self.send_msg(msg, method_str='info')

    def send_warning(self, msg):
        self.send_msg(msg, method_str='warning')

    def send_error(self, msg):
        self.send_msg(msg, method_str='error')

    def get_stderr(self):
        if not self.stderr_path:
            return None
        try:
            return open(self.stderr_path, 'rU').read()
        except IOError, e:
            self.send_error('Could not open stderr file')
            raise e

    def run(self):
        self.send_info('Starting process with following command:\n\t'
                '{0}'.format(' '.join(self.cmd)))
        if self.stdout_path:
            sout = open(self.stdout_path, 'w')
        else:
            sout = subprocess.PIPE
        if self.stderr_path:
            serr = open(self.stderr_path, 'w')
        else:
            serr = subprocess.PIPE
        p = subprocess.Popen(self.cmd,
                stdout = sout,
                stderr = serr,
                shell = False)
        stdout, stderr = p.communicate()
        exit_code = p.wait()
        if self.stdout_path:
            sout.close()
        if self.stderr_path:
            serr.close()
        if exit_code != 0:
            if self.stderr_path:
                stderr = open(self.stderr_path, 'rU').read()
            send_error('execution failed')
            raise WorkerExecutionError('{0} ({1}) failed. stderr:\n{2}'.format(
                self.name, self.pid, stderr))
        results = {'exit_code': exit_code}
        self.queue.put(results)

    def finish(self, **kwargs):
        results = self.queue.get()
        self.subprocess_exit_code = results['exit_code']
        try:
            self._post_process(**kwargs)
        except Exception, e:
            self.send_error('Error during post-processing')
            raise e
        self.finished = True

    def _post_process(self, **kwargs):
        pass


##############################################################################
## msBayes class for generating prior files

class MsBayesWorker(Worker):
    count = 0
    def __init__(self,
            sample_size,
            config_path,
            output_dir,
            output_prefix = 'prior',
            exe_path = None,
            model_index = None,
            sort_index = None,
            report_parameters = True,
            seed = None,
            **kwargs):
        Worker.__init__(self, **kwargs)
        self.__class__.count += 1
        self.name = 'MsBayesWorker-' + str(self.count)
        self.sample_size = int(sample_size)
        self.config_path = expand_path(config_path)
        self.output_dir = expand_path(output_dir)
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
        self.out_path = os.path.join(
                self.output_dir,
                '{0}-{1}-{2}.txt'.format(
                        output_prefix,
                        self.sample_size,
                        self.seed))
        self._update_cmd()

    def _update_cmd(self):
        cmd = [self.exe_path,
               '-r', str(self.sample_size),
               '-c', self.config_path,
               '-o', self.out_path,
               '-S', str(self.seed),]
        if self.sort_index != None:
            cmd.extend(['-s', str(self.sort_index)])
        if self.model_index != None:
            cmd.extend(['-m', str(self.model_index)])
        if self.report_parameters:
            cmd.append('-p')
        self.cmd = cmd
        
class MsRejectWorker(Worker):
    count = 0
    def __init__(self,
            exe_path,
            prior_path,
            out_path,
            tolerance,
            stats,
            **kwargs):
        Worker.__init__(self, **kwargs)
        self.__class__.count += 1
        self.name = 'MsRejectWorker-' + str(self.count)
        self.exe_path = e
        self.prior_path = prior_path
        self.out_path = out_path
        self.tolerance = tolerance
        self.stats = stats
        self.kill_received = False

    def run(self):
        time.sleep(120)

if __name__ == '__main__':
    jobs = []
    for i in range(5):
        p = MsBayesWorker(
                exe_path = 'msbayes.pl',
                config_path = 'conf',
                out_path = 'prior',
                seed = 234234)
        jobs.append(p)
        p.start()
    for j in jobs:
        j.join()

