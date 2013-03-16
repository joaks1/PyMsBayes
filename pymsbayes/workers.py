#! /usr/bin/env python

import os
import sys
import multiprocessing
import subprocess
import time

from pymsbayes.utils import BIN_DIR, get_random_int, expand_path
from pymsbayes.utils.errors import WorkerExecutionError
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)
_LOCK = multiprocessing.Lock()

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

    def finish(self):
        results = self.queue.get()
        self.subprocess_exit_code = results['exit_code']
        self.finished = True

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

