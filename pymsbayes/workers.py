#! /usr/bin/env python

import os
import sys
import multiprocessing
import subprocess
import time

from pymsbayes.utils import BIN_DIR, get_random_int
from pymsbayes.utils.errors import MsBayesExecutionError
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

    def run(self):
        pass

class MsBayesWorker(Worker):
    count = 0
    def __init__(self,
            sample_size,
            config_path,
            out_path,
            exe_path = None,
            model_index = None,
            sort_index = None,
            report_parameters = False,
            seed = None,
            **kwargs):
        Worker.__init__(self, **kwargs)
        self.__class__.count += 1
        self.name = 'MsBayesWorker-' + str(self.count)
        self.sample_size = int(sample_size)
        self.config_path = config_path
        self.out_path = out_path
        self.exe_path = exe_path
        if not self.exe_path:
            self.exe_path = os.path.join(BIN_DIR, 'msbayes.pl')
        self.model_index = int(model_index)
        self.sort_index = int(sort_index)
        self.report_parameters = report_parameters
        self.seed = int(seed)
        if not self.seed:
            self.seed = get_random_int()
        self.cmd = []
        self.update_cmd()
        self.kill_received = False
        self.finished = False
        self.msbayes_stdout = None
        self.msbayes_stderr = None
        self.msbayes_exit_code = None

    def update_cmd(self):
        cmd = [self.exe_path,
               '-r', str(self.sample_size),
               '-c', self.config_path,
               '-o', self.out_path,
               '-S', str(self.seed),]
        if self.sort_index:
            cmd.extend(['-s', str(self.sort_index)])
        if self.model_index:
            cmd.extend(['-m', str(self.model_index)])
        if self.report_parameters:
            cmd.append('-p')
        self.cmd = cmd

    def run(self):
        self.update_cmd()
        self.send_info('Starting msbayes.pl with following command:\n\t'
                '{0}'.format(' '.join(self.cmd)))
        p = subprocess.Popen(self.cmd,
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE,
                shell = False)
        self.msbayes_stdout, self.msbayes_stderr = p.communicate()
        self.msbayes_exit_code = p.wait()
        

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

