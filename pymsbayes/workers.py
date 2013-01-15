#! /usr/bin/env python

import os
import sys
import multiprocessing
import time

from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__, 'debug')
_LOCK = multiprocessing.Lock()

class Worker(multiprocessing.Process):
    def __init__(self,
            process_id,
            log = None,
            lock = None):
        multiprocessing.Process.__init__(self)
        self.process_id = process_id
        self.log = log
        if not log:
            self.log = _LOG
        self.lock = lock
        if not lock:
            self.lock = _LOCK

    def send_msg(self, msg, method_str='info'):
        msg = '{0}: {1}: {2}'.format(self.name, self.process_id, msg)
        self.lock.acquire()
        try:
            getattr(self.log, method_str)(msg)
        finally:
            self.lock.release()

    def send_debug(self, msg):
        self.send_msg(msg, method_str='debug')

    def send_info(self, msg):
        self.send_msg(msg, method_str='info')

    def send_warn(self, msg):
        self.send_msg(msg, method_str='warn')

    def run(self):
        pass

class MsBayesWorker(Worker):
    def __init__(self,
            process_id,
            exe_path,
            config_path,
            out_path,
            seed,
            log = None,
            lock = None,
            report_parameters = False):
        Worker.__init__(self,
                process_id = process_id,
                log = log,
                lock = lock)
        self.exe_path = exe_path
        self.config_path = config_path
        self.out_path = out_path
        self.seed = seed
        self.kill_received = False

    def run(self):
        for i in range(4):
            self.send_warn('working...')
            time.sleep(30)
        

class MsRejectWorker(Worker):
    def __init__(self,
            process_id,
            exe_path,
            prior_path,
            out_path,
            tolerance,
            stats,
            log = None,
            lock = None):
        Worker.__init__(self,
                process_id = process_id,
                log = log,
                lock = lock)
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
                process_id = 'msbw-{0}'.format(i),
                exe_path = 'msbayes.pl',
                config_path = 'conf',
                out_path = 'prior',
                seed = 234234)
        jobs.append(p)
        p.start()
    for j in jobs:
        j.join()

