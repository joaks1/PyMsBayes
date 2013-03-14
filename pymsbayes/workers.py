#! /usr/bin/env python

import os
import sys
import multiprocessing
import time

from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__, 'debug')
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
            exe_path,
            config_path,
            out_path,
            seed,
            report_parameters = False,
            **kwargs):
        Worker.__init__(self, **kwargs)
        self.__class__.count += 1
        self.name = 'MsBayesWorker-' + str(self.count)
        self.exe_path = exe_path
        self.config_path = config_path
        self.out_path = out_path
        self.seed = seed
        self.kill_received = False

    def run(self):
        for i in range(4):
            self.send_warning('working...')
            time.sleep(2)
        

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

