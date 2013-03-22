#! /usr/bin/env python

import os
import sys
import multiprocessing
import Queue

from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)
_LOCK = multiprocessing.Lock()

class Manager(multiprocessing.Process):
    count = 0
    def __init__(self, **kwargs):
        multiprocessing.Process.__init__(self)
        self.__class__.count += 1
        self.name = 'Manager-' + str(self.count)
        self.log = kwargs.get('log', _LOG)
        self.lock = kwargs.get('lock', _LOCK)
        self.work_queue = kwargs.get('job_queue', multiprocessing.Queue())
        self.result_queue = kwargs.get('result_queue', multiprocessing.Queue())
        self.killed = False

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
        while not self.killed:
            try:
                worker = self.work_queue.get_nowait()
            except Queue.Empty:
                break
            self.send_info('starting worker {0}'.format(worker.name))
            worker.start()
            self.send_info('worker {0} finished'.format(worker.name))
            self.result_queue.put(worker)
        if self.killed:
            self.send_error('manager was killed!')

