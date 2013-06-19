#! /usr/bin/env python

import os
import sys
import multiprocessing
import Queue

from pymsbayes.utils import WORK_FORCE
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)
_LOCK = multiprocessing.Lock()

class Manager(multiprocessing.Process):
    count = 0
    def __init__(self,
            work_queue = None,
            result_queue = None,
            log = None,
            lock = None):
        multiprocessing.Process.__init__(self)
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        if not work_queue:
            work_queue = WORK_FORCE
        if not result_queue:
            result_queue = multiprocessing.Queue()
        self.work_queue = work_queue
        self.result_queue = result_queue
        if not log:
            log = _LOG
        self.log = log
        if not lock:
            lock = _LOCK
        self.lock = lock
        self.killed = False

    def compose_msg(self, msg):
        return '{0} ({1}): {2}'.format(self.name, self.pid, msg)

    def send_msg(self, msg, method_str='info'):
        self.lock.acquire()
        try:
            getattr(self.log, method_str)(self.compose_msg(msg))
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
        while not self.killed:
            try:
                worker = self.work_queue.get(block=True, timeout=0.1)
                # without blocking processes were stopping when the queue
                # was not empty, and without timeout, the processes would
                # hang waiting for jobs.
            except Queue.Empty:
                break
            self.send_info('starting worker {0}'.format(worker.name))
            worker.start()
            self.send_info('worker {0} finished'.format(worker.name))
            self.result_queue.put(worker)
        if self.killed:
            self.send_error('manager was killed!')

