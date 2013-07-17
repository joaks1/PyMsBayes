#! /usr/bin/env python

import os
import sys
import multiprocessing
import Queue
import time

from pymsbayes.utils import WORK_FORCE
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)
_LOCK = multiprocessing.Lock()

class Manager(multiprocessing.Process):
    count = 0
    def __init__(self,
            work_queue = None,
            result_queue = None,
            get_timeout = 0.4,
            put_timeout = 0.2,
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
        self.get_timeout = get_timeout
        self.put_timeout = put_timeout
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

    def _get_worker(self):
        worker = None
        try:
            self.send_debug('getting worker')
            worker = self.work_queue.get(block=True, timeout=self.get_timeout)
            self.send_debug('received worker {0}'.format(
                    getattr(worker, 'name', 'nameless')))
            # without blocking processes were stopping when the queue
            # was not empty, and without timeout, the processes would
            # hang waiting for jobs.
        except Queue.Empty:
            time.sleep(0.2)
            if not self.work_queue.empty():
                self.send_warning('raised Queue.Empty, but queue is '
                        'not empty... trying again')
                return self._get_worker()
            else:
                self.send_info('work queue is empty')
        return worker

    def _put_worker(self, worker):
        try:
            self.send_debug('returning worker {0}'.format(
                    getattr(worker, 'name', 'nameless')))
            self.result_queue.put(worker, block=True, timeout=self.put_timeout)
            self.send_debug('worker {0} returned'.format(
                    getattr(worker, 'name', 'nameless')))
        except Queue.Full, e:
            time.sleep(0.2)
            if not self.result_queue.full():
                self.send_warning('raised Queue.Full, but queue is '
                        'not full... trying again')
                self._put_worker(worker)
            else:
                self.send_error('result queue is full... aborting')
            self.killed = True
            raise e

    def run(self):
        self.send_debug('starting run')
        while not self.killed:
            worker = self._get_worker()
            if worker is None:
                break
            self.send_info('starting worker {0}'.format(
                    getattr(worker, 'name', 'nameless')))
            worker.start()
            self.send_info('worker {0} finished'.format(
                    getattr(worker, 'name', 'nameless')))
            self._put_worker(worker)
        if self.killed:
            self.send_error('manager was killed!')
        self.send_debug('end run')

