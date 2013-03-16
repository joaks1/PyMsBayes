#! /usr/bin/env python

import unittest
import os

from pymsbayes import workers
from pymsbayes.test.support import package_paths
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class MsBayesWorkerTestCase(PyMsBayesTestCase):
    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path('4pairs_1locus.cfg')
        self.prior_dir = self.get_test_subdir(prefix='prior_files')

    def tearDown(self):
        self.tear_down()

    def test_simple(self):
        w = workers.MsBayesWorker(
                sample_size = 10,
                config_path = self.cfg_path,
                output_dir = self.prior_dir,
                output_prefix = self.test_id +'-test_prior')
        self.assertIsInstance(w, workers.MsBayesWorker)
        self.assertFalse(w.finished)
        w.start()
        w.join()
        w.finish()
        self.assertFalse(w.is_alive())
        _LOG.debug('{0}\n{1}\n{2}'.format(w.finished,
            w.subprocess_exit_code,
            w.exitcode))


if __name__ == '__main__':
    unittest.main()

