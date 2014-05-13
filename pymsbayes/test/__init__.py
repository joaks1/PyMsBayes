#!/usr/bin/env python

import unittest
import re
import os
import sys
import logging
import glob

from pymsbayes.utils.oproperty import OProperty as property
from pymsbayes.test.support import package_paths
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class TestLevel:
    FAST, NORMAL, SLOW, EXHAUSTIVE = 0, 10, 20, 30
    def name(i):
        if i <= TestLevel.FAST:
            return "FAST"
        if i <= TestLevel.NORMAL:
            return "NORMAL"
        if i <= TestLevel.SLOW:
            return "SLOW"
        return "EXHAUSTIVE"

    name = staticmethod(name)

    @classmethod
    def name_to_int(cls, l):
        try:
            return int(l)
        except:
            pass
        l = l.upper()
        if l == "FAST":
            return cls.FAST
        if l == "NORMAL":
            return cls.NORMAL
        if l == "SLOW":
            return cls.SLOW
        if l == "EXHAUSTIVE":
            return cls.EXHAUSTIVE
        raise ValueError('TestLevel {0} unrecognized'.format(l))

    @classmethod
    def test_skip_msg(cls, log, module_name, level, message=None):
        if message is None:
            message = "tests skipped"
        log.warning(
            'Running in {0} Testing Level. Skipping {1} tests in {2}: '
                    '{3}'.format(
                        cls.name(cls.get_current_level()),
                        cls.name(level),
                        module_name,
                        message))
    
    @classmethod
    def get_current_level(cls):
        l = os.environ.get("PYMSBAYES_TESTING_LEVEL")
        if l is None:
            return cls.FAST
        try:
            return cls.name_to_int(l)
        except:
            _LOG.warn('the value {0} for PYMSBAYES_TESTING_LEVEL is not recognized.'
                      'Using NORMAL level'.format(l))
        return cls.FAST
    
    @classmethod
    def test_enabled(cls, level, log=None, module_name="", message=None):
        tl = cls.get_current_level()
        if level > tl:
            if log:
                cls.test_skip_msg(log, module_name, level, message)
            return False
        return True

def get_test_file_names():
    if TestLevel.get_current_level() < TestLevel.NORMAL:
        test_file_paths = glob.glob(os.path.join(package_paths.TEST_DIR,
                "test_end_user*.py"))
    else:
        test_file_paths = glob.glob(os.path.join(package_paths.TEST_DIR,
                "test_*.py"))
    tests = []
    for path in test_file_paths:
        test_file_name = os.path.basename(path).split(".")[0]
        test = ".".join(["pymsbayes.test", test_file_name])
        sys.stdout.write("%s\n" % test)
        tests.append(test)
    return tests

def get_unittest_suite():
    test_file_names = get_test_file_names()
    tests = unittest.defaultTestLoader.loadTestsFromNames(test_file_names)
    suite = unittest.TestSuite(tests)
    return suite

def run():
    test_runner = unittest.TextTestRunner()
    test_runner.run(get_unittest_suite())

def load_tests(loader, tests, pattern):
    # override default unittest test loader
    return get_unittest_suite()

if __name__ == '__main__':
    run()

