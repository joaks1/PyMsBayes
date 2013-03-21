#!/usr/bin/env python

import unittest
import re
import os
import sys
import logging
import glob

from pymsbayes.test.support import package_paths
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

def get_test_file_names():
    test_file_paths = glob.glob(os.path.join(package_paths.TEST_DIR, "test_*"))
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

    def name_to_int(l):
        try:
            return int(l)
        except:
            pass
        l = l.upper()
        if l == "FAST":
            return TestLevel.FAST
        if l == "NORMAL":
            return TestLevel.NORMAL
        if l == "SLOW":
            return TestLevel.SLOW
        if l == "EXHAUSTIVE":
            return TestLevel.EXHAUSTIVE
        raise ValueError('TestLevel {0} unrecognized'.format(l))

    name_to_int = staticmethod(name_to_int)

def test_skip_msg(log, module_name, message=None, level=TestLevel.FAST):
    if message is None:
        message = "tests skipped"
    log.warning('\n'
        'Running in {0} Testing Level. Skipping {1} tests in {2}: {3}\n'.format(
                TestLevel.name(get_testing_level()),
                TestLevel.name(level),
                module_name,
                message))

def get_testing_level():
    l = os.environ.get("PYMSBAYES_TESTING_LEVEL")
    if l is None:
        return TestLevel.NORMAL
    try:
        return TestLevel.name_to_int(l)
    except:
        _LOG.warn('the value {0} for PYMSBAYES_TESTING_LEVEL is not recognized.'
                  'Using NORMAL level'.format(l))
    return TestLevel.NORMAL

def test_enabled(level, log=None, module_name="", message=None):
    tl = get_testing_level()
    if level > tl:
        if log:
            test_skip_msg(log, module_name, message, level)
        return False
    return True

if __name__ == '__main__':
    run()

