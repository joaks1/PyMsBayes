#! /usr/bin/env python

import sys
import os
import platform
import random
import multiprocessing

WORK_FORCE = multiprocessing.Queue()
GLOBAL_RNG = random.Random()
PLATFORM = platform.system().lower()

PACKAGE_DIR = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
BASE_DIR = os.path.abspath(os.path.dirname(PACKAGE_DIR))
SCRIPTS_DIR = os.path.join(BASE_DIR, "scripts")
BIN_DIR = None
if PLATFORM == 'linux':
    BIN_DIR = os.path.join(BASE_DIR, "bin", "linux")
elif PLATFORM == 'darwin':
    BIN_DIR = os.path.join(BASE_DIR, "bin", "mac")
elif PLATFORM == 'windows':
    BIN_DIR = os.path.join(BASE_DIR, "bin", "win")

def get_tool_path_mapping():
    tool_path_map = {}
    tool_path_map['abcestimator'] = os.path.join(BIN_DIR, 'ABCestimator')
    tool_path_map['eureject'] = os.path.join(BIN_DIR, 'eureject')
    tool_path_map['msreject'] = os.path.join(BIN_DIR, 'msReject')
    tool_path_map['regress_cli'] = os.path.join(BIN_DIR, 'regress_cli.r')
    tool_path_map['obssumstats'] = os.path.join(BIN_DIR, 'new', 'obsSumStats.pl')
    tool_path_map['msbayes'] = os.path.join(BIN_DIR, 'new', 'msbayes.pl')
    tool_path_map['msbayes-old'] = os.path.join(BIN_DIR, 'old', 'msbayes.pl')
    for name, path in tool_path_map.iteritems():
        if not os.path.exists(path):
            raise Exception('The path {0!r} for {1!r} does not exist'.format(
                    path, name))
    return tool_path_map

TOOL_PATH_MAP = get_tool_path_mapping()

def get_tool_path(name):
    if not TOOL_PATH_MAP.has_key(name.lower()):
        raise Exception(
                '{0!r} is not a valid tool. Valid tools include: {1}'.format(
                        name.lower(), ', '.join(TOOL_PATH_MAP.keys())))
    return TOOL_PATH_MAP[name.lower()]

