#! /usr/bin/env python

import sys
import os
import platform
import random
import multiprocessing
import time
import atexit
try:
    from guppy import hpy
    HEAPY = True
except ImportError:
    HEAPY = False

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
    tool_path_map['obssumstats'] = os.path.join(BIN_DIR, 'obsSumStats.pl')
    tool_path_map['dpp-msbayes'] = os.path.join(BIN_DIR, 'dpp-msbayes.pl')
    tool_path_map['msbayes'] = os.path.join(BIN_DIR, 'msbayes.pl')
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


class MSBAYES_SORT_INDEX(object):
    valid_values = list(range(12))
    _default = 0
    _i = _default

    @classmethod
    def set_index(cls, i):
        i = int(i)
        if i not in cls.valid_values:
            raise Exception('MsBayesSortIndex {0} is not valid; valid options '
                    'are: {1}'.format(i,
                            ','.join([str(x) for x in cls.valid_values])))
        cls._i = i
        
    @classmethod
    def current_value(cls):
        return cls._i

    @classmethod
    def reset_default(cls):
        cls.set_index(cls._default)

MEMORY_LOGGING_ENV_VAR = "PYMSBAYES_MEMORY_LOGGING_FREQUENCY"
MEMORY_LOGGING_FREQUENCY = float(os.environ.get(MEMORY_LOGGING_ENV_VAR, 0))
DUMP_DEBUG_INFO_ENV_VAR = "PYMSBAYES_DUMP_DEBUG_INFO"
DUMP_DEBUG_INFO = bool(int(os.environ.get(DUMP_DEBUG_INFO_ENV_VAR, 0)))
_LAST_MEMORY_LOG_TIME = time.time()
_MEMORY_LOG_PATH = 'pymsbayes-memory.log'
TRACE_LINES_INTO = []

def memory_trace_lines(frame, event, arg):
    if event != 'line':
        return
    co = frame.f_code
    func_name = co.co_name
    func_line_num = frame.f_lineno
    file_name = co.co_filename
    hp = hpy().heap()
    with open(_MEMORY_LOG_PATH, 'a') as out:
        out.write('{0}: {1}: {2}\n{3}\n\n'.format(
                file_name,
                func_name,
                func_line_num,
                hp))
    return

def memory_trace_calls(frame, event, arg):
    global _LAST_MEMORY_LOG_TIME
    t = time.time()
    if t - _LAST_MEMORY_LOG_TIME < MEMORY_LOGGING_FREQUENCY:
        return
    _LAST_MEMORY_LOG_TIME = t
    if event != 'call':
        return
    co = frame.f_code
    func_name = co.co_name
    if func_name == 'write':
        return
    func_line_num = frame.f_lineno
    file_name = co.co_filename
    caller = frame.f_back
    caller_line_num = caller.f_lineno
    caller_file_name = caller.f_code.co_filename
    hp = hpy().heap()
    with open(_MEMORY_LOG_PATH, 'a') as out:
        out.write('Call to {0}: {1}: {2} (from {3}: {4})\n'.format(
            file_name,
            func_name,
            func_line_num,
            caller_file_name,
            caller_line_num))
        dump_memory_info(out)
        out.write('\n')
    if func_name in TRACE_LINES_INTO:
        memory_trace_lines
    return

def dump_memory_info(stream = None):
    close = False
    if not stream:
        stream = open(_MEMORY_LOG_PATH, 'a')
        close = True
    if HEAPY:
        hp = hpy().heap()
        stream.write('heap:\n')
        stream.write('{0}\n'.format(hp))
        stream.write('heap.byrcs:\n')
        stream.write('{0}\n'.format(hp.byrcs))
        stream.write('heap[0].byrcs:\n')
        stream.write('{0}\n'.format(hp[0].byrcs))
        stream.write('heap[0].byid:\n')
        stream.write('{0}\n'.format(hp[0].byid))
        stream.write('heap[0].byvia:\n')
        stream.write('{0}\n'.format(hp[0].byvia))
        stream.write('heap[0].byrcs[0].referrers.byrcs:\n')
        stream.write('{0}\n'.format(hp[0].byrcs[0].referrers.byrcs))
        stream.write('heap[0].byrcs[0].referrers.byrcs[0].referents:\n')
        stream.write('{0}\n'.format(hp[0].byrcs[0].referrers.byrcs[0].referents))
        stream.write('heap[0].byrcs[0].referrers.byrcs[0].referents.byvia:\n')
        stream.write('{0}\n'.format(hp[0].byrcs[0].referrers.byrcs[0].referents.byvia))
    else:
        stream.write('hpy (heapy) from the guppy package is not available\n')
    if close:
        stream.close()

def set_memory_trace():
    if HEAPY and (MEMORY_LOGGING_FREQUENCY > 0):
        sys.settrace(memory_trace_calls)

