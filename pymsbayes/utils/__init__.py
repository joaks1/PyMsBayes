#! /usr/bin/env python

import sys
import os
import platform
import random
import subprocess
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

PACKAGE_DIR = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
BASE_DIR = os.path.abspath(os.path.dirname(PACKAGE_DIR))
SCRIPTS_DIR = os.path.join(BASE_DIR, "scripts")


class ToolPathManager(object):
    _ignore_internal_tools = False 
    _package_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
    _base_dir = os.path.abspath(os.path.dirname(_package_dir))
    _platform = platform.system().lower()
    _bin_dir = None
    if _platform == 'linux':
        _bin_dir = os.path.join(_base_dir, "bin", "linux")
    elif _platform == 'darwin':
        _bin_dir = os.path.join(_base_dir, "bin", "mac")
    elif _platform == 'windows':
        _bin_dir = os.path.join(_base_dir, "bin", "win")
    class ToolNotFoundError(Exception):
        def __init__(self, *args, **kwargs):
            Exception.__init__(self, *args, **kwargs)

    @classmethod
    def set_to_ignore_internal_tools(cls):
        cls._ignore_internal_tools = True

    @classmethod
    def set_to_prefer_internal_tools(cls):
        cls._ignore_internal_tools = False

    @classmethod
    def ignoring_internal_tools(cls):
        return cls._ignore_internal_tools

    @classmethod
    def get_tool_path(cls, name):
        """
        The primary package-wide method for finding exe paths for subprocess
        calls. The priority of the search is (1) use full path if given, (2)
        check internal bundled tools, (3) check system. The logic is as
        follows:
            - if `name` is a full path and executable, `name` is returned.
            - next, if `cls._ignore_internal_tools` is `False`, the bundled
              tool bin is checked for `name`, and if an exe is found, the full
              path is returned.
            - last, check the PATH for `name`, and if found return `name`.
        """
        # top priority given to full paths
        if cls.is_executable(name):
            return name
        # next check bundled tools if allowed
        if not cls._ignore_internal_tools:
            internal_path = os.path.join(cls._bin_dir, name)
            if cls.is_executable(internal_path):
                return internal_path
        # last check system
        external_path = cls.get_external_tool(name)
        if external_path:
            return external_path
        # not found: raise error here so that the multiprocessing
        # machinery does not get flooded with jobs doomed to die.
        raise cls.ToolNotFoundError('Unable to find executable '
                '{0!r}'.format(name))

    @classmethod
    def get_external_tool(cls, path):
        """
        Uses `subprocess.check_call` to check system for `exe_name`. If found,
        `exe_name` is returned, else `None` is returned.

        """
        try:
            p = subprocess.Popen([path],
                    stdout = subprocess.PIPE,
                    stderr = subprocess.PIPE)
            p.terminate()
        except subprocess.CalledProcessError:
            return path
        except OSError:
            return None
        return path

    @classmethod
    def is_executable(cls, path):
        return os.path.isfile(path) and bool(cls.get_external_tool(path))

    @classmethod
    def which(cls, exe):
        if cls.is_executable(exe):
            return exe
        name = os.path.basename(exe)
        for p in os.environ['PATH'].split(os.pathsep):
            p = p.strip('"')
            exe_path = os.path.join(p, name)
            if cls.is_executable(exe_path):
                return exe_path
        return None


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

