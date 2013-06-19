#! /usr/bin/env python

import os
import sys
from gzip import GzipFile
from __builtin__ import file
from __builtin__ import open as fopen

from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)


def expand_path(path):
    return os.path.abspath(os.path.expandvars(os.path.expanduser(path)))

def open(*args, **kwargs):
    return FileStream(*args, **kwargs)

def process_file_arg(file_arg, mode='rU', compresslevel = None):
    if isinstance(file_arg, str):
        fp = expand_path(file_arg)
        gzipped = False
        if os.path.isfile(fp):
            gzipped = is_gzipped(fp)
        if compresslevel:
            return GzipFileStream(fp, mode, compresslevel), True
        elif gzipped:
            return GzipFileStream(fp, mode), True
        else:
            return open(fp, mode), True
    return file_arg, False

def is_gzipped(file_path):
    with open(expand_path(file_path)) as fs:
        l = fs.next()
        if l.startswith("\x1f\x8b"):
            return True
        return False

class FileStream(file):
    open_files = set()
    def __init__(self, *args, **kwargs):
        file.__init__(self, *args, **kwargs)
        self.__class__.open_files.add(self)

    def close(self):
        file.close(self)
        self.__class__.open_files.remove(self)

class GzipFileStream(GzipFile):
    open_files = set()
    def __init__(self, *args, **kwargs):
        GzipFile.__init__(self, *args, **kwargs)
        self.__class__.open_files.add(self)

    def close(self):
        GzipFile.close(self)
        self.__class__.open_files.remove(self)

