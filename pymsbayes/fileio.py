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
    if kwargs.has_key('name'):
        kwargs['name'] = expand_path(kwargs['name'])
        return FileStream(*args, **kwargs)
    path = expand_path(args[0])
    return FileStream(path, *args[1:], **kwargs)

def process_file_arg(file_arg, mode='rU', compresslevel = None):
    close = False
    file_stream = file_arg
    if isinstance(file_arg, str):
        if compresslevel:
            file_stream = GzipFileStream(expand_path(file_arg), mode,
                    compresslevel)
        else:
            file_stream = open(file_arg, mode)
        close = True
    return file_stream, close

def is_gzip(file_arg):
    fs, close = process_file_arg(file_arg)
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

