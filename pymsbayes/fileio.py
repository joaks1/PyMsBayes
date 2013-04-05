#! /usr/bin/env python

import os
import sys
from __builtin__ import file
from __builtin__ import open as fopen

from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)


def expand_path(path):
    return os.path.abspath(os.path.expandvars(os.path.expanduser(path)))

def open(*args, **kwargs):
    path = expand_path(args[0])
    return FileStream(path, *args[1:], **kwargs)

def process_file_arg(file_arg, mode='rU'):
    close = False
    file_stream = file_arg
    if isinstance(file_arg, str):
        file_stream = open(file_arg, mode)
        close = True
    return file_stream, close

class FileStream(file):
    open_files = set()
    def __init__(self, *args, **kwargs):
        file.__init__(self, *args, **kwargs)
        self.__class__.open_files.add(self)

    def close(self):
        file.close(self)
        self.__class__.open_files.remove(self)
