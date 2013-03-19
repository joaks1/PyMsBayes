#! /usr/bin/env python

import sys
import os
import errno
import random
import string
import stat

from pymsbayes.utils import GLOBAL_RNG

def mkdr(path):
    """
    Creates directory `path`, but suppresses error if `path` already exists.
    """
    try:
        os.makedirs(path)
    except OSError, e:
        if e.errno == errno.EEXIST:
            pass
        else:
            raise

def random_str(length=8,
        char_pool=string.ascii_letters + string.digits):
    return ''.join(random.choice(char_pool) for i in range(length))

def get_random_int():
    return GLOBAL_RNG.randint(1, 999999999)

def expand_path(path):
    return os.path.abspath(os.path.expandvars(os.path.expanduser(path)))

def process_file_arg(file_arg, mode='rU'):
    close = False
    file_stream = file_arg
    if isinstance(file_arg, str):
        file_stream = open(file_arg, mode)
        close = True
    return file_stream, close

def get_indices_of_patterns(target_list, regex_list):
    indices = []
    for regex in regex_list:
        indices.extend([i for i, e in enumerate(target_list) if regex.match(e)])
    return sorted(indices)

def get_indices_of_strings(target_list, string_list):
    indices = []
    for s in string_list:
        indices.extend([i for i, e in enumerate(target_list) if s.strip() == e.strip()])
    return sorted(indices)
    
def reduce_columns(in_file, out_file, column_indices, sep='\t',
        extra_tab=False):
    in_stream, close_in = process_file_arg(in_file, 'rU')
    out_stream, close_out = process_file_arg(out_file, 'w')
    line_iter = iter(in_stream)
    for line_num, line in enumerate(line_iter):
        l = line.strip().split(sep)
        new_line = [l[i] for i in column_indices]
        if extra_tab:
            out_stream.write('%s\t\n' % sep.join(new_line))
        else:
            out_stream.write('%s\n' % sep.join(new_line))
    if close_in:
        in_stream.close()
    if close_out:
        out_stream.close()

def list_splitter(l, n, by_size=False):
    """
    Returns generator that yields list `l` as `n` sublists, or as `n`-sized
    sublists if `by_size` is True.
    """
    if by_size:
        for i in range(0, len(l), n):
            yield l[i:i+n]
    else:
        for i in range(0, len(l), (len(l)/int(n))):
            yield l[i:i+(len(l)/int(n))]

def whereis(file_name):
    """
    Returns the first absolute path to `file_name` encountered in $PATH.
    Returns `None` if `file_name` is not found in $PATH.
    """
    paths = os.environ.get('PATH', '').split(':')
    for path in paths:
        abs_path = os.path.join(path, file_name)
        if os.path.exists(abs_path) and not os.path.isdir(abs_path):
            return abs_path
            break
    return None

def is_file(path):
    if not path:
        return False
    if not os.path.isfile(path):
        return False
    return True

def is_executable(path):
    is_f = is_file(path)
    if not is_f:
        return False
    if (os.stat(path)[stat.ST_MODE] & (stat.S_IXUSR|stat.S_IXGRP|stat.S_IXOTH)) == 0:
        return False
    return True

