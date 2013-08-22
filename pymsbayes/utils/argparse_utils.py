#! /usr/bin/env python

import os
import argparse

from pymsbayes.fileio import expand_path
from pymsbayes.config import MsBayesConfig
from pymsbayes.utils.functions import is_file, is_dir

def arg_is_path(path):
    try:
        if not os.path.exists(path):
            raise
    except:
        msg = 'path {0!r} does not exist'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return expand_path(path)

def arg_is_file(path):
    try:
        if not is_file(path):
            raise
    except:
        msg = '{0!r} is not a file'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return expand_path(path)

def arg_is_config(path):
    try:
        if not MsBayesConfig.is_config(path):
            raise
    except:
        msg = '{0!r} is not an msBayes config file'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return expand_path(path)

def arg_is_dir(path):
    try:
        if not is_dir(path):
            raise
    except:
        msg = '{0!r} is not a directory'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return expand_path(path)

