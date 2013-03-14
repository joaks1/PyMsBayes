#! /usr/bin/env python

import sys
import os
import errno
import random
import platform
import string

GLOBAL_RNG = random.Random()
PLATFORM = platform.system().lower()

PACKAGE_DIR = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
BASE_DIR = os.path.abspath(os.path.dirname(PACKAGE_DIR))
SCRIPTS_DIR = os.path.join(BASE_DIR, "scripts")
BIN_DIR = None
if PLATFORM == 'linux':
    BIN_DIR = os.path.join(BASE_DIR, "bin", "linux", "bin")
elif PLATFORM == 'darwin':
    BIN_DIR = os.path.join(base_dir, "bin", "mac", "bin")
elif PLATFORM == 'windows':
    BIN_DIR = os.path.join(base_dir, "bin", "win", "bin")

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

