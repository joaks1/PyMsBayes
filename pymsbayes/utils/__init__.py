#! /usr/bin/env python

import sys
import os
import errno
import random
import platform
import string

GLOBAL_RNG = random.Random()
PLATFORM = platform.system().lower()

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

def random_str(self, length=8,
        char_pool=string.ascii_letters + string.digits):
    return ''.join(random.choice(char_pool) for i in range(length))

