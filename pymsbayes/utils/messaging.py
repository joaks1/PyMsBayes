#! /usr/bin/env python

import sys
import os
import logging

_LOGGING_LEVEL_ENV_VAR = "PYMSBAYES_LOGGING_LEVEL"

def get_logging_level(level=None):
    if level:
        if level.upper() == "NOTSET":
            return logging.NOTSET
        elif level.upper() == "DEBUG":
            return logging.DEBUG
        elif level.upper() == "INFO":
            return logging.INFO
        elif level.upper() == "WARNING":
            return logging.WARNING
        elif level.upper() == "ERROR":
            return logging.ERROR
        elif level.upper() == "CRITICAL":
            return logging.CRITICAL
        else:
            return logging.WARNING
    else:
        return logging.WARNING

def get_env_logging_level():
    return get_logging_level(os.environ.get(_LOGGING_LEVEL_ENV_VAR, None))

def get_logger(name = 'pymsbayes', level=None):
    if level:
        l = get_logging_level(level)
    else:
        l = get_env_logging_level()
    log = logging.getLogger(name)
    log.setLevel(l)
    h = logging.StreamHandler()
    h.setLevel(l)
    log.addHandler(h)
    return log
