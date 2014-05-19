#! /usr/bin/env python

import sys
import os
import logging

class LoggingControl(object):
    _internal_level = None
    logging_level_env_var = "PYMSBAYES_LOGGING_LEVEL"

    @classmethod
    def get_logging_level(cls, level = None):
        if not level:
            level = os.environ.get(cls.logging_level_env_var, None)
            # allow internal override of env logging level
            if cls._internal_level:
                level = cls._internal_level
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
        return logging.WARNING

    @classmethod
    def set_logging_level(cls, level):
        cls._internal_level = level
    
    @classmethod
    def get_logger(cls, name, level = None):
        l = cls.get_logging_level(level)
        log = logging.getLogger(name)
        log.setLevel(l)
        h = logging.StreamHandler()
        h.setLevel(l)
        log.addHandler(h)
        return log

get_logger = LoggingControl.get_logger

class InfoLogger(object):
    def __init__(self, path):
        self.path = path

    def write(self, msg, log_func = None):
        out = open(self.path, 'a')
        out.write(msg + os.linesep)
        out.close()
        if log_func:
            log_func(msg)

