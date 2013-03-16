#! /usr/bin/env python

class PyMsBayesError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class TempFSError(PyMsBayesError):
    def __init__(self, *args, **kwargs):
        PyMsBayesError.__init__(self, *args, **kwargs)

class WorkerExecutionError(PyMsBayesError):
    def __init__(self, *args, **kwargs):
        PyMsBayesError.__init__(self, *args, **kwargs)

