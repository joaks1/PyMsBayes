#! /usr/bin/env python

class PyMsBayesError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class ArgumentError(PyMsBayesError):
    def __init__(self, *args, **kwargs):
        PyMsBayesError.__init__(self, *args, **kwargs)

class TempFSError(PyMsBayesError):
    def __init__(self, *args, **kwargs):
        PyMsBayesError.__init__(self, *args, **kwargs)

class WorkerExecutionError(PyMsBayesError):
    def __init__(self, *args, **kwargs):
        PyMsBayesError.__init__(self, *args, **kwargs)

class PriorMergeError(PyMsBayesError):
    def __init__(self, *args, **kwargs):
        PyMsBayesError.__init__(self, *args, **kwargs)

class SummaryFileParsingError(PyMsBayesError):
    def __init__(self, *args, **kwargs):
        PyMsBayesError.__init__(self, *args, **kwargs)

class ParameterParsingError(PyMsBayesError):
    def __init__(self, *args, **kwargs):
        PyMsBayesError.__init__(self, *args, **kwargs)

class SampleTableRowError(PyMsBayesError):
    def __init__(self, *args, **kwargs):
        PyMsBayesError.__init__(self, *args, **kwargs)

class SampleTableError(PyMsBayesError):
    def __init__(self, *args, **kwargs):
        PyMsBayesError.__init__(self, *args, **kwargs)
