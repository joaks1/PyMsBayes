#! /usr/bin/env python

import os
import argparse

from pymsbayes.fileio import expand_path
from pymsbayes.config import MsBayesConfig
from pymsbayes.utils.functions import is_file, is_dir, is_executable

class SmartHelpFormatter(argparse.HelpFormatter):
    '''
    A class to allow customizable line breaks for an argument help message
    on a per argument basis.
    '''

    def _split_lines(self, text, width):
        if text.startswith('r|'):
            return text[2:].splitlines()
        return argparse.HelpFormatter._split_lines(self, text, width)

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

def arg_is_executable(path):
    try:
        p = ToolPathManager.get_external_tool(path)
        if not p:
            raise
    except:
        msg = '{0!r} is not an executable'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return p

def arg_is_nonnegative_int(i):
    try:
        if int(i) < 0:
            raise
    except:
        msg = '{0!r} is not a non-negative integer'.format(i)
        raise argparse.ArgumentTypeError(msg)
    return int(i)

def arg_is_positive_int(i):
    try:
        if int(i) < 1:
            raise
    except:
        msg = '{0!r} is not a positive integer'.format(i)
        raise argparse.ArgumentTypeError(msg)
    return int(i)

def arg_is_positive_float(i):
    try:
        if float(i) <= 0.0:
            raise
    except:
        msg = '{0!r} is not a positive real number'.format(i)
        raise argparse.ArgumentTypeError(msg)
    return float(i)

def get_sort_index_help_message():
    return (
'''r|The sorting index used by
`dpp-msbayes.pl`/`msbayes.pl` and `obsSumStats.pl`
scripts to determine how the summary statistic vectors
calculated from the alignments of the observed and
simulated data are to be grouped and sorted.
The default is %(default)s.
0:    Do not group or sort. The identity and order of
      the summary statistics of each alignment are
      maintained and compared when calculating
      Euclidean distance.
1-7:  **NOTE**, options 1-7 all re-sort the summary
      statistics in some way, and thus compare the
      statistics from *different* alignments when
      calculating the Euclidean distance.  This is not
      valid and these options should *NOT* be used.
      They are maintained for backwards compatibility
      with the original msBayes.
8-11: All of these options group the summary
      statistics from multiple loci by taxon and then
      calculate moments of each statistic across the
      loci for each taxon, and then use these moments
      to calculate Euclidean distance. The order of
      the taxa is maintained, and so this is valid,
      but you are losing a lot of information
      contained in your loci by simply taking the mean
      (option 11) across them. If you have A LOT of
      loci, this sacrifice might be necessary to
      reduce the number of summary statistics.
      **NOTE**, options 8-10 are NOT well tested.
      8:  Use the first 4 moments (mean, variance,
          skewness, and kurtosis) of each statistic.
      9:  Use the first 3 moments (mean, variance,
          and skewness) of each statistic.
      10: Use the first 2 moments (mean and variance)
          of each statistic.
      11: Use the first 1 moment (mean) of each
          statistic.'''
)

