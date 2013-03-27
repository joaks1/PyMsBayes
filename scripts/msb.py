#! /usr/bin/env python

"""
Main CLI for PyMsBayes package.
"""

import os
import sys
import argparse

from pymsbayes.utils import WORK_FORCE, GLOBAL_RNG
from pymsbayes.utils.functions import is_file, is_dir, expand_path
from pymsbayes.utils.messaging import get_logger

_program_info = {
    'name': os.path.basename(__file__),
    'author': 'Jamie Oaks',
    'version': 'Version 0.1.0',
    'description': __doc__,
    'copyright': 'Copyright (C) 2013 Jamie Oaks',
    'license': 'GNU GPL version 3 or later',}

_LOG = get_logger(__name__, 'INFO')


def arg_is_file(path):
    if not is_file(path):
        msg = '{0!r} is not a file'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return expand_path(path)

def arg_is_dir(path):
    if not is_dir(path):
        msg = '{0!r} is not a directory'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return expand_path(path)

def main_cli():
    description = '{name} {version}'.format(**_program_info)
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('-o', '--observed-config',
            action = 'store',
            type = arg_is_file,
            required = True,
            help = ('The msBayes config file that specifies the observed '
                    'data. If used in combination with `-r` this config will '
                    'be used to simulate pseudo-observed data. If analyzing '
                    'real data, do not use the `-r` option, and the data '
                    'files specified in the config must exist and contain the '
                    'sequence data.'))
    parser.add_argument('-p', '--prior-config',
            nargs = '+',
            type = arg_is_file,
            required = True,
            help = ('One or more config files to be used to generate prior '
                    'samples. If more than one config is specified, they '
                    'should be separated by spaces.'))
    parser.add_argument('-r', '--reps',
            action = 'store',
            type = int,
            help = ('This option has two effects. First, it signifies that '
                    'the analysis will be simulation based (i.e., no real '
                    'data will be used). Second, it specifies how many '
                    'simulation replicates to perform (i.e., how many data '
                    'sets to simulate and analyze).'))
    parser.add_argument('--output-dir',
            action = 'store',
            type = arg_is_dir,
            help = ('The directory in which all output files will be written. '
                    'The default is to use the directory of the observed '
                    'config file.'))
    parser.add_argument('-s', '--stat-prefixes',
            nargs = '*',
            type = str,
            help = ('Prefixes of summary statistics to use in the analyses. '
                    'The prefixes should be separated by spaces. '
                    'Default: `-s pi wattTheta pi.net tajD.denom`.'))
    parser.add_argument('-c', '--continuous-prefixes',
            nargs = '*',
            type = str,
            help = ('Prefixes of continuous parameters on which regression '
                    'adjustment should be performed. '
                    'The prefixes should be separated by spaces. '
                    'Default: `-c PRI.E.t PRI.omega`.'))
    parser.add_argument('-d', '--discrete-prefixes',
            nargs = '*',
            type = str,
            help = ('Prefixes of discrete parameters on which regression '
                    'adjustment should be performed. '
                    'The prefixes should be separated by spaces. '
                    'Default: `-c PRI.model PRI.Psi`.'))
    parser.add_argument('-m', '--merge-priors',
            action = 'store_true',
            default = 'False',
            help = ('Merge the priors generated from the specified configs '
                    'together before analysis. Otherwise, the observed data '
                    'sets will be analyzed under each prior separately.'))
    parser.add_argument('-k', '--keep-priors',
            action = 'store_true',
            default = 'False',
            help = 'Keep prior-sample files.')
    parser.add_argument('--keep-temps',
            action = 'store_true',
            default = 'False',
            help = 'Keep all temporary files.')
    parser.add_argument('--seed',
            action = 'store',
            type = int,
            help = 'Random number seed to use for the analysis.')
    parser.add_argument('--dry-run',
            action = 'store_true',
            default = 'False',
            help = 'Report configuration of the analysis and exit.')
    parser.add_argument('--version',
            action = 'version',
            version = '%(prog)s ' + _program_info['version'],
            help = 'Report version and exit.')
    parser.add_argument('-v', '--verbose',
            action = 'store_true',
            default = 'False',
            help = 'Run with verbose messaging.')
    parser.add_argument('--debug',
            action = 'store_true',
            default = 'False',
            help = 'Run in debugging mode.')

    args = parser.parse_args()
    print args

if __name__ == '__main__':
    main_cli()

