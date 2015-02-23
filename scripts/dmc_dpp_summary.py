#!/usr/bin/env python

"""
Dirichlet process calculations.
"""

import os
import sys
import argparse
import random
import multiprocessing

from pymsbayes.utils import argparse_utils

_program_info = {
    'name': os.path.basename(__file__),
    'author': 'Jamie Oaks',
    'version': 'Version 0.1.0',
    'description': __doc__,
    'copyright': 'Copyright (C) 2013 Jamie Oaks',
    'license': 'GNU GPL version 3 or later',}

def main():
    keys_to_print = {
            'ncats': 'expected number of categories',
            'concentration': 'concentration parameter',
    }
    parameter_options = ['concentration', 'ncats']
    description = '{name} {version}'.format(**_program_info)
    parser = argparse.ArgumentParser(description = description,
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--version',
            action='version',
            version='%(prog)s ' + _program_info['version'],
            help='report version and exit')
    parser.add_argument('parameter',
            choices=parameter_options,
            nargs=1,
            help = ('The parameter provided. The two options are:\n'
                    '`concentration`: the concentration parameter\n'
                        '\tof the Dirichlet process.\n'
                    '`ncats`: the expected (mean) number of\n'
                        '\tcategories for the dirichlet\n'
                        '\tprocess.\n'
                    'You provide one of these two parameters along\n'
                    'with the number of elements (taxon pairs),\n'
                    'and this program calculates and returns the\n'
                    'other one accordingly.'))
    parser.add_argument('parameter_value',
            metavar='X',
            type=argparse_utils.arg_is_positive_float,
            help=('Value of the parameter'))
    parser.add_argument('--shape',
            required = False,
            type=argparse_utils.arg_is_positive_float,
            help = ('Shape parameter of a gamma hyperprior on the\n'
                    'concentration parameter of the Dirichlet\n'
                    'process. If provided, the program will\n'
                    'calculate a corresponding scale parameter\n'
                    'for the gamma hyperprior such that the\n'
                    'mean of the gamma hyperprior is equal to\n'
                    'the reported concentration parameter and the\n'
                    'prior expectation for the number of\n'
                    'categories is equal to `ncats`.'))
    parser.add_argument('--reps',
            action = 'store',
            type = int,
            required = False,
            help = ('The number of simulation replicates to use for\n'
                    'estimating the probability of the number of categories.\n'
                    'By default, no simulations are run and the probabilities\n'
                    'are not estimated or reported.'))
    parser.add_argument('--np',
            action = 'store',
            type = argparse_utils.arg_is_positive_int,
            default = multiprocessing.cpu_count(),
            help = ('The maximum number of processes to run in parallel. The\n'
                    'default is the number of CPUs available on the machine.\n'
                    'This option is only used if `--reps` is specified.'))
    parser.add_argument('--seed',
            action = 'store',
            type = argparse_utils.arg_is_positive_int,
            help = 'Random number seed to use for the analysis.')
    parser.add_argument('num_elements',
            metavar='N',
            type=argparse_utils.arg_is_nonnegative_int,
            help='Number of elements (i.e., number of taxon pairs).')

    args = parser.parse_args()

    from pymsbayes.utils.messaging import (LoggingControl,
            InfoLogger)

    LoggingControl.set_logging_level("WARNING")


    from pymsbayes.teams import DppSimTeam
    from pymsbayes.utils import GLOBAL_RNG, probability
    from pymsbayes.utils.stats import Partition

    p = Partition('0' * args.num_elements)

    results = dict(zip(parameter_options,
            [None for k in parameter_options]))
    args.parameter = args.parameter[0]
    if args.parameter == 'concentration':
        results['concentration'] = args.parameter_value
        results['ncats'] = p.get_dpp_expected_num_cats(args.parameter_value)

    elif args.parameter == 'ncats':
        if args.parameter_value > args.num_elements:
            sys.stderr.write('ERROR: `ncats` cannot be greater than the number '
                    'of elements\n')
            sys.exit(1)
        elif args.parameter_value < 1.0:
            sys.stderr.write('ERROR: `ncats` cannot be less than 1\n')
            sys.exit(1)
        results['ncats'] = args.parameter_value
        results['concentration'] = p.get_dpp_concentration(args.parameter_value)

    else:
        raise Exception('parameter option {0} is not valid'.format(
                args.parameter))

    alpha = results['concentration']
    if args.shape:
        results['shape'] = args.shape
        results['scale'] = results['concentration'] / args.shape
        parameter_options.extend(['shape', 'scale'])
        alpha = probability.GammaDistribution(
                shape = results['shape'],
                scale = results['scale'])
        
    sys.stdout.write('number of elements = {0}\n'.format(args.num_elements))
    for key in parameter_options:
        sys.stdout.write('{0} = {1}\n'.format(
                keys_to_print.get(key, key),
                results[key]))

    if args.reps:
        sys.stderr.write('\nStarting simulations to estimate probabilities...\n')
        if not args.seed:
            args.seed = random.randint(1, 999999999)
        sys.stderr.write('Using seed {0}\n\n'.format(args.seed))
        GLOBAL_RNG.seed(args.seed)

        sim_team = DppSimTeam(
                alpha = alpha,
                num_elements = args.num_elements,
                base_distribution = None,
                num_samples = args.reps,
                num_processors = args.np)
        sim_team.start()

        sys.stderr.write('Estimated probabilities of the number of categories:\n')
        for k, prob in sim_team.psi_probs.iteritems():
            sys.stdout.write('\tp(ncats = {0}) = {1:.4f} (n = {2})\n'.format(
                    k,
                    prob,
                    p.number_of_partitions_into_k_subsets(k)))

if __name__ == '__main__':
    main()

