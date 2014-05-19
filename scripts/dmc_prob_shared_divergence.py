#! /usr/bin/env python

"""
CLI for estimating the probability of given taxon pairs co-diverging.
This script is deprecated and superseded by `dmc_posterior_probs.py`.
"""

import os
import sys
import math
import random
import multiprocessing
import argparse

from pymsbayes.utils import argparse_utils

_program_info = {
    'name': os.path.basename(__file__),
    'author': 'Jamie Oaks',
    'version': 'Version 0.1.1',
    'description': __doc__,
    'copyright': 'Copyright (C) 2013 Jamie Oaks',
    'license': 'GNU GPL version 3 or later',}

def main_cli():
    description = '{name} {version}'.format(**_program_info)
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('div_model_path',
            metavar='DIV-MODEL-RESULTS-FILE',
            type=argparse_utils.arg_is_file,
            help = ('Path to divergence model results file (i.e., '
                    '`*-div-model-results.txt`).'))
    parser.add_argument('-i', '--taxon-indices',
            nargs = '+',
            type = argparse_utils.arg_is_positive_int,
            required = True,
            help = ('Two or more space-separated indices of taxa for which to '
                    'calculate the probability of them co-diverging. Indices '
                    'correspond to the line in the sample table of the config, '
                    'starting at 1 for the first line of the table. At least '
                    'two indices are required.'))
    parser.add_argument('-c', '--config',
            type = argparse_utils.arg_is_config,
            help = ('msBayes config file to be used to estimate prior '
                    'probability via simulations. If provided, the '
                    'posterior and prior probability and bayes factor is '
                    'reported. If not provided, only the posterior '
                    'probability is reported.'))
    parser.add_argument('-n', '--num-prior-samples',
            action = 'store',
            type = argparse_utils.arg_is_positive_int,
            default = 100000,
            help = ('The number of prior samples to simulate for estimating'
                    'prior probabilities. Only used if a config file is '
                    'provided with the `-c` argument.'))
    parser.add_argument('--np',
            action = 'store',
            type = argparse_utils.arg_is_positive_int,
            default = multiprocessing.cpu_count(),
            help = ('The maximum number of processes to run in parallel for '
                    'prior simulations. The default is the number of CPUs '
                    'available on the machine. This option is only relevant '
                    'if a config file is provided using the `-c` argument.'))
    parser.add_argument('--seed',
            action = 'store',
            type = argparse_utils.arg_is_positive_int,
            help = ('Random number seed to use for simulations. This option '
                    'is only relevant if a config file is provided using the '
                    '`-c` argument.'))
    parser.add_argument('--version',
            action = 'version',
            version = '%(prog)s ' + _program_info['version'],
            help = 'Report version and exit.')
    parser.add_argument('--quiet',
            action = 'store_true',
            help = 'Run without verbose messaging.')
    parser.add_argument('--debug',
            action = 'store_true',
            help = 'Run in debugging mode.')

    args = parser.parse_args()

    ##########################################################################
    ## handle args

    from pymsbayes.utils.messaging import (LoggingControl,
            InfoLogger)

    LoggingControl.set_logging_level("INFO")
    if args.quiet:
        LoggingControl.set_logging_level("WARNING")
    if args.debug:
        LoggingControl.set_logging_level("DEBUG")
    log = LoggingControl.get_logger(__name__)

    from pymsbayes import config
    from pymsbayes.teams import ModelProbabilityEstimatorTeam
    from pymsbayes.utils import parsing, GLOBAL_RNG

    if len(args.taxon_indices) < 2:
        log.error('At least two taxon indices are required')
        sys.exit(1)
    if not args.seed:
        args.seed = random.randint(1, 999999999)
    GLOBAL_RNG.seed(args.seed)

    div_models = parsing.OrderedDivergenceModelCollection(
            div_model_results_path = args.div_model_path)
    for i in args.taxon_indices:
        if ((i < 1) or (i > div_models.npairs)):
            log.error('taxon index {0} is out of bounds'.format(i))
            sys.exit(1)
    args.taxon_indices = [i - 1 for i in args.taxon_indices]
    prob_shared_div = div_models.prob_of_shared_divergence(args.taxon_indices)

    if args.config:
        prob_estimator_team = ModelProbabilityEstimatorTeam(
                config_paths = [args.config],
                num_samples = args.num_prior_samples,
                num_processors = args.np)
        prob_estimator_team.start()
        prior_prob = prob_estimator_team.shared_div_probs[args.config][
                len(args.taxon_indices)]
        bf = ((prob_shared_div / (1 - prob_shared_div)) /
                (prior_prob / (1 - prior_prob)))

    sys.stdout.write('posterior probability = {0}\n'.format(prob_shared_div))
    if args.config:
        sys.stdout.write('prior probability = {0}\n'.format(prior_prob))
        sys.stdout.write('Bayes factor = {0}\n'.format(bf))
        sys.stdout.write('2ln(Bayes factor) = {0}\n'.format(2 * math.log(bf)))
if __name__ == '__main__':
    main_cli()

