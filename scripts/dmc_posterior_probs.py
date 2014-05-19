#! /usr/bin/env python

"""
CLI for estimating the probability of given taxon pairs co-diverging.
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
    parser.add_argument('config',
            metavar='CONFIG-FILE',
            type = argparse_utils.arg_is_config,
            help = ('msBayes config file used to estimate the posterior '
                    'sample.'))
    parser.add_argument('posterior_sample_path',
            metavar='POSTERIOR-SAMPLE-FILE',
            type=argparse_utils.arg_is_file,
            help = ('Path to posterior sample file (i.e., '
                    '`*-posterior-sample.txt`).'))
    parser.add_argument('-e', '--expression',
            dest = 'expressions',
            action = 'append',
            metavar = 'TAXON-INDEX-EXPRESSION',
            type = str,
            required = True,
            help = ('A conditional expression of divergence times based on '
                    'the taxon-pair indices for which to calculate the '
                    'posterior probability of being true. Indices correspond '
                    'to the order that pairs of taxa appear in the sample '
                    'table of the config, starting at 0 for the first '
                    'taxon-pair to appear in the table (starting from the '
                    'top). E.g., `-e "0 == 3 == 4"` would request the '
                    'proportion of times the 1st, 4th, and 5th taxon-pairs '
                    '(in order of appearance in the sample table of the '
                    'config) share the same divergence time in the '
                    'posterior sample, whereas `-e "0 > 1" would request the '
                    'proportion of times the the 1st taxon-pair diverged '
                    'further back in time than the 2nd taxon-pair in the '
                    'posterior sample.'))
    parser.add_argument('-n', '--num-prior-samples',
            action = 'store',
            type = argparse_utils.arg_is_positive_int,
            help = ('The number of prior samples to simulate for estimating'
                    'prior probabilities; prior probabilities and Bayes '
                    'factors will be reported. The default is to only report '
                    'posterior probabilities.'))
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
    from pymsbayes.teams import DivModelSimulatorTeam
    from pymsbayes.utils import parsing, stats, GLOBAL_RNG

    if not args.seed:
        args.seed = random.randint(1, 999999999)
    GLOBAL_RNG.seed(args.seed)

    cfg = config.MsBayesConfig(args.config)

    evaluators = []
    for exp in args.expressions:
        evaluators.append(stats.ListConditionEvaluator(exp,
                index_labels = cfg.taxa))

    div_models = parsing.get_partitions_from_posterior_sample_file(
            args.posterior_sample_path)

    sim_team = None
    if args.num_prior_samples:
        sim_team = DivModelSimulatorTeam(
                config_paths = [args.config],
                num_samples = args.num_prior_samples,
                num_processors = args.np)
        sim_team.start()

    for e in evaluators:
        title = '{0} --- {1}:'.format(e.expression,
                e.pretty_expression)
        section_title = '\n{0}\n{1}\n'.format(title, '-' * len(title))
        sys.stdout.write('{0}'.format(section_title))
        prob_shared_div = div_models.get_condition_frequency(e)
        sys.stdout.write('posterior probability = {0}\n'.format(prob_shared_div))
        if sim_team:
            prior_prob = sim_team.div_models[
                    args.config].get_condition_frequency(e)
            bf = ((prob_shared_div / (1 - prob_shared_div)) /
                    (prior_prob / (1 - prior_prob)))
            sys.stdout.write('prior probability = {0}\n'.format(prior_prob))
            sys.stdout.write('Bayes factor = {0}\n'.format(bf))
            sys.stdout.write('2ln(Bayes factor) = {0}\n'.format(2 * math.log(bf)))
        sys.stdout.write('\n')

if __name__ == '__main__':
    main_cli()

