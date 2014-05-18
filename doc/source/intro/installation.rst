.. _installation:

Installation
============

You should have |python27|_ and |git|_ installed.

First, clone the Git repository::

    $ git clone git://github.com/joaks1/PyMsBayes.git

Next, let's make sure |pmb|_ works on your maching::

    $ cd PyMsBayes
    $ python setup.py test

This will run a small number of tests to ensure the tools bundled with |pmb|_
work on your computer. If these tests pass, you are good to go and can install
|pmb|_::

    $ sudo python setup.py develop

If the install was successful, you should be able to call up the help menu of
the main program of |pmb|_::

    $ dmc.py -h

You should see output that looks like::

    usage: dmc.py [-h] -o OBSERVED_CONFIGS [OBSERVED_CONFIGS ...] -p PRIOR_CONFIGS
                  [PRIOR_CONFIGS ...] [-r REPS] [-n NUM_PRIOR_SAMPLES]
                  [--prior-batch-size PRIOR_BATCH_SIZE] [--generate-samples-only]
                  [--num-posterior-samples NUM_POSTERIOR_SAMPLES]
                  [--num-standardizing-samples NUM_STANDARDIZING_SAMPLES]
                  [--np NP] [--output-dir OUTPUT_DIR] [--temp-dir TEMP_DIR]
                  [--staging-dir STAGING_DIR]
                  [-s [STAT_PREFIXES [STAT_PREFIXES ...]]] [-b BANDWIDTH]
                  [-q NUM_POSTERIOR_QUANTILES]
                  [--reporting-frequency REPORTING_FREQUENCY]
                  [--sort-index {0,1,2,3,4,5,6,7,8,9,10,11}]
                  [--no-global-estimate] [--compress] [--keep-temps] [--seed SEED]
                  [--output-prefix OUTPUT_PREFIX] [--data-key-path DATA_KEY_PATH]
                  [--start-from-simulation-index START_FROM_SIMULATION_INDEX]
                  [--start-from-observed-index START_FROM_OBSERVED_INDEX]
                  [--dry-run] [--version] [--quiet] [--debug]
    
    dmc.py Version 0.2.0
    
    optional arguments:
      -h, --help            show this help message and exit
      -o OBSERVED_CONFIGS [OBSERVED_CONFIGS ...], --observed-configs OBSERVED_CONFIGS [OBSERVED_CONFIGS ...]
                            One or more msBayes config files to be used to either
                            calculate or simulate observed summary statistics. If
                            used in combination with `-r` each config will be used
                            to simulate pseudo-observed data. If analyzing real
                            data, do not use the `-r` option, and the fasta files
                            specified within the config must exist and contain the
                            sequence data.
    
