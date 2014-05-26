.. _installation:

************
Installation
************

You should have |python27|_ installed.

First, if you have |git|_ installed, clone the Git repository::

    $ git clone https://github.com/joaks1/PyMsBayes.git

If you do not have |git|_, you can download a snapshot of the repository as a
tar or zip archive by clicking the respective folder icons at the top of this
page. Either way, move into the downloaded directory::

    $ cd PyMsBayes

Next, let's make sure the bundled |dpp-msbayes|_ tools work on your machine::

    $ python setup.py test

This will run a small number of tests to ensure the tools bundled with |pmb|_
work on your computer. If these tests pass, you are good to go and can install
|pmb|_::

    $ sudo python setup.py install

If any of the tests fail and/or you get the following error message when you
try to install |pmb|_::

    **********************************************************************
    ****************************** WARNING *******************************
    The bundled `dpp-msbayes` tools are not being installed, because some
    of them are not executable on this system. The `PyMsBayes` package and
    scripts will still be installed, however, you will need to build and
    install `dpp-msbayes` (https://github.com/joaks1/dpp-msbayes) yourself
    in order to use them.  Sorry for the inconvenience.
    **********************************************************************

You will need to build and install |dpp-msbayes|_ from |dpp-msbayes-url|, and
install |abctb|_ from |abctb_url|. After you install |dpp-msbayes|_ and
|abctb|_ to your PATH, if you already ran the ``sudo python setup.py install``
command above, you are good to go.  If you haven't, ``cd`` back to the
``PyMsBayes`` directory and run this command (Note: you will still get the
above warning message, but you can ignore it because the |dpp-msbayes|_ tools
are already on your system).

If the install was successful, you should be able to call up the help menu of
the main program of |pmb|_::

    $ dmc.py -h

You should see output that begins like::

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
    
    dmc.py Version 0.2.4
    
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
    
