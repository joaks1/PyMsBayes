.. role:: bolditalic
.. role:: hlight 
.. role:: codehlight 

.. contents:: 
    :local:
    :depth: 3

.. _multiple_model_analysis:

*******************************************
An Empirical Analysis under Multiple Models
*******************************************

For this tutorial, I will assume that you have already run through
the :ref:`simple empirical analysis<simple_empirical_analysis>`,
so we will go into a lot less detail here.

In the tutorial for the :ref:`simple empirical
analysis<simple_empirical_analysis>`,
we performed a basic analysis of the lizard sequence data under
the model specified in the ``dpp-simple.cfg`` configuration file.
We did this by entering the following at the command line
from within the |lizard-config-dir| directory:

.. parsed-literal::

    $ |dmc| -o dpp-simple.cfg -p dpp-simple.cfg -n 5000

One nice feature of |ldmc| is that it is very easy to perform analyses
of the same data under an arbitrary number of different models.
For example, we can easily expand on our example above to analyze the lizard
data under two different models: the |dpp-msbayes|_ model specified in the
|lizard-config-dir|\ ``/dpp-simple.cfg`` configuration file and an |msbayes|_
model specified in the |lizard-config-dir|\ ``msbayes.cfg`` file:

.. parsed-literal::

    $ |dmc| -o dpp-simple.cfg -p dpp-simple.cfg msbayes.cfg -n 5000

Notice, we now specify both of the model configuration files for the prior
model (``-p``) option.
Also, because both of the configuration files share the same sample table,
it does not matter which one is specified for the observed-data configuration
option (``-o``) (i.e., either way the observed summary statistics will be
calculated from the same sequence alignments).

What this new command accomplishes is that |ldmc| will now analyze the
data under both models simultaneously.
This will take a little longer to run than the previous analysis, but should
finish in less than a few minutes or so on a modern laptop.

The output
==========

Because we analyzed the data under two models, we will now have
more output.
There should be a new directory named |result-dir| that was created when you
ran the analysis.
If you still had the |result-dir| in the |lizard-config-dir| from the previous
analysis, then the new results will be output to a new directory named
|result-dir|\ ``-0``.

The info file
-------------

Our new |info-path|, will look something
like::

    [pymsbayes]
    	version = Version 0.2.4
    	output_directory = /home/jamie/software/dev/pymsbayes-tutorial-data/lizards/configs/pymsbayes-results-0
    	temp_directory = /home/jamie/software/dev/pymsbayes-tutorial-data/lizards/configs/pymsbayes-results-0/temp-files-5ym88h
    	sort_index = 0
    	simulation_reps = 0
    	seed = 192883961
    	num_processors = 8
    	num_prior_samples = 5000
    	num_standardizing_samples = 5000
    	bandwidth = 0.002
    	posterior_quantiles = 1000
    	posterior_sample_size = 1000
    	stat_patterns = ^\s*pi\.\d+\s*$, ^\s*wattTheta\.\d+\s*$, ^\s*pi\.net\.\d+\s*$, ^\s*tajD\.denom\.\d+\s*$
    	num_taxon_pairs = 3
    	dry_run = False
    	[[tool_paths]]
    		dpp_msbayes = /home/jamie/software/dev/PyMsBayes/bin/linux/dpp-msbayes.pl
    		msbayes = /home/jamie/software/dev/PyMsBayes/bin/linux/msbayes.pl
    		eureject = /home/jamie/software/dev/PyMsBayes/bin/linux/eureject
    		abcestimator = /home/jamie/software/dev/PyMsBayes/bin/linux/ABCestimator
    	[[observed_configs]]
    		1 = ../dpp-simple.cfg
    	[[observed_paths]]
    		1 = observed-summary-stats/observed-1.txt
    	[[prior_configs]]
    		1 = ../dpp-simple.cfg
    		2 = ../msbayes.cfg
    	[[run_stats]]
    		start_time = 2015-02-11 19:40:17.874783
    		stop_time = 2015-02-11 19:41:33.645795
    		total_duration = 0:01:15.771012

This is very similar to our first analysis, except that there are now two prior
configs listed.


The data and model keys and corresponding directories
-----------------------------------------------------

If we look within the new |output-dir| directory, we will 
see that the ``model-key.txt`` now contains::

    m1 = ../../dpp-simple.cfg
    m2 = ../../msbayes.cfg

This tells us that we will find the results of the analysis under
the ``dpp-simple.cfg`` model within the |output-dir|\ ``/d1/m1``
directory, and will find the results of the analysis under
the ``msbayes.cfg`` model within the |output-dir|\ ``/d1/m2`` directory.

Model-averaged results
----------------------

In addition to the results of the analyses under both models, there is also
results in the |output-dir|\ ``/d1/m12-combined`` directory that are averaged
over both models.
This kind of model averaged results were encouraged by :cite:`Hickerson2013`,
however, :cite:`Oaks2014reply` found such model-averaged analyses to perform
quite poorly.
Nonetheless, the option is yours to explore these types of model-averaging
analyses.

Rather than summarize and plot the results of this very short example
multi-model analysis, let's go ahead and look at some provided output from the
same analysis run for many more samples from the prior models.

Output from a longer analysis
=============================

In the tutorial data you downloaded, you can find the output of the same
multi-model analysis we ran above, but run for much longer.
You can find these results in the directory, |lizard-example-results|_.

These results were generated via the following command, which took
less than a day to run on a laptop::

    $ dmc.py --np 8 \
        -o dpp-simple.cfg \
        -p dpp-simple.cfg msbayes.cfg \
        -n 10000000 \
        --prior-batch-size 12500 \
        --num-posterior-samples 1000 \
        --num-standardizing-samples 100000 \
        -q 1000 \
        --reporting-frequency 200000 \
        --compress \
        --seed 845225390

You can see the full bash script that was used to run this analysis
in the file |lizard-example-qsub|_.
We can see that this analysis will draw 10 million samples from
both prior models (``-n 10000000``), and report the results
every 200,000 samples (``--reporting-frequency 200000``).
Furthermore, the first 100,000 samples from each prior will be used to
standardize the summary statistics calculated from all the prior simulations
and the observed data (``num-standardizing-samples 100000``).
I chose the number of processes (``--np 8``) and the size of the "batches" in
which the prior samples will be generated (``--prior-batch-size 12500``) such
that their product (100,000 samples) is a common factor of the number of
standardizing samples (100,000), the reporting frequency (200,000), and the
total number of samples (10 million).
This is not mandatory, but it will maximize the efficiency of the
multi-processing.
I also specify a seed for the random number generator (``--seed 845225390``) so
that I, or anyone else, can replicate the results.
The seed is always reported in the ``pymsbayes-info.txt`` file, so this is not
necessary unless you are trying to replicate previous results.
But, it is **always** good to report the seed (or at least make it available)
so that others can replicate your work.

Additional results reported during analysis
-------------------------------------------

If we look in the |lizard-example-results|\ ``/``\ |result-dir| directory, we
will find that most of the output (e.g., the info file, key files, directory
structure) is very similar to the short analysis we ran above.
The main difference is in the number of result files that are in the 
|lizard-example-results|\ ``/``\ |output-dir|\ ``/d1/m1``
,
|lizard-example-results|\ ``/``\ |output-dir|\ ``/d1/m2``
, and 
|lizard-example-results|\ ``/``\ |output-dir|\ ``/d1/m12-combined``
directories.
As you can see, we have results files for every 200,000 prior samples.

The trace file
--------------

There is also an additional "trace" file:

*   ``d1-m1-s1-trace.txt``

Which shows the posterior means of various estimates as the prior samples
accumulated.
This file is created whenever the ``--reporting-frequency`` option is used to
report results as the analysis progresses.
This file is useful for checking to see if estimates stabilized as the number
of prior samples evaluated increased, and to see if estimates from multiple,
independent analyses "converged" to similar estimates.

Plotting the results
--------------------

...

