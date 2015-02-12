.. role:: bolditalic
.. role:: hlight 
.. role:: codehlight 

.. _parameter_key:

****************************************************************
A key to the output of model parameters and associated statstics
****************************************************************

Below is a key to all of the parameter and statistics you will see in the
output.
In the ``posterior-sample.txt``, all column headers that start with ``PRI.``
designate a parameter or a summary statistic of parameters.
All other column headers designate summary statistics calculated from the
simulated data, which are used to calculate the distance from the summary
statistics calculated from the observed data.


.. contents:: 
    :local:
    :depth: 2

The parameters
==============


PRI.t.1, PRI.t.2, ...
---------------------

These are the divergence times of the populations pairs.
The pairs are numbered in the order that they appear (from top to bottom) in
the configuration file.

PRI.d1theta.1, PRId.1theta.2, ...
---------------------------------

The effective population size of one of the descendant populations of
each pair of populations.
The pairs are numbered in the order that they appear (from top to bottom) in
the configuration file.

PRI.d2theta.1, PRI.d2theta.2, ...
---------------------------------

The effective population size of the other descendant population of each pair
of populations.
The pairs are numbered in the order that they appear (from top to bottom) in
the configuration file.

PRI.atheta.1, PRI.atheta.2, ...
-------------------------------

The effective population size of the ancestral population of each pair of
populations.
The pairs are numbered in the order that they appear (from top to bottom) in
the configuration file.

PRI.div.model
-------------

The model of divergence (see :ref:`the section on divergence
models<comparative_divergence_models>`).
In the ``posterior-sample.txt`` file, only indices of the
divergence models are given, which is not very useful.
However, the ``div-model-results.txt`` file contains all of the information
regarding the posterior sample of divergence models.

PRI.model
---------

The prior model.
This will not be meaningful if only analyzing data under a single prior model.
However, if you average over multiple models, the estimated posterior
probability of each model is reported in the ``model-results.txt`` file.


The statistics
==============

:hlight:`NOTE:` the variables below are :hlight:`NOT` parameters of the model.
They are statistics that summarize the values of parameters.

PRI.Psi
-------

The number of divergence events, or the number of divergence-time parameters in
the divergence model.

PRI.E.t
-------

The mean of the divergence times across all pairs of populations.
:hlight:`NOTE`, this is the mean across all the population pairs, not the
divergence-time parameters (which can be shared across pairs).

PRI.var.t
---------

The variance of the divergence times across all pairs of populations.
:hlight:`NOTE`, this is the variance across all the population pairs, not
the divergence-time parameters (which can be shared across pairs).

.. _pri_omega:

PRI.omega
---------

The dispersion index (i.e., variance / mean) of the divergence times across all
pairs of populations.
:hlight:`NOTE`, this is the dispersion index across all the population pairs,
not the divergence-time parameters (which can be shared across pairs).

PRI.cv
------

The coefficient of variation (CV; i.e., standard deviation / mean) of the
divergence times across all pairs of populations.
:hlight:`NOTE`, this is the CV across all the population pairs, not the
divergence-time parameters (which can be shared across pairs).

