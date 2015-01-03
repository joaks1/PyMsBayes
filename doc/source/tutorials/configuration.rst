.. role:: bolditalic
.. role:: hlight 

.. _config:

**********************
The Configuration File
**********************

Configuration files are used to control analysis settings for |dpp-msbayes|_
and |msbayes|_.
|pmb|_ allows you to use |dpp-msbayes|_ and |msbayes|_ configuration files (and
their respective models) interchangeably.
Here is an example of a configuration file for |dpp-msbayes|_::

    concentrationShape = 1000.0
    concentrationScale = 0.00437
    thetaShape = 4.0
    thetaScale = 0.001
    ancestralThetaShape = 0
    ancestralThetaScale = 0
    thetaParameters = 000
    tauShape = 1.0
    tauScale = 0.02
    timeInSubsPerSite = 1
    bottleProportionShapeA = 0
    bottleProportionShapeB = 0
    bottleProportionShared = 0
    migrationShape = 0
    migrationScale = 0
    numTauClasses = 0
    
    BEGIN SAMPLE_TBL
    species-1	locus-1	1.0	1.0	10	8	32.42	389	0.27	0.24	0.26	species-1-locus-1.fasta
    species-1	locus-2	1.0	1.0	8	6	5.51	500	0.25	0.22	0.24	species-1-locus-2.fasta
    species-1	locus-3	1.0	1.0	6	8	8.38	524	0.26	0.23	0.26	species-1-locus-3.fasta
    species-1	locus-4	1.0	1.0	8	10	5.20	345	0.25	0.23	0.24	species-1-locus-4.fasta
    species-1	locus-5	1.0	1.0	8	8	29.59	417	0.27	0.23	0.21	species-1-locus-5.fasta
    species-1	locus-mt	0.25	4.0	5	5	8.15	600	0.22	0.24	0.27	species-1-locus-mt.fasta
    species-2	locus-1	1.0	1.0	6	10	7.53	400	0.25	0.24	0.26	species-2-locus-1.fasta
    species-2	locus-3	1.0	1.0	10	8	11.14	550	0.27	0.22	0.24	species-2-locus-3.fasta
    species-2	locus-4	1.0	1.0	8	8	9.39	350	0.24	0.24	0.23	species-2-locus-4.fasta
    species-2	locus-5	1.0	1.0	10	10	13.32	450	0.26	0.24	0.22	species-2-locus-5.fasta
    species-2	locus-mt	0.25	4.0	4	5	7.59	549	0.23	0.26	0.23	species-2-locus-mt.fasta
    species-3	locus-1	1.0	1.0	10	6	17.03	367	0.25	0.23	0.27	species-3-locus-1.fasta
    species-3	locus-3	1.0	1.0	8	10	59.17	541	0.26	0.22	0.25	species-3-locus-3.fasta
    species-3	locus-4	1.0	1.0	6	8	6.90	333	0.28	0.23	0.21	species-3-locus-4.fasta
    species-3	locus-mt	0.25	4.0	5	4	11.42	587	0.22	0.22	0.25	species-3-locus-mt.fasta
    END SAMPLE_TBL

A configuration file has two parts:

#. The preamble_ containing keyword arguments for specifying prior probability
   distributions of model parameters.
#. The sample_table_ containing information and parameter values for each sequence
   alignment. This section is delimited by the ``BEGIN SAMPLE_TBL`` and ``END
   SAMPLE_TBL`` lines.


.. _preamble:

The Preamble
============

There are two main differences to the preamble between |pmb|_ and
|msbayes|_.

#. The |dpp-msbayes|_ preamble gives you more (and different) options to
   control the number of parameters in the model and the priors on those
   parameters.
#. The |msbayes|_- and |dpp-msbayes|_-style preambles are processed differently
   in |pmb|_: The keywords in the preamble are case insensitive, and if any
   unrecognized keywords are encountered, an error is reported and the process
   exits (i.e., crashes).
   
   In |msbayes|_, the keywords are case-sensitive and, more importantly, the
   the default value for a setting is :hlight:`quietly` used when any options
   are mis-typed.
   For example, if you specify ``uppertheta = 0.01`` in an |msbayes|_ config
   file (note the lower case "t"), it will quietly use the default setting for
   ``upperTheta`` and :hlight:`not` report any warning or error.
   So, if you prefer the model implemented in |msbayes|_, |pmb|_ is useful as a
   "safer" interface for it.

Below we walk through all of the preamble options for the |dpp-msbayes|_
preamble.
If you want to use the |msbayes|_ model with |pmb|_, please see the
`msBayes documentation
<https://docs.google.com/document/d/15heQlz60cGe6GWKcXqf1AIYMBZ6p2sKuGct-EoEHiNU/edit>`_
for information about the preamble.


.. contents:: Preamble Elements
    :local:


.. _concentration_parameter_prior:

``concentrationShape`` / ``concentrationScale``
-----------------------------------------------
    
The shape and scale parameters of a gamma-distributed prior on the
concentration parameter of the Dirichlet-process prior on divergence models.


.. _theta_parameterization:

``thetaParameters``
-------------------

|dpp-msbayes|_ gives you full control over the parameterization of the
population sizes for each pair of populations, which we will refer to as
:math:`\ancestralTheta{}`, :math:`\descendantTheta{1}{}`, and
:math:`\descendantTheta{2}{}` for the ancestral, and two descendant
populations.

This setting is controlled by a sequence of three integers that always starts
with ``0`` and increments by 1 whenever a free parameter is added.
The first two integers represent :math:`\descendantTheta{1}{}` and
:math:`\descendantTheta{2}{}`, and the last integer represents
:math:`\ancestralTheta{}`.

``000`` is one extreme, where :math:`\ancestralTheta{}`,
:math:`\descendantTheta{1}{}`, and :math:`\descendantTheta{2}{}` are all
constrained to be equal for each population pair (population sizes will still
vary among pairs).

``012`` is the other extreme, where :math:`\ancestralTheta{}`,
:math:`\descendantTheta{1}{}`, and :math:`\descendantTheta{2}{}` are all
estimated as independent parameters for each pair of populations.
This is most similar to the original |msbayes|_
However, the descendant population sizes are constrained to be negatively
correlated in |msbayes|_ (see :cite:`Oaks2012` and :cite:`Oaks2014dpp`).

Another example is ``001``: the descendant populations share the same size
parameter, but the ancestral population size is free to vary.

For ``011`` and ``010``, one of the descendant population is constrained to the
same size as the ancestral, and the other is free to vary.


.. _theta_prior:

``thetaShape`` / ``thetaScale``
-------------------------------

These settings define the shape and scale parameters of a gamma prior on the
effective population sizes. Population sizes are scaled by the per-site
mutation rate (:math:`\mu`): :math:`4N_e\mu`.


.. _ancestral_theta_prior:

``ancestralThetaShape`` / ``ancestralThetaScale``
-------------------------------------------------

If these settings are both provided and both are positive, they define the
shape and scale parameters of a gamma prior on the effective population size of
ancestral populations.

If they are excluded, or both are zero, the ``thetaShape`` and ``thetaScale``
settings are used for the gamma prior on ancestral population size parameters.

.. _divergence_time_prior:

``tauShape`` / ``tauScale``
---------------------------

These settings define the shape and scale parameters of a gamma prior
on divergence times. 
See the timescale_setting_ ``timeInSubsPerSite`` setting for the information on
the time units.


.. _timescale_setting:

``timeInSubsPerSite``
---------------------

This setting controls the time scale of the model and has two settings:

* ``timeInSubsPerSite = 1``: Time units are in expected substitutions per site.
  For example, a divergence of 0.05 means that, on average, 5% of sites have
  changed since the populations diverged (so you expect 10% divergence between
  the populations since the population divergence).
  Thus, you can convert these units to the number of generations by dividing by
  the mutation rate.

* ``timeInSubsPerSite = 0``: Time units are in coalescent units,
  :math:`\globalcoalunit` generations, where :math:`\globalpopsize` is the size of a
  constant reference population based on the mean of the theta_prior_
  (defined by settings ``thetaShape`` and ``thetaScale``).

  If we use :math:`\globaltheta` to represent the mean of the theta prior, then
  
  .. math::
      \globalcoalunit = \frac{\globaltheta}{\mutationRate},

  where :math:`\mutationRate` is the per-site mutation rate.  Thus, you can
  convert these ":math:`\globalcoalunit` generations" units to the number of
  generations by assuming a mutation rate and multiplying by
  :math:`(\globaltheta/\mutationRate)`.  See :cite:`Oaks2014dpp` for more
  details.

  .. note::
      :hlight:`Why use the mean of the prior on theta to scale time?`
      I have no idea.
      This is legacy from |msbayes|_, and is the default setting.
      However, I strongly discourage using this time scale, because it makes it
      very difficult to compare results across analyses with different settings
      for the theta_prior_.
      It also requires you to re-scale the divergence_time_prior_ every time
      you change the theta_prior_.
      Scaling time by the expected substitutions per site is much more straight
      forward.


``bottleProportionShapeA`` / ``bottleProportionShapeB``
-------------------------------------------------------

If both are positive, these settings define the shape parameters alpha and
beta, respectively, of a beta prior on the magnitude of a post-divergence
bottleneck in each of the descendant populations.

The bottleneck magnitude is the proportion of the effective population size
that remains following the bottleneck.
For example, a value of 0.95 would mean that bottleneck reduces the effective
population size by 5%.

If either or both are zero or less, there is no post-divergence population
bottleneck in the descendant populations (i.e., the bottleneck-magnitude
parameters, along with the timing of each bottleneck, are removed from the
model).

.. note::
    There are also parameters in the model for the timing of the end of the
    bottleneck (it begins at speciation in forward time). There is one of these
    parameters for each pair of populations (i.e., the descendant populations
    of each pair share the same bottleneck-end-time parameter).
    Thus if either or both of the
    ``bottleProportionShapeA``/``bottleProportionShapeB`` settings are zero or
    less, you are also removing these bottleneck timing parameters from the
    model.
    This means you are removing :math:`3\npairs{}` parameters from the model,
    where ":math:`\npairs{}`" is the number of pairs of populations.


``bottleProportionShared``
--------------------------

If ``bottleProportionShared = 0``, then there are two free bottleneck-magnitude
parameters for each population pair (one for each descendant population).
If ``bottleProportionShared = 1``, then there is one bottleneck-magnitude
parameter for each population pair (i.e., the descendant populations of each
pair share the same bottleneck magnitude; the bottleneck magnitude still varies
among the pairs).

.. note::
    This setting is overridden if either or both of the
    ``bottleProportionShapeA`` and ``bottleProportionShapeB`` settings is zero
    or less (because then there is no bottleneck at all).


``migrationShape`` / ``migrationScale``
---------------------------------------

These settings define the shape and scale parameters of a gamma prior on the
symmetric migration between the descendant populations of each pair (in units
of the number of gene copies per generation).

If either or both settings are zero or less, there is no migration in the
model.


``recombinationShape`` / ``recombinationScale``
-----------------------------------------------

These settings define the shape and scale parameters of a gamma prior
on the intragenic recombination rate.

If either or both are zero or less, there is no recombination in the model.

.. note::
    I recommend :hlight:`not` including intragenic recombination in the model,
    because the current implementation is very inefficient and poorly tested.


``numTauClasses``
-----------------

If this setting is zero (the default), the number of divergence events is free
to vary according to the Dirichlet process prior on divergence models.

If it is greater than zero, then the model is constrained to ``numTauClasses``
divergence events.
This is useful for simulation-based power analyses, but should not be used for
empirical analyses.


.. _sample_table:

The Sample Table
================

Blah

