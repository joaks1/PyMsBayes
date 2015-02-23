.. role:: bolditalic
.. role:: hlight 
.. role:: codehlight 

.. _caveats:

***************************
Caveats and Recommendations
***************************


Be skeptical of strongly supported results
==========================================

Both |dpp-msbayes|_ and |msbayes|_ attempt to approximate the posterior of a
very parameter rich and highly stochastic model.
Furthermore, both methods use relatively little information from the data
when doing so (i.e., the sequence alignments are reduced to a small number
of simple statistics).
Under these conditions we should expect there to be quite a bit of posterior
uncertainty about divergence models and the number of divergence events.
A lack of posterior uncertainty (i.e., strongly supported results) is likely a
red flag that the prior on divergence times is strongly weighting the posterior
probabilities of the divergence models (i.e., the strong support is being
driven by the prior on the divergence times, and **not** the data).
See the :ref:`the section about the prior on divergence
times<prior_on_divergence_times>` for more details.


Issues of ABC model choice
==========================

ABC approximations of posterior probabilities (and associated Bayes factors)
can be inaccurate when the summary statistics are insufficient across models
(which will always be the case for |dpp-msbayes|_ and |msbayes|_ analyses)
:cite:`Robert2011`.
Thus, ABC model choice should really be treated as a means of exploring your
data, rather than a rigorous evaluation of hypotheses.


Assess the sensitivity of results to prior assumptions
======================================================

Given some of the pitfalls discussed above, it is always a good idea to
assess the sensitivity of your results to your prior assumptions.
Luckily, |pmb|_ makes this very easy for you, because it allows
you to analyze data under an arbitrary number of models very
easily.
See the :ref:`multi-model tutorial<multiple_model_analysis>` about how to do
this.


Simulation-based analyses can illuminate temporal resolution
============================================================

As shown in :cite:`Oaks2012`, :cite:`Oaks2014reply`, and :cite:`Oaks2014dpp`,
simulation-based analyses provide insight about the time scale at which the
method can detect random variation in divergence times from your data.
This can be very useful in guiding your interpretation of the results.


Do not scramble the summary statistics among taxa!
==================================================

See :ref:`the section about the re-sorting across taxa<sorting_taxa>` that is
done in |msbayes|_ (and is still an option in |dpp-msbayes|_).
This is mathematically incorrect, introduces bias, and should not be done.


The method might be inappropriate for comparing many taxa
=========================================================

See the note about :ref:`the number of divergence models as the number of taxa
increase<too_many_taxa>`.
Once the number of taxa gets above 15 or so, the number of possible models of
divergence gets enormous, and might be beyond the capacity of the simple
rejection sampling algorithm implemented in |dpp-msbayes|_ and |msbayes|_.
If you use either method with this many pairs of taxa, you should run multiple
replicates, each with large numbers of samples from the prior, to make sure
your estimates are stabilizing as the samples increase within each run, and
converging to similar values across runs.


Please choose your priors
=========================

Please do not use the default priors. You know your system infinitely better
than I do, and thus you should be the one to choose the assumptions about prior
uncertainty.

