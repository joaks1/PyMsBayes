.. _background:

**********
Background
**********

.. _comparative_divergence_models:

Comparative divergence models
=============================

Biogeographers often seek to explain diversity on historical events.
That is, using genetic data from contemporary populations, we would like to
infer diversification patterns and see if they support patterns predicted by
past event(s) of interest (e.g., islands fragmented by rises in sea level,
changes in climate fragmenting communities into refugia, etc.).
For example, if an event split a community of species 260,000 years ago, we
might expect the divergences across multiple species co-distributed across the
barrier created by the event (the ominous "black rectange" below) to be
temporally clustered.
More specifically, let's say we are interested in investigating three species
of lizards that are co-distributed across the putative barrier.
In order to infer the affect of the historical event on diversification, we
want to compare, across the three species, the timing of the divergence between
the populations on oppositie sides of the putative barrier.
If the historical event caused divergence, we would expect that each of the
three pairs of lizard populations (or some subset of them) diverged at the same
time, as shown in the divergence_model_111_ figure.

.. _divergence_model_111:
.. figure:: /_static/div-model-cartoon-111.png
   :align: center
   :width: 600 px
   :figwidth: 60 %
   :alt: divergence model 111
   
   A cartoon showing three pairs of lizard populations that co-diverge due to
   an event 260,000 years ago.

We can think of this as a particular *divergence model* where all three pairs
of populations share the same divergence-time parameter.
If we give the divergence-time parameter the index "1", we can use the notation
:math:`\divModel{1} = 111` to show that this divergence model assigns
population pairs 1, 2, and 3 to divergence-time parameter 1.
However, this is only one possible divergence model, and happens to be the most
constrained.
With three population pairs, there are 4 other possible models of divergence (5
total possible models).
Three of these models have two divergence-time parameters.
We can assign population-pair 1 to a second divergence-time parameter to get
divergence model :math:`\divModel{2} = 211`, as shown in the divergence_model_211_
figure.

.. _divergence_model_211:
.. figure:: /_static/div-model-cartoon-211.png
   :align: center
   :width: 600 px
   :figwidth: 60 %
   :alt: divergence model 211
   
   A cartoon showing population-pair 1 assigned to divergence-time parameter 2,
   and population-pairs 2 and 3 assigned to divergence-time parameter 1.

We can also assign population-pair 2 to divergence-time parameter 2 to get
divergence model :math:`\divModel{3} = 121`, as shown in the divergence_model_121_
figure.

.. _divergence_model_121:
.. figure:: /_static/div-model-cartoon-121.png
   :align: center
   :width: 600 px
   :figwidth: 60 %
   :alt: divergence model 121
   
   A cartoon showing population-pair 2 assigned to divergence-time parameter 2,
   and population-pairs 1 and 3 assigned to divergence-time parameter 1.

And for the last possible divergence model with two divergence-time parameters,
we assign population-pair 3 to divergence-time parameter 2 to get divergence
model :math:`\divModel{4} = 112`, as shown in the divergence_model_112_ figure.

.. _divergence_model_112:
.. figure:: /_static/div-model-cartoon-112.png
   :align: center
   :width: 600 px
   :figwidth: 60 %
   :alt: divergence model 112
   
   A cartoon showing population-pair 3 assigned to divergence-time parameter 2,
   and population-pairs 1 and 2 assigned to divergence-time parameter 1.

Finally, we can add a third divergence-time parameter so that each pair of
populations is assigned to its own divergence-time parameter (divergence model
:math:`\divModel{5} = 123`), as shown in the divergence_model_123_ figure.
This is the most general model of divergence, and has no co-divergence among
taxa.
Biogeographically, we can think of each free divergence-time parameter
as a "divergence event" during which one or more pairs of populations
can diverge.

.. _divergence_model_123:
.. figure:: /_static/div-model-cartoon-123.png
   :align: center
   :width: 600 px
   :figwidth: 60 %
   :alt: divergence model 123
   
   A cartoon showing the most general model of divergence where all three
   pairs of lizard populations diverge at unique times.

Being energetic herpetologists, we go out and sample individuals from each of
the lizard populations, and from those individuals collect DNA sequence data
from one or more orthologous loci per pair of populations.
You can find our example sequence data in fasta format in the |lizard-seq-dir|_
directory.
We know that the sequences of a locus are related by a genealogy,
and that the shape of this genealogy is governed by demographic processes.
We also know that the genetic variation we see in the data accumulated as the
sequences evolved via mutational processes along the genealogy.
We can modify our cartoon of model :math:`M_5 = 123` to better
represent this, as shown in figure pop_divergence_model_123_.

.. _pop_divergence_model_123:
.. figure:: /_static/pop-div-model-cartoon-123.png
   :align: center
   :width: 600 px
   :figwidth: 60 %
   :alt: divergence model 123
   
   A cartoon showing the most general model of divergence where all three
   pairs of lizard populations diverge at unique times.

Next, let's jump to the ":ref:`bayesian_divergence_model_choice`" section to
see how we can use the information in the sequence data to infer the temporal
distribution of the population divergences across the three lizard species.

.. _bayesian_divergence_model_choice:

Bayesian divergence-model choice
================================

In the figures above, we used :math:`\divTimeMap{1}, \divTimeMap{2},` and
:math:`\divTimeMap{3}` to represent the divergence times of the three
pairs of lizard populations. Now, let's use :math:`\divTimeMapVector`
to represent all three divergence times; that is, 
:math:`\divTimeMapVector = \divTimeMap{1}, \divTimeMap{2}, \divTimeMap{3}`.
The number of unique divergence-time values (i.e., the number of free
divergence-time parameters) within :math:`\divTimeMapVector`, and the
assignment of the lizard species to these values, depends on the divergence
model.
For example, for model :math:`\divModel{1}` in Figure divergence_model_111_
above, the divergence times would be 
:math:`\divTimeMapVector = 260, 260, 260`
(in thousands of years).
For model :math:`\divModel{5}` in Figure divergence_model_123_ above, the
divergence times would be 
:math:`\divTimeMapVector = 260, 96, 397`.
In order to learn about the affect the "black rectangle" had on the
diversification of these lizard populations, it would be ideal if we could
jointly infer the divergence model and the divergence times from the DNA
sequence data we collected.

In order to do this, we need to assume a probabilistic evolutionary model
that gave rise to the data we collected.
If we assume a Markov chain model of nucleotide substitution, we can calculate
the probability of the sequence data given the genealogies and a set of
parameter values for the substitution model.
Both |dpp-msbayes|_ and |msbayes|_ assume an HKY85 model of nucleotide
substitution :cite:`HKY`.
If we further assume a coalescent model of ancestral processes, we can
calculate the probability of the genealogies given the parameter values
for the sizes of the populations.
For simplicity, let's lump all the parameters of the substitution and
coalescent models for all three pairs of lizard populations into
:math:`\demographicParamVector`.
Let's also use :math:`\alignmentVector` to represent all of our sequence
alignments (which are in the directory |lizard-seq-dir|_).
Lastly, let's use :math:`\geneTreeVector` to represent all of the gene trees
(one for each alignment) that relate the sequences in our alignments.
If we make assumptions about the relative rates of mutations and the relative
generation times among the three lizard species, we can calculate the posterior
probability distribution of the divergence times (and other nuisance
parameters) given the data and one of the models of divergence using Bayes
rule:

.. math::
    :label: postdensity

    p(\divTimeMapVector, \geneTreeVector, \demographicParamVector \given
    \alignmentVector, \divModel{1}) = \frac{p(\alignmentVector \given
    \divTimeMapVector, \geneTreeVector, \demographicParamVector,
    \divModel{1})p(\divTimeMapVector,\geneTreeVector,\demographicParamVector
    \given \divModel{1})}{p(\alignmentVector\given\divModel{1})}

The denominator of Bayes' rule (Equation :eq:`postdensity`) is the marginal
probability of the data under divergence model :math:`\divModel{1}`, a.k.a the
marginal likelihood of divergence model :math:`\divModel{1}`.
This is equal to the integral over the entire parameter space of model
:math:`\divModel{1}` of the likelihood density weighted by the prior density:

.. math::
    :label: marginallike

    p(\alignmentVector \given \divModel{1}) =
    \int_{\divTimeMapVector}
    \int_{\geneTreeVector}
    \int_{\demographicParamVector}
    p(\alignmentVector \given \divTimeMapVector, \geneTreeVector,
    \demographicParamVector, \divModel{1})
    p(\divTimeMapVector, \geneTreeVector, \demographicParamVector, \given
    \divModel{1})
    d\divTimeMapVector
    d\geneTreeVector
    d\demographicParamVector

You can think of this as the "average" likelihood of divergence model
:math:`\divModel{1}`, and this average is weighted by the prior over the
entire space of the model.
If we calculate the marginal likelihood of all five possible divergence
models, we can use Bayes' rule again to calculate the posterior probability
of divergence model :math:`\divModel{1}` given our sequence data:

.. math::
    :label: postmass1

    p(\divModel{1} \given \alignmentVector) = \frac{ p(\alignmentVector \given
    \divModel{1}) p(\divModel{1}) }{
    p(\alignmentVector \given \divModel{1}) p(\divModel{1}) +
    p(\alignmentVector \given \divModel{2}) p(\divModel{2}) +
    p(\alignmentVector \given \divModel{3}) p(\divModel{3}) +
    p(\alignmentVector \given \divModel{4}) p(\divModel{4}) +
    p(\alignmentVector \given \divModel{5}) p(\divModel{5}) }

Or, more generally, we can calculate the posterior probability of any
divergence model ":math:`i`" using:

.. math::
    :label: postmass

    p(\divModel{i} \given \alignmentVector) = \frac{ p(\alignmentVector \given
    \divModel{i}) p(\divModel{i}) }{ \sum_{i} p(\alignmentVector \given
    \divModel{i}) p(\divModel{i}) }

This is essentially the relative marginal likelihood of the model (it is
exactly that if assume equal prior mass for each divergence model).
We can combine Equations :eq:`postdensity` and :eq:`postmass` to better
represent that we will be jointly inferring the posterior probabilities of
divergence models and the posterior densities of the divergence models'
parameters:

.. math::
    :label: jointpost

    p(\divTimeMapVector, \geneTreeVector, \demographicParamVector, \divModel{i}
    \given \alignmentVector) = \frac{p(\alignmentVector \given
    \divTimeMapVector, \geneTreeVector, \demographicParamVector,
    \divModel{i})p(\divTimeMapVector,\geneTreeVector,\demographicParamVector
    \given \divModel{i})p(\divModel{i})}{p(\alignmentVector)}

By jointly sampling over the posterior of all the divergence models, Equation
:eq:`jointpost` will also give us model-averaged estimates of the divergence
times for each of our pairs of populations (i.e., we get estimates of
divergence times that account for uncertainty in divergence models).

The key take home here is that the *marginal* likelihoods are the "guts" of
Bayesian model choice, as shown in Equation :eq:`postmass`.
I.e., it is the *marginal* probability of our data under a given model that
updates our prior expectation and informs the posterior probability of that
model.
As you might expect, because the marginal likelihoods are weighted by the
priors on parameters, the posterior probabilities of the models can be quite
sensitive to these priors.
NOTE, it is important to realize here that the posterior probability of the
models can be very sensitive to the priors on the *parameters*, not just the
priors on the *models* themselves.
Thus, we have to choose the priors on parameters carefully, and should always
assess the sensitivity of our results to differences in these prior
assumptions.
We will discuss how the choice of prior distribution on divergence times can
have a major affect on posterior probabilities of divergence models for both
|dpp-msbayes|_ and |msbayes|_ in the ":ref:`prior_on_divergence_times`"
section.
But first, let's talk about how we will approximate the posterior in Equation
:eq:`jointpost`.

We cannot calculate all of the integrals in Equation :eq:`marginallike`
exactly, so we will need to use a numerical integration algorithm to
approximate the posterior.
Furthermore, to avoid deriving and calculating the likelihood function, we will
use approximate likelihoods for our numerical integration algorithm.
(Digression: this is why I do not like the term "approximate Bayesian
computation." This describes *all* Bayesian applications except for trivial
models where the posterior can be solved exactly. "Approximate-likelihood
Bayesian computation" would be much more useful, but then we would lose the
beloved acronym ABC.)

Approximate-*likelihood* Bayesian computation
=============================================

We will use a simple Monte Carlo rejection algorithm based on approximate
likelihoods to approximate the posterior in Equation :eq:`postmass`.
Approximate-likelihood techniques use simulations to avoid calculating the
likelihood function.
The idea is very simple: given values for all the parameters in the model, we
simulate a dataset with the same "dimensions" as the observed data (i.e., the
same number of sequence alignments with same number of rows and columns), and
compare the simulated dataset to the observed data.
The closer to the observed data, the higher the likelihood for the set of
parameter values.
If we did this many times, randomly drawing the set of parameter values from
the prior distribution each time, and only retained the samples that produced
datasets that matched our observed sequence alignments (or sufficient summary
statistics of those alignments) exactly, this would be equivalent to an
exact-likelihood Bayesian approach.
However, the sun would probably burn out while we waited to run enough
simulations to collect a decent number of posterior samples in this way.
So, to make things more computationally tractable, we will introduce
two sources of approximation:

#. We will reduce our observed and simulated datasets down to a set of
   insufficient statistics. This adds a "fudge" factor to the method, because
   we are throwing away information in our data when we do this.
#. We will retain simulations that produce values of these insufficient
   statistics that are "close enough" to the values calculated from our
   observed data. This "wiggle room" (tolerance) around the observed summary
   statistics is another source of approximation.

For illustration purposes, let's assume we reduce our dataset for the three
pairs of lizard populations into one summary statistic per species; perhaps its
the average sequence divergence between the two populations.
Then, we will simulate lots of datasets under the model (each time based on a
set of parameter values drawn from the prior distribution) and reduce each of
them to the same three summary statistics.
Lastly, we retain the sets of parameter values that produced summary statistics
that fall within the "good enough" zone around our observed data.
An example of this is animated in the rejection_sampling_ gif below.

.. _rejection_sampling:
.. figure:: /_static/rejection-sampling.gif
   :align: center
   :width: 600 px
   :figwidth: 60 %
   :alt: rejection sampler gif
   
   An illustration of a Monte Carlo rejection sampler.

This animation begins with a blue dot representing the values of the three
summary statistics calculated from the observed sequence alignments.
Next, a grey sphere illustrates the "good enough" zone.
Then, we see black points accumulate, which represent the values of the three
summary statistics calculated from datasets that were simulated under sets of
parameter values drawn randomly from the prior.
Lastly, we see the retained sample of points that fell within our "good enough"
zone; this is our sample from the approximate posterior.

So, how do we decide how large the "good enough" zone is? Well, the smaller the
better, but this is governed arbitrarily by computational limitations.
What we do is simulate as many datasets as we are willing to wait for and then
select the desired number of them that produced summary statistics closest to
the observed summary statistics.
For example, we might draw 10 million sets of parameter values from the prior,
and keep the 10,000 sets that produced summary statistics nearest to the
observed statistics; that's our approximate posterior sample.
In this example, the radius of the "good enough" space is determined by the
distance between the observed summary statistics and the 10,001st nearest
simulated summary statistics.
Again, this is arbitrary; drawing 100 million samples and keeping the closest
10,000 would be better.
However, we can get a sense of whether we have evaluated a sufficient number of
samples from the prior by:

#. Keeping track of the parameter estimates as we accumulate samples and watch
   for them to stabilize.
#. Running multiple, independent analyses to make sure the estimates stabilize
   to similar values each time.

.. _prior_on_divergence_times:

Prior on divergence times
=========================

As mentioned in section ":ref:`bayesian_divergence_model_choice`", the prior
distribution used for divergence-time parameters can have a very large affect
on the posterior probabilities of the divergence models, due to how the priors
weight the marginal likelihoods of the models
:cite:`Oaks2014reply` :cite:`Oaks2014dpp`.
So, we have to take care when we choose a probability distribution to represent
our prior knowledge about the divergence times of our three lizard population
pairs.
This is because this prior has a strong influence on the marginal likelihoods
of the divergence models.
As we add divergence-time parameters to a divergence model, the model is forced
to integrate over a *much* greater parameter space.
For example, let's consider a uniform prior of 0 to 5 million years on
divergence times.
The :math:`\divModel{1}` model above only has a single divergence-time
parameter, and so would have to integrate over a single dimension from 0 to 5.
The :math:`\divModel{2}`, :math:`\divModel{3}`, and :math:`\divModel{4}` models
each have 2 divergence-time parameters, and so have to integrate over a
:math:`5 \times 5` square.
The :math:`\divModel{3}` model has three divergence-time parameters, and so has
to integrate over a :math:`5 \times 5 \times 5` cube.
Now imagine we were comparing 20 pairs of populations; the most general model
would integrate over a :math:`5^{20}` multidimensional space!!

If using a uniform distribution to represent our prior uncertainty, we
necessarily have to put a lot of prior density in unlikely regions of parameter
space to avoid excluding the true divergence times before we even start the
analysis :cite:`Oaks2014reply`.
For example, we might suspect all three pairs of lizard populations diverged
within the last 5 million years.
However, to feel confident that we are not excluding the (unknown) true values
of the divergence times *a priori*, we might need to specify a prior of 0 to 10
million years.
The consequence of this is that we are placing the same amount of prior density
between 5--10 million years as we are between 0--5, even though we suspect the
former is quite improbable *a priori*.
So, why does this matter?
Well, if we were correct *a priori*, and the likelihood of all three species
diverging between 5--10 million years is nearly zero, we have imposed a very
strong "penalty" for models with more divergence-time parameters.
The :math:`\divModel{3}` will integrate over a :math:`5 \times 5 \times 5` cube
with very small likelihood, but a lot of prior weight, which will result in a
very small marginal (or "average") likelihood, and thus a small posterior
probability.
Again, imagine the marginal likelihood of the most general model if we were
comparing 20 lizard species!!
The :math:`\divModel{1}` might have the largest marginal likelihood (even if it
does not explain the data very well) simply because it is "averaged" over less
space with small likelihood and large prior density.

.. _likelihood_surface:
.. figure:: /_static/marginal-plot-3d.png
   :align: center
   :width: 600 px
   :figwidth: 60 %
   :alt: likelihood surface plot
   
   The likelihood surface of a divergence model with two divergence-time
   parameters.
   The white line shows the likelihood of the co-divergence (1-parameter)
   model, and the red dashed line shows the outline of a uniform prior.
   Despite capturing much less of the likelihood density, the constrained
   1-parameter model has a larger *marginal* likelihood in this example.

If we use a uniform prior in this case, we will likely end up with strong
posterior support for a model with shared divergence times, even if the three
pairs of lizard populations diverged at quite different times.
|msbayes|_ uses a uniform prior on divergence times, and this is a key reason
it will often support models of highly clustered divergences even when taxa
diverge randomly over quite broad timescales; see :cite:`Oaks2012` and
:cite:`Oaks2014reply` for more details.

A simple solution to this problem is to use a more flexible prior on divergence
times that allows us to better represent our prior uncertainty.
In this example, we would like to specify a prior that places most of the prior
density on divergence times between 0--5 million years, but allows for a tail
with low density to capture our prior uncertainty up to 10 million years.
If we look at just one divergence-time dimension (Figure gamma_prior_), we can
see in the gamma_prior_ figure below that a gamma probability distribution
works quite well for this; |dpp-msbayes|_ uses a gamma prior on divergence
times.

.. _gamma_prior:
.. figure:: /_static/marginal-plot-2d.png
   :align: center
   :width: 600 px
   :figwidth: 60 %
   :alt: gamma prior plot
   
   The flexibility of a gamma distribution (blue) to better represent prior
   knowledge about divergence times. The black line represents the likelihood
   density, and the red line is a uniform prior.

From a lot of analyses of simulated and empirical data, I have found that by
placing much less prior weight in unlikely regions of parameters space, gamma
priors on divergence times is much less likely to spuriously support models of
shared divergences across taxa :cite:`Oaks2014dpp`.

A non-parametric approach to divergence models
==============================================

In addition to placing priors on all of the parameters of the divergence
models, we also have to place a prior on the divergence models themselves.
This can be a bit tricky, because there can be *a lot* of divergence models.
In our example of :math:`\npairs{} = 3` lizard species above, we saw there were
five possible models of divergence (i.e., there were five possible ways to
assign the three species to possible divergence-time parameters):
There was only one way to assign the species to both one and three divergence
events,
and there were three ways to assign the three species to two divergence events.
More generally, the number of ways to assign :math:`\npairs{}` taxa to
:math:`n` divergence events is the
`Stirling number of the second kind
<http://en.wikipedia.org/wiki/Stirling_numbers_of_the_second_kind>`_.
Taking this a step further, there can be anywhere from :math:`1` to
:math:`\npairs{}` divergence events, and so to calculate
the total number of possible divergence models, we need to calculate
the
`Stirling number of the second kind
<http://en.wikipedia.org/wiki/Stirling_numbers_of_the_second_kind>`_
for :math:`1,2,\ldots,\npairs{}` divergence events and sum them all up
(this is the `Bell number <http://en.wikipedia.org/wiki/Bell_number>`_
:cite:`Bell1934`).
For 3, 5, 10, and 20 taxa, there are 5, 52, 115975, and 51724158235372
possible divergence models, respectively.
The number of possible models quickly explodes as we compare more taxa!
So, how do we put a prior on all of them?!

|msbayes|_ avoids this problem by assigning equal prior probability to all
possible *number* of divergence events (divergence-time parameters).
However, it is important to realize that this strategy can create a **very**
non-uniform prior on the divergence *models*.
This is because there are *many* more ways to assign our taxa to intermediate
numbers of divergence events.
For example, if we are comparing 10 taxa, Figure number_of_models_, shows the
number of possible assignments of those taxa to :math:`1,2,\ldots,10`
divergence events (i.e., the number of possible divergence models with
:math:`1,2,\ldots,10` divergence-time parameters).

.. _number_of_models:
.. figure:: /_static/number-of-div-models-10.png
   :align: center
   :width: 600 px
   :figwidth: 60 %
   :alt: number-of-models plot
   
   The number of divergence models for 10 taxa.

As a result, if we look at the prior probability of each divergence model with
:math:`1,2,\ldots,10` divergence-time parameters, we see that it is *very*
non-uniform (Figure probability_of_models_).

.. _probability_of_models:
.. figure:: /_static/prob-of-div-models-10.png
   :align: center
   :width: 600 px
   :figwidth: 60 %
   :alt: probability-of-models plot
   
   The average prior probability of a divergence model with 1 to 10 divergence
   events.

