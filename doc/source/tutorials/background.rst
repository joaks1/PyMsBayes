.. _background:

**********
Background
**********

Comparative divergence models
=============================

Biogeographers often seek to explain diversity on historical events.
That is, using genetic data from contemporary populations, we would like to
infer diversification patterns and see if they support patterns predicted by
past event(s) of interest (e.g., islands fragmented by rises in sea level,
changes in climate fragmenting communities into refugia, etc.).
For example, if an event (a black rectangle) split a community of species
260,000 years ago, we might expect the divergences across multiple species
co-distributed across the barrier to be temporally clustered.
More specifically, let's say we are interested in investigating three pairs of
lizard populations, where the populations of each pair occur on opposite sides
of the putative barrier.
If the historical event caused divergence, we would expect that the pairs of
lizard populations (or some subset of them) diverged at the same time, as shown
in the divergence_model_111_ figure.

.. _divergence_model_111:
.. figure:: /_static/div-model-cartoon-111.png
   :align: center
   :width: 600 px
   :figwidth: 60 %
   :alt: divergence model 111
   
   A cartoon showing three population-pairs of lizards that co-diverge due to an
   event 260,000 years ago.

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
We can assign population pair 1 to a second divergence-time parameter to get
divergence model :math:`\divModel{2} = 211`, as shown in the divergence_model_211_
figure.

.. _divergence_model_211:
.. figure:: /_static/div-model-cartoon-211.png
   :align: center
   :width: 600 px
   :figwidth: 60 %
   :alt: divergence model 211
   
   A cartoon showing population-pair 1 assigned to divergence-time parameter 2,
   and population pairs 2 and 3 assigned to divergence-time parameter 1.

We can also assign population pair 2 to divergence-time parameter 2 to get
divergence model :math:`\divModel{3} = 121`, as shown in the divergence_model_121_
figure.

.. _divergence_model_121:
.. figure:: /_static/div-model-cartoon-121.png
   :align: center
   :width: 600 px
   :figwidth: 60 %
   :alt: divergence model 121
   
   A cartoon showing population-pair 2 assigned to divergence-time parameter 2,
   and population pairs 1 and 3 assigned to divergence-time parameter 1.

And for the last possible divergence model with two divergence-time parameters,
we assign population pair 3 to divergence-time parameter 2 to get divergence
model :math:`\divModel{4} = 112`, as shown in the divergence_model_112_ figure.

.. _divergence_model_112:
.. figure:: /_static/div-model-cartoon-112.png
   :align: center
   :width: 600 px
   :figwidth: 60 %
   :alt: divergence model 112
   
   A cartoon showing population-pair 3 assigned to divergence-time parameter 2,
   and population pairs 1 and 2 assigned to divergence-time parameter 1.

Finally, we can add a third divergence-time parameter so that each pair of
populations is assigned to its own divergence-time parameter (divergence model
:math:`\divModel{5} = 123`), as shown in the divergence_model_123_ figure.
This is the most general model of divergence.

.. _divergence_model_123:
.. figure:: /_static/div-model-cartoon-123.png
   :align: center
   :width: 600 px
   :figwidth: 60 %
   :alt: divergence model 123
   
   A cartoon showing the most general model of divergence where all three
   population-pairs of lizards that diverge at unique times.

Being energetic herpetologists, we go out and sample individuals from each of
the lizard populations, and from those individuals collect DNA sequence data
from one or more orthologous loci per pair of populations.
You can find our sequence data in fasta format in the |lizard-seq-dir|_
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
   population-pairs of lizards that diverge at unique times.

We can use these data to infer the temporal distribution of the population
divergences across the three lizard species.

Bayesian divergence-model choice
================================

In the figures above, we used :math:`\divTimeMap{1}, \divTimeMap{2},` and
:math:`\divTimeMap{3}` to represent the divergence times of the three
pairs of lizard populations. Now, let's use :math:`\divTimeMapVector`
to represent all three divergence times; that is, 
:math:`\divTimeMapVector = \divTimeMap{1}, \divTimeMap{2}, \divTimeMap{3}`.
The number of unique divergence-time values (parameters) within
:math:`\divTimeMapVector`, and the assignment of the lizard species to these
values, depends on the divergence model.
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
If we assume a Markov chain model of nucleotide substitution, we can
calculate the probability of the sequence data given the genealogies and a set
of parameter values for the substitution model.
If we further assume a coalescent model of ancestral processes, we can
calculate the probability of the genealogies given the parameter values
for the sizes of the populations.
Let's lump the mutational and demographic parameters of the substitution and
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
    p(\alignmentVector \given \divModel{2}) p(\divModel{2}) +
    p(\alignmentVector \given \divModel{5}) p(\divModel{5}) }

Or, more generally, we can calculate the posterior probability of any
divergence model ":math:`i`" using:

.. math::
    :label: postmass

    p(\divModel{i} \given \alignmentVector) = \frac{ p(\alignmentVector \given
    \divModel{i}) p(\divModel{i}) }{ \sum_{i} p(\alignmentVector \given
    \divModel{i}) p(\divModel{i}) }

This is essentially the relative marginal likelihood of the model (it is
exactly that if assume equal prior mass for each divergence model). Thus, the
marginal likelihoods are the "guts" of Bayesian model choice.
Because we are sampling over all the divergence models, Equation :eq:`postmass`
will also give us model-averaged estimates of the divergence times for each of
our pairs of populations (i.e., we get estimates of divergence times that
account for uncertainty in divergence models).

As you might expect, because the marginal likelihoods are weighted by the
priors on parameters, it can be quite sensitive to the priors.
Thus, we have to choose the priors on parameters carefully, and should always
assess the sensitivity of our results to differences in these prior
assumptions.
So, now seems like as good of time as any to discuss priors.

Prior on divergence times
-------------------------

We have to choose a probability distribution to represent our prior knowledge
about the divergence times of our three lizard population pairs.
It turns out that is choice is very important for the model described above.
This is because this prior has a strong influence on the marginal likelihoods
of the divergence models.
As we add divergence-time parameters to a divergence model, it has to integrate
over a *much* greater parameter space.
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
analysis.
For example, we might suspect all three pairs of lizard populations diverged
within the last 5 million years.
However, to feel confident that we are not excluding the true values of the
divergence times *a priori*, we might need to specify a prior of 0 to 10
million years.
The consequence of this is that we are placing the same amount prior density
between 5--10 million years as we are between 0--5, even though we suspect
the former is quite improbable *a priori*
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
The :math:`\divModel{1}` will have the largest marginal likelihood (even if it
does not explain the data very well) because it is "averaged" over less space
with small likelihood and large prior density.

.. _likelihood_surface:
.. figure:: /_static/marginal-plot-3d.png
   :align: center
   :width: 600 px
   :figwidth: 60 %
   :alt: likelihood surface plot
   
   The likelihood surface of a divergence model with two divergence-time
   parameters.
   The white line shows the likelihood of the 1-parameter constrained model,
   and the red dashed line shows the outline of a uniform prior.
   Despite capturing much less of the prior density, the constrained
   1-parameter model has a larger *marginal* likelihood in this example.

If we use a uniform prior in this case, we will likely end up with strong
posterior support for a model with shared divergence times, even if the three
pairs of lizard populations diverged at quite different times.
|msbayes|_ uses a uniform prior on divergence times, and this is a key reason
it will often support models of highly clustered divergences even when taxa
diverge randomly over quite broad timescales; see :cite:`Oaks2012` and
:cite:`Oaks2014reply` for more details.

A simple solution to this problem is to use a more flexible prior on divergence
times that allows us to better represent our prior knowledge.
In this example, we would like to specify a prior that places most of the prior
density on divergence times between 0--5 million years, but allows for a tail
with low density to capture our prior uncertainty up to 10 million years.

.. _gamma_prior:
.. figure:: /_static/marginal-plot-2d.png
   :align: center
   :width: 600 px
   :figwidth: 60 %
   :alt: gamma prior plot
   
   The flexibility of a gamma distribution (blue) to better represent prior
   knowledge about divergence times. The black line represents the likelihood
   density, and the red line is a uniform prior.

However, we cannot calculate all of those integrals exactly, so we will need to
use a numerical integration algorithm to approximate the posterior.
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
summary statistics plotted in three dimensions.
Next, a grey sphere illustrates the "good enough" zone.
Then, we see black points accumulate, which represent the values of the three
summary statistics calculated from datasets simulated under sets of parameter
values drawn randomly from the prior.
Lastly, we see the retained sample of points that fell within our "good enough"
zone; this is our sample from the approximate posterior.

So, how do we decide how large the "good enough" zone is? Well, the smaller the
better, but this is governed arbitrarily by computational limitations.
What we do is simulate as many datasets as we are willing to wait for and then
select the desired number of them that produced summary statistics closest to
the observed summary statistics.
For example, we might draw 10 million sets of parameter values from the prior,
and keep the 10,000 sets that produced summary statistics nearest to the
observed statistics as our approximate posterior sample.
Thus, the radius of the "good enough" space is determined by the distance
between the observed summary statistics and the 10,001st nearest simulated
summary statistics.
Again, this is arbitrary; drawing 100 million samples and keeping the closest
10,000 would be better.
However, we can get a sense of whether we have evaluated a sufficient number of
samples from the prior by keeping track of the parameter estimates as we
accumulate samples and watch for them to stabilize.

