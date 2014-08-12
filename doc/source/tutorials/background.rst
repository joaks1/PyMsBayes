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
For example, if an event (a black rectange) split a community of species
260,000 years ago, we might expect to find that divergences across multiple
pairs of populations co-distributed across the barrier were temporally
clustered, as shown in the divergence_model_111_ figure.

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
:math:`\divModel{1} = 111` to show that this divergence model assigns population pairs
1, 2, and 3 to divergence-time parameter 1.
However, this is only one possible divergence model, and is the most
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
:math:`\divModel{5} = 123`), as shown in the divergence_model_123_ figure. This is the
most general model of divergence.

.. _divergence_model_123:
.. figure:: /_static/div-model-cartoon-123.png
   :align: center
   :width: 600 px
   :figwidth: 60 %
   :alt: divergence model 123
   
   A cartoon showing the most general model of divergence where all three
   population-pairs of lizards that diverge at unique times.

If we go out and sample individuals from each of the lizard populations, and
from those individuals collect DNA sequence data from one or more orthologous
loci per pair of populations, we can use these data to infer the temporal
distribution of the population divergences across the three lizard species.
We know that the sequences of a locus of a pair are related by a genealogy, and
that the shape of this genealogy is governed by demographic processes.
We also know that the variation we see in the contemporary sequences
accumulated as the sequences evolved via mutational processes along the
genealogy.
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
Let's also use :math:`\alignmentVector` to represent the sequence alignments
of all the loci we sequenced for the individuals we sampled from all three
pairs of lizard populations.
If we make assumptions about the relative rates of mutations and the relative
generation times among the three lizard species, we can calculate the posterior
probability distribution of the divergence times (and other nuisance
parameters) given the data and one of the models of divergence using
Bayes rule:

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
:math:`\divModel{1}`, where this average is weighted by the prior over the
entire space of the model.
If we calculate this marginal likelihood of all five possible divergence
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

However, we cannot calculate all of those integrals exactly, so we will need to
use a numerical integration tool to approximate the posterior.
For example, we could derive the likelihood term in Equation :eq:`postdensity`
and use Markov chain Monte Carlo to approximate the posterior.
However, we can avoid deriving and calculating the likelihood if we use
approximate likelihoods.
Thus we will use a numerical integration algorithm that uses approximate
likelihoods to approximate the posterior.
