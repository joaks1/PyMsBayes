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

.. divergence_model_111:
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
:math:`M_1 = 111` to show that this divergence model assigns population pairs
1, 2, and 3 to divergence-time parameter 1.
However, this is only one possible divergence model, and is the most
constrained.
With three population pairs, there are 4 other possible models of divergence (5
total possible models).
Three of these models have two divergence-time parameters.
We can assign population pair 1 to a second divergence-time parameter to get
divergence model :math:`M_2 = 211`, as shown in the divergence_model_211_
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
divergence model :math:`M_3 = 121`, as shown in the divergence_model_121_
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
model :math:`M_4 = 112`, as shown in the divergence_model_112_ figure.

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
:math:`M_5 = 123`), as shown in the divergence_model_123_ figure. This is the
most general model of divergence.

.. _divergence_model_123:
.. figure:: /_static/div-model-cartoon-123.png
   :align: center
   :width: 600 px
   :figwidth: 60 %
   :alt: divergence model 123
   
   A cartoon showing the most general model of divergence where all three
   population-pairs of lizards that diverge at unique times.

Bayesian divergence-model choice
================================

If we go out and sample individuals from each of the lizard populations, and
from those individuals collect DNA sequence data from one or more orthologous
loci per pair of populations, we can use these data to infer the temporal
distribution of the population divergences across the three lizard species.
We know that the sequences of a locus of a pair are related by a genealogy, and
that the shape of this genealogy is governed by demographic processes.
We also know that the variation we see in the contemporary sequences
accumulated as the sequences evolved via mutational processes along the
genealogy.
So, if we use the coalescent to model the demographic processes controlling
the branching rate of the gene tree, and a Markov chain to model the
mutational processes, and make an assumption about the relative mutation
rates among the three lizard species, we can estimate the relative times
that the three pairs of populations diverged.
So, we can modify our cartoon of model :math:`M_5 = 123` to better
represent data and model.
