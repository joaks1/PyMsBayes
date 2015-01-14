.. role:: bolditalic
.. role:: hlight 
.. role:: codehlight 

.. _selecting_priors:

****************
Selecting Priors
****************

.. _gamma_intro:

An introduction to the gamma probability distribution
=====================================================

|dpp-msbayes| uses the gamma probability distribution for priors on many of the
parameters in the model.
Thus, before we dive into all the priors of the |dpp-msbayes| model, it will be
helpful to

#. :ref:`introduce some basic characteristics of gamma
   distributions<gamma_basics>`
#. :ref:`introduce some graphical tools to help us visualize and select gamma
   priors<gamma_in_r>`


.. _gamma_basics:

Gamma basics
------------

All of the gamma priors in |dpp-msbayes| have two parameters: The shape
(:math:`\gshape`) and scale (:math:`\gscale`) parameters.
The mean and variance of a gamma distribution are :math:`\gshape\gscale` and
:math:`\gshape\gscale^2`, respectively.

As you might have guessed, shape parameter controls the shape of the
distribution, while the scale parameter controls the scale.
You can think of it this way: all gamma distributions with the same value of
the shape parameter have the same shape, and differences among them in the
scale parameter just "re-scales" the x-axis.
When the shape is :math:`\gshape = 1`, the gamma becomes an exponential
distribution with a mean of :math:`\gscale`.
When the shape is less than or equal to 1, the mode, or "peak", of the
distribution is at zero.
When the shape is greater than 1, the mode is greater than zero
(:math:`(\gshape-1)\gscale`).

When thinking about the gamma distribution as a prior, the shape parameter
represents how much *a prior* knowledge we have about a parameter.
A larger value means we are more certain about the value of the parameter *a
priori*, whereas a smaller means we are less certain.
The scale parameter simply "shifts" the distribution along the x-axis to "fit"
our prior expectations.


.. _gamma_in_r:

Using R to help visualize gamma priors
--------------------------------------

Let's say we have strong prior knowledge that the value of a parameter is
around 10.0, and we are pretty sure that it's less than 15.0.
Because, we want the prior to be centered near 10.0, we can take advantage of
the fact that the mean of a gamma distribution should be around
:math:`\gshape\gscale = 10.0`.
We can use some quick-and-dirty R to get a rough idea of what a gamma prior
with this mean and a shape of 1.0 (i.e., an exponential prior) will look like:

.. code-block:: r

    > x = rgamma(100000, shape=1.0, scale=10.0)
    > hist(x)

This will look something like:

.. _gamma_1_10_hist:
.. figure:: /_static/gamma_1_10_hist.png
    :align: center
    :width: 600 px
    :figwidth: 60 %
    :alt: gamma(1, 10) histogram

    Gamma(1, 10) histogram
   
With a little more typing in the R prompt, we can get a more accurate
illustration of this prior:

    
.. code-block:: r

    > mx = qgamma(0.999, shape=1.0, scale=10.0)
    > x = seq(from=0, to=mx, by=mx/1000)
    > dens = dgamma(x, shape=1.0, scale=10.0)
    > plot(x, dens, type='l')

which will give us something like:

.. _gamma_1_10_plot:
.. figure:: /_static/gamma_1_10_plot.png
    :align: center
    :width: 600 px
    :figwidth: 60 %
    :alt: gamma(1, 10) plot

    Gamma(1, 10)

We said that we were quite confident the parameter is less than 15.0.
It looks like this prior puts too much prior probability on values greater than
15.
We can use R to see that the prior probability of this parameter being greater
than 15.0 is about 0.22:

.. code-block:: r

    > pgamma(15.0, shape=1.0, scale=10.0, lower.tail=F)
    [1] 0.2231302

So, it looks like we need to increase our prior knowledge. Let's try a shape
parameter of 2.0 (we need to adjust the scale parameter to 5.0 to keep the mean
of 10.0):


.. code-block:: r

    > pgamma(15.0, shape=2.0, scale=5.0, lower.tail=F)
    [1] 0.1991483

Hmmm... Still too much prior probability on values over 15. Let's try
a shape of 10.0:

.. code-block:: r

    > pgamma(15.0, shape=10.0, scale=1.0, lower.tail=F)
    [1] 0.06985366
    
Let's assume this fits our prior expectation pretty well (i.e., we want to
state *a priori* that the probability of the prior being less than 15 is 0.93).
Let's take a look at this gamma prior with a shape and mean of 10.0:

.. code-block:: r

    > mx = qgamma(0.999, shape=10.0, scale=1.0)
    > x = seq(from=0, to=mx, by=mx/1000)
    > dens = dgamma(x, shape=10.0, scale=1.0)
    > plot(x, dens, type='l')

.. _gamma_10_1_plot:
.. figure:: /_static/gamma_10_1_plot.png
    :align: center
    :width: 600 px
    :figwidth: 60 %
    :alt: gamma(10, 1) plot

    Gamma(10, 1)


.. contents:: Priors 
    :local:


.. _concentration_parameter:

Concentration parameter of the Dirichlet process
================================================

We have to choose a gamma-distributed prior for the concentration parameter
(:math:`\alpha`) of the Dirichlet process controlling the assignment of taxa to
divergence events.
From the ":ref:`dpp`" section, we know that as the concentration parameter
decreases, we are putting more prior probability on models of divergence that
are more clustered (i.e., models with fewer shared divergence events).
Alternatively, as we increase :math:`\alpha`, we place more prior probability
on divergence models with less co-divergence among taxa.

|pmb|_ comes with a program named |ldppsum| that helps us to
get a feel for this.
Let's say we have sequence data from 10 taxa, and we want to know what
value of the concentration parameter corresponds with a prior mean
of 5 divergence events.
We can use |ldppsum| to calclulate this by typing:

.. parsed-literal::

    $ |dppsum| ncats 5 10

The output should look like::
    
    number of elements = 10
    concentration parameter = 3.30149636133
    expected number of categories = 5.0

This tells thus that divergence models generated under a Dirichlet process with
10 taxa and a concentration parameter of about 3.3 will have 5 divergence
events (parameters) on average.
We can confirm this by typing:
    
.. parsed-literal::

    $ |dppsum| concentration 3.3 10

Which reports::

    number of elements = 10
    concentration parameter = 3.3
    expected number of categories = 4.99909319002

Ok, that's useful, but what about the probability of other numbers of events?
Well, we can use the ``--reps`` option to tell |dppsum| to use simulations to
estimate such probabilities:

.. parsed-literal::

    $ |dppsum| ncats 5 10 --reps 10000

This generates 10000 random divergence models under a Dirichlet process prior,
and reports the estimated prior probabilites of all possible values for the
number of divergence event. It also reports the exact value of the number of
possible divergence models for each number of divergence events::

    number of elements = 10
    concentration parameter = 3.30149636133
    expected number of categories = 5.0
    
    Starting simulations to estimate probabilities...
    Using seed 436471208
    
    Estimated probabilities of the number of categories:
    	p(ncats = 1) = 0.0024 (n = 1)
    	p(ncats = 2) = 0.0280 (n = 511)
    	p(ncats = 3) = 0.1012 (n = 9330)
    	p(ncats = 4) = 0.2240 (n = 34105)
    	p(ncats = 5) = 0.2912 (n = 42525)
    	p(ncats = 6) = 0.2048 (n = 22827)
    	p(ncats = 7) = 0.1080 (n = 5880)
    	p(ncats = 8) = 0.0348 (n = 750)
    	p(ncats = 9) = 0.0048 (n = 45)
    	p(ncats = 10) = 0.0008 (n = 1)

For example, this output tells us that the prior probability of a divergence
model with 2 divergence-time parameters, under aDirichlet process with 10 taxa
and a concentation parameter of about 3.3, is approximately 0.028.
It also tells us that there are 511 possible divergence models with 2
divergence events (i.e., 511 different ways of assigning our taxa to 2
divergence events).

This is all well and good, but in our |dpp-msbayes|_ :ref:`configuration
file<config>`, we need to specify the shape and scale parameters for a gamma
prior on the concentration parameter.
No problem, |dppsum| can help us with that too.
If we want to essentially fix the concentration parameter to 3.3, we can
specify a very large shape parameter for the gamma prior:

.. parsed-literal::

    $ |dppsum| ncats 5 10 --reps 10000 --shape 1000

The output will be something like:

.. parsed-literal::

    number of elements = 10
    concentration parameter = 3.30149636133
    expected number of categories = 5.0
    shape = :codehlight:`1000.0`
    scale = :codehlight:`0.00330149636133`
    
    Starting simulations to estimate probabilities...
    Using seed 428982720
    
    Estimated probabilities of the number of categories:
    	p(ncats = 1) = 0.0012 (n = 1)
    	p(ncats = 2) = 0.0292 (n = 511)
    	p(ncats = 3) = 0.1104 (n = 9330)
    	p(ncats = 4) = 0.2296 (n = 34105)
    	p(ncats = 5) = 0.2756 (n = 42525)
    	p(ncats = 6) = 0.2076 (n = 22827)
    	p(ncats = 7) = 0.1048 (n = 5880)
    	p(ncats = 8) = 0.0328 (n = 750)
    	p(ncats = 9) = 0.0080 (n = 45)
    	p(ncats = 10) = 0.0008 (n = 1)

As you can see, aside from some estimation error due to a finite number of
simulation replicates, the probabilities are nearly identical to our previous
command.
Notice that |dppsum| now reports the ``shape`` and ``scale`` parameters
(highlighted above); these correspond the shape and scale parameters of a gamma
prior on the concentration parameter.
So, if we put the following in our configuration file:

.. parsed-literal::

    concentrationShape = :codehlight:`1000.0`
    concentrationScale = :codehlight:`0.00330149636133`

we will be using a Dirichlet process prior with a (nearly) fixed concentration
parameter of 3.3, which, on average, yields divergence models with 5 divergence
events.

