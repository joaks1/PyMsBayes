---
layout: docs
title: Requirements
permalink: docs/requirements/
---

PyMsBayes has only been tested on version 2.7 of Python. The only Python
dependency is ConfigObj (<http://www.voidspace.org.uk/python/configobj.html>).
This dependency should get installed automatically when you follow the
installation instructions for PyMsBayes below. However, if it does not, it can
easily be installed using `pip install configobj` or `easy_install configobj`;
see <http://www.voidspace.org.uk/python/configobj.html#installing> for more
details.

Optional additions
------------------

If you wish to use any of the plotting features of `PyMsBayes`, you will need
to install `matplotlib` (<http://matplotlib.org/>).

For post-hoc regression adjustment of posterior samples, the general linear
model (GLM) method of Leuenberger and Wegmann (2010) is the default, and does
not require any extra dependencies.  However, if you wish to use the weighted
local-linear regression (multinomial logistic regression for discrete
parameters) method of Beaumont et al. (2002) via the API, then you will need to
install the following R packages: VGAM, KernSmooth, locfit, and optparse.

