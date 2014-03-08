Overview
========

PyMsBayes is a multi-processing Python API for approximate Bayesian computation
(ABC), and provides a multi-processing interface to the comparative
phylogeographic software package, msBayes. The package is written by [Jamie
Oaks](http://www.phyletica.com).

Binaries of the following tools come bundled with PyMsBayes:

 *  dpp-msbayes (<https://github.com/joaks1/dpp-msbayes.git>)
 *  msBayes (<http://msbayes.sourceforge.net/>)
 *  ABCtoolbox (<http://www.cmpg.iee.unibe.ch/content/softwares__services/computer_programs/abctoolbox/index_eng.html>)
 *  EuReject of the ABACUS package (<https://github.com/joaks1/abacus.git>)

Requirements
============

PyMsBayes has only been tested on version 2.7 of Python.

Installation
============

You should have Python 2.7 installed.

Clone Git repository and install using:

    $ git clone git://github.com/joaks1/PyMsBayes.git
    $ cd PyMsBayes
    $ python setup.py develop

If the install fails due to a permission error, try:

    $ sudo python setup.py develop

To run the test suite, use:

    $ python setup.py test

Citing PyMsBayes
================

If you publish results obtained using this software, please cite the software
and the appropriate citations for the bundled tools that were used:

 *  dpp-msbayes and euReject:

    > Oaks, J. R. (2014). An Improved Approximate-Bayesian Model-choice Method
    > for Estimating Shared Evolutionary History. arXiv:1402.6303 [q-bio:PE].
    > <http://arxiv.org/abs/1402.6303>

 *  msBayes:

    > Huang, W., N. Takebayashi, Y. Qi, and M. J. Hickerson, 2011.
    > MTML-msBayes: Approximate Bayesian comparative phylogeographic
    > inference from multiple taxa and multiple loci with rate
    > heterogeneity. BMC Bioinformatics 12:1.

 *  ABCtoolbox:

    > Wegmann, D., C. Leuenberger, S. Neuenschwander, and L.
    > Excoffier, 2010. ABCtoolbox: a versatile toolkit for approximate Bayesian
    > computations. BMC Bioinformatics 11:116.

Acknowledgements
================

This software greatly benefited from funding provided to [Jamie
Oaks](http://www.phyletica.com) from the National Science Foundation (DEB
1011423 and DBI 1308885), University of Kansas (KU) Office of Graduate Studies,
Society of Systematic Biologists, Sigma Xi Scientific Research Society, KU
Ecology and Evolutionary Biology Department, and the KU Biodiversity Institute.

Some of the functions for calculating properties of Dirichlet processes were
adapted from functions within `util.h` from the software package
[`DPPDiv`](http://phylo.bio.ku.edu/content/tracy-heath-dppdiv) version 1.0b
(Copyright Tracy Heath, Mark Holder, and John Huelsenback; licensed under GPL
v3; <http://phylo.bio.ku.edu/content/tracy-heath-dppdiv>).

License
=======

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.

See "LICENSE.txt" for full terms and conditions of usage.

