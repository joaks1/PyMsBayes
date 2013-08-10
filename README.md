Overview
========

PyMsBayes is a multi-processing Python API for approximate Bayesian computation
(ABC), and provides a multi-processing interface to the comparative
phylogeographic software package, msBayes.

Binaries of the following tools come bundled with PyMsBayes:

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
    $ python setup.py install

If the install fails due to a permission error, try:

    $ sudo python setup.py install

To run the test suite, use:

    $ python setup.py test

If you plan to develop the code, install via:

    $ python setup.py develop

Citing PyMsBayes
================

If you publish results obtained using this software, please cite the software
and the appropriate citations for the bundled tools that were used:

 *  msBayes:

    > Huang, W., N. Takebayashi, Y. Qi, and M. J. Hickerson, 2011.
    > MTML-msBayes: Approximate Bayesian comparative phylogeographic
    > inference from multiple taxa and multiple loci with rate
    > heterogeneity. BMC Bioinformatics 12:1.

 *  ABCtoolbox:

    > Wegmann, D., C. Leuenberger, S. Neuenschwander, and L.
    > Excoffier, 2010. ABCtoolbox: a versatile toolkit for approximate Bayesian
    > computations. BMC Bioinformatics 11:116.

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

