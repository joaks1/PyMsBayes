********
Overview
********

|pmb|_ provides a multi-processing interface to the comparative phylogeographic
software packages |dpp-msbayes|_ and |msbayes|_.
More generally, |pmb|_ can serve as a multi-processing Python API for inferring
comparative models of diversification via approximate Bayesian computation
(ABC).
For more details about the methods implemented in |pmb|_ checkout :ref:`the
background section<background>`.
You should also be aware of :ref:`relevant caveats<caveats>` before using the
method.
The package is written by |jro|_.

|pmb|_ is essentially a multi-processing wrapper around some popular ABC
tools. Executables of the following tools come bundled with PyMsBayes:

 *  |dpp-msbayes|_
 *  |msbayes|_
 *  |abctb|_
 *  |eureject|_ of the |abacus|_ package

Executables for Linux and Mac are included. The Linux executables are
statically linked and should run on any 32-bit or 64-bit Linux system. The Mac
executables are universal (built for architectures ppc, i386, and x86_64) and
*mostly* statically-linked, and so *should* work on any Mac. But, Macs are
weird, so no promises. See the ":ref:`installation`" section for more details
about installing the bundled tools.
