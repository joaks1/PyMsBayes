********
Overview
********

|pmb|_ is a multi-processing Python API for approximate Bayesian computation
(ABC), and provides a multi-processing interface to the comparative
phylogeographic software package, |msbayes|_. The package is written by |jro|_.

PyMsBayes is essentially a multi-processing wrapper around some popular ABC
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
