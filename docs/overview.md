---
layout: docs
title: Overview
permalink: docs/overview/
---

PyMsBayes is a multi-processing Python API for approximate Bayesian computation
(ABC), and provides a multi-processing interface to the comparative
phylogeographic software package, msBayes. The package is written by [Jamie
Oaks](http://www.phyletica.com).

PyMsBayes is essentially a multi-processing wrapper around some popular ABC
tools. Executables of the following tools come bundled with PyMsBayes:

 *  dpp-msbayes (<https://github.com/joaks1/dpp-msbayes.git>)
 *  msBayes (<http://msbayes.sourceforge.net/>)
 *  ABCtoolbox (<http://www.cmpg.iee.unibe.ch/content/softwares__services/computer_programs/abctoolbox/index_eng.html>)
 *  EuReject of the ABACUS package (<https://github.com/joaks1/abacus.git>)

Executables for Linux and Mac are included. The Linux executables are
statically linked and should run on any 64-bit Linux system. The Mac
executables are universal (built for architectures ppc, i386, and x86_64) and
*mostly* statically-linked, and so *should* work on any Mac. But, Macs are
weird, so no promises. If you are having problems, please let me know.


