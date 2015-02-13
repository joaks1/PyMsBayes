.. role:: bolditalic
.. role:: hlight 
.. role:: codehlight 

.. contents:: 
    :local:
    :depth: 3

.. _simulation_analysis:

**************************************
A Simulation-based Validation Analysis
**************************************

.. parsed-literal::

    $ |dmc| -o simulated-data-configs/exponential-02.cfg simulated-data-configs/exponential-04.cfg simulated-data-configs/exponential-06.cfg -p dpp-simple.cfg msbayes.cfg -r 100 -n 5000

    total_duration = 0:08:56.954086

Add ``--no-global-estimate``? 900 rejection/regressions reduced to 600
