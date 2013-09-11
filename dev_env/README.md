Info
====

This directory contains the shell script `pymsbayesdev.sh` for setting up
useful environmental variables for developing `PyMsBayes`. If you are going to
be working on `PyMsBayes` code, I recommend putting this script in your
favorite location (creatively `/path/to/` in the example below) and adding the
following to your `.bashrc` file:

    function pymsbayesdev () {
        source /path/to/pymsbayesdev.sh
        pymsbayes_dev $@
    }

This will allow you to toggle the development environment on/off with simply

    $ pymsbayesdev on
    $ pymsbayesdev off

Also, if you just type

    $ pymsbayesdev

it will show you the current state of `PyMsBayes`-related environmental
variables, without making any changes to them.

If you do not want to add a function to your `.bashrc` you can also simply add:

    source /path/to/pymsbayesdev.sh

and use the sourced `pymsbayes_dev` function just as `pymsbayesdev` was used in
the examples above. I like wrapping this up in a function so that the 99% of
the time I open a shell with no intention of working on `PyMsBayes` code, I am
not unnecessarily sourcing all the functions in `pymsbayesdev.sh`.

Alternatively, you can leave your `.bashrc` file as is, and simply type in

    $ source /path/to/pymsbayesdev.sh

at the prompt whenever you want to activate the `PyMsBayes` dev environment.
I'm just too lazy for that.

