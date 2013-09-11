#!/bin/sh

###############################################################################
# functions to control pymsbayes dev environment

function echo_current_pymsbayes_env () {
    echo "current pymsbayes environment:"
    echo "    PYMSBAYES_DEVELOPER: $PYMSBAYES_DEVELOPER"
    echo "    PYMSBAYES_LOGGING_LEVEL: $PYMSBAYES_LOGGING_LEVEL"
    echo "    PYMSBAYES_TESTING_LEVEL: $PYMSBAYES_TESTING_LEVEL"
    echo "    PYMSBAYES_DUMP_DEBUG_INFO: $PYMSBAYES_DUMP_DEBUG_INFO"
    echo "    PYMSBAYES_MEMORY_LOGGING_FREQUENCY: $PYMSBAYES_MEMORY_LOGGING_FREQUENCY"
}
function activate_pymsbayes_dev_env () {
    echo "activating pymsbayes development environment..."
    PYMSBAYES_DEVELOPER=1
    PYMSBAYES_LOGGING_LEVEL=debug
    PYMSBAYES_TESTING_LEVEL=EXHAUSTIVE
    PYMSBAYES_DUMP_DEBUG_INFO=0
    PYMSBAYES_MEMORY_LOGGING_FREQUENCY=0
    export PYMSBAYES_DEVELOPER PYMSBAYES_LOGGING_LEVEL PYMSBAYES_TESTING_LEVEL \
            PYMSBAYES_DUMP_DEBUG_INFO PYMSBAYES_MEMORY_LOGGING_FREQUENCY
}
function deactivate_pymsbayes_dev_env () {
    echo "deactivating pymsbayes development environment..."
    unset PYMSBAYES_DEVELOPER
    unset PYMSBAYES_LOGGING_LEVEL
    unset PYMSBAYES_TESTING_LEVEL
    unset PYMSBAYES_DUMP_DEBUG_INFO
    unset PYMSBAYES_MEMORY_LOGGING_FREQUENCY
    export PYMSBAYES_DEVELOPER PYMSBAYES_LOGGING_LEVEL PYMSBAYES_TESTING_LEVEL \
            PYMSBAYES_DUMP_DEBUG_INFO PYMSBAYES_MEMORY_LOGGING_FREQUENCY
}
function pymsbayes_dev_usage () {
    echo "usage: siftdev [on] [off]\n"
    echo "  -h|--help  show help message and exit.\n"
    echo "  on         activate sate development environment.\n"
    echo "  off        deactivate sate development environment.\n"
    echo "             if no arguments provided, report state of sate environment.\n"
}
function pymsbayes_dev () {
    if [ $# -gt 1 ]
        then
        echo "function takes at most one argument."
        pymsbayes_dev_usage
    elif [ $# -eq 0 ]
        then
        echo_current_pymsbayes_env
    else
        case $1 in
        on|On|ON)
            activate_pymsbayes_dev_env
            echo_current_pymsbayes_env
            shift
            ;;
        off|Off|OFF)
            deactivate_pymsbayes_dev_env
            echo_current_pymsbayes_env
            shift
            ;;            
        -h|--help)
            pymsbayes_dev_usage
            shift
            ;;
        *)
            echo "unrecognized arguments: $@"
            pymsbayes_dev_usage
            break
            ;;
        esac
    fi
}
