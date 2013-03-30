#! /usr/bin/env python

import os
import sys
import glob

BIN_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT_DIR = os.path.abspath(os.path.dirname(BIN_DIR))
CONFIG_DIR = os.path.abspath(os.path.join(PROJECT_DIR, 'configs'))
RESULT_DIR = os.path.abspath(os.path.join(PROJECT_DIR, 'results'))

def main():
    sys.stdout.write("%s" % PROJECT_DIR)

if __name__ == '__main__':
    main()

