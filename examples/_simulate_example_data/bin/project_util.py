#! /usr/bin/env python

import os
import sys

BIN_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT_DIR = os.path.dirname(BIN_DIR)
EXAMPLE_DIR = os.path.dirname(PROJECT_DIR)
LIZARD_DIR = os.path.join(EXAMPLE_DIR, 'lizards')
LIZARD_SEQ_DIR = os.path.join(LIZARD_DIR, 'sequences')
LIZARD_CONFIG_DIR = os.path.join(LIZARD_DIR, 'configs')
DATA_DIR = os.path.join(PROJECT_DIR, 'data')
SP_TREE_DIR = os.path.join(DATA_DIR, 'species-tree')
SP_TREE_PATH = os.path.join(SP_TREE_DIR, 'species-tree.nex')
SP_TREE_MT_PATH = os.path.join(SP_TREE_DIR, 'species-tree-mt.nex')
GENE_TREE_DIR = os.path.join(DATA_DIR, 'gene-trees')
SEQ_DIR = os.path.join(DATA_DIR, 'sequences')

def main():
    sys.stdout.write("{0}".format(PROJECT_DIR))

if __name__ == '__main__':
    main()

