#! /usr/bin/env python

import os
import sys
import multiprocessing

from pymsbayes.utils import messaging
from pymsbayes import workers
from pymsbayes import manager
import project_util

_LOG = messaging.get_logger(__name__)

def simulate_alignments(
        num_individuals_per_pop = 5,
        hky_kappa = 10,
        alignment_lengths = [600, 400, 500, 550, 350, 450],
        ):
    if not os.path.exists(project_util.GENE_TREE_DIR):
        os.mkdir(project_util.GENE_TREE_DIR)
    if not os.path.exists(project_util.SEQ_DIR):
        os.mkdir(project_util.SEQ_DIR)
    tree_workers = []
    seq_workers = []
    for i, align_len in enumerate(alignment_lengths):
        tree_path = os.path.join(project_util.GENE_TREE_DIR,
                'gene-tree-locus-{0}.nex'.format(i))
        seq_path = os.path.join(project_util.SEQ_DIR,
                'locus-{0}.nex'.format(i))
        num_locus_copies = 2
        sp_tree_path = project_util.SP_TREE_PATH
        if i == 0:
            tree_path = os.path.join(project_util.GENE_TREE_DIR,
                    'gene-tree-locus-mt.nex')
            seq_path = os.path.join(project_util.SEQ_DIR,
                    'locus-mt.nex')
            num_locus_copies = 1
            # Using a 4-fold faster mutation rate for a mitochondrial locus
            # (all branches of species tree are multiplied by 4) so that the
            # difference in ploidy cancels out. E.g., 2Nu = 0.5N4u
            sp_tree_path = project_util.SP_TREE_MT_PATH
        if os.path.exists(tree_path):
            raise Exception('script has already been run: output files already '
                    'exist... aborting')
        tw = workers.GeneTreeSimWorker(args = [
                '-n', '1',
                '-t', num_individuals_per_pop * num_locus_copies,
                '-o', tree_path,
                sp_tree_path])
        sw = workers.SequenceSimWorker(args = [
                '-m', 'HKY,1,{0}'.format(hky_kappa),
                '-n', align_len,
                tree_path],
                output_path = seq_path)
        tree_workers.append(tw)
        seq_workers.append(sw)
    tree_workers = manager.Manager.run_workers(tree_workers,
            num_processors = multiprocessing.cpu_count())
    seq_workers = manager.Manager.run_workers(seq_workers,
            num_processors = multiprocessing.cpu_count())
    return 

def main():
    simulate_alignments()


if __name__ == '__main__':
    main()

