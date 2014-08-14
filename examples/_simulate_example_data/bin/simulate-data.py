#! /usr/bin/env python

import os
import sys
import re
import subprocess
import multiprocessing

from seqsift.utils import dataio
from seqsift.seqops import seqstats

from pymsbayes.utils import messaging, ToolPathManager, tempfs, functions, stats
from pymsbayes import workers
from pymsbayes import manager
import project_util

_LOG = messaging.get_logger(__name__)

number_pattern_string = r'[\d\.Ee\-\+]+'
paup_hky_string = (
        '^\s*Estimated\s+base\s+frequencies\s+=\s+'
        'A:(?P<a>{0})\s+C:(?P<c>{0})\s+'
        'G:(?P<g>{0})\s+T:(?P<t>{0})'
        '\n\s*Estimated\s+ti/tv\s+ratio\s+=\s+{0}\s*'
        '\(kappa\s+=\s+(?P<kappa>{0})\)\s*$'.format(number_pattern_string)
        )
paup_hky_pattern = re.compile(paup_hky_string, re.MULTILINE)

def write_paup_hky_file(path, nexus_data_path):
    with open(path, 'w') as out:
        out.write("#NEXUS\n\n"
                "BEGIN PAUP;\n\t"
                "set warnreset=no; set increase=auto; set warnroot=no;\n\t"
                "execute {0};\n\t"
                "DSet distance=LOGDET objective=ME base=equal rates=equal pinv=0\n\t"
                "subst=all negbrlen=setzero;\n\t"
                "NJ showtree=no breakties=random;\n"
                "END;\n\n"
                "BEGIN PAUP;\n\t"
                "Default lscores longfmt=yes;\n\t"
                "Set criterion=like;\n\t"
                "[!** Calculating HKY **]\n\t"
                "lscores 1/ nst=2 base=est tratio=est rates=equal pinv=0;\n"
                "END;\n".format(nexus_data_path)
                )

def parse_paup_log_file(path):
    with open(path, 'r') as in_stream:
        log = in_stream.read()
        match_iter = paup_hky_pattern.finditer(log)
        i = -1
        for i, m in enumerate(match_iter):
            pass
        if i < 0:
            return {'kappa': 1.0,
                    'a': 0.25,
                    'c': 0.25,
                    'g': 0.25,
                    't': 0.25}
        elif i > 0:
            raise Exception('found multiple sets of HKY parameters in '
                    '{0!r}'.format(path))
    return m.groupdict()

class PaupWorker(object):
    count = 0
    def __init__(self,
            nex_path,
            exe_path = None,
            stdout_path = None,
            stderr_path = None,
            subprocess_kwargs = {},
            tag = None):
        self.stdout_path = stdout_path
        self.stderr_path = stderr_path
        self.tag = tag
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        if not exe_path:
            exe_path = 'paup'
        self.exe_path = ToolPathManager.get_external_tool(exe_path)
        self.nex_path = nex_path
        self.cmd = [self.exe_path, '-n', self.nex_path]
        self.subprocess_kwargs = subprocess_kwargs

    def get_stderr(self):
        if not self.stderr_path:
            return self.stderr
        try:
            se = open(self.stderr_path, 'r')
        except IOError, e:
            _LOG.error('Could not open stderr file')
            raise e
        msg = se.read()
        se.close()
        return msg

    def start(self):
        _LOG.debug('Starting process with following command:\n\t'
                '{0}'.format(' '.join(self.cmd)))
        if self.stdout_path:
            sout = open(self.stdout_path, 'w')
        else:
            sout = subprocess.PIPE
        if self.stderr_path:
            serr = open(self.stderr_path, 'w')
        else:
            serr = subprocess.PIPE
        self.process = subprocess.Popen(self.cmd,
                stdout = sout,
                stderr = serr,
                shell = False,
                **self.subprocess_kwargs)
        self.stdout, self.stderr = self.process.communicate()
        self.exit_code = self.process.wait()
        if hasattr(sout, 'close'):
            sout.close()
        if hasattr(serr, 'close'):
            serr.close()
        # if self.exit_code != 0:
        #         _LOG.error('execution failed')
        #         raise WorkerExecutionError('{0} failed.\ninvocation:\n{1}\n'
        #                 'stderr:\n{2}'.format(
        #                 self.name, ' '.join([str(x) for x in self.cmd]),
        #                 self.get_stderr()))
        self.finished = True

def check_for_output(path):
    if os.path.exists(path):
        raise Exception('script has already been run: output files already '
                'exist... aborting')

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
    seq_file_paths = []
    for i, align_len in enumerate(alignment_lengths):
        locus = i
        num_locus_copies = 2
        sp_tree_path = project_util.SP_TREE_PATH
        if i == 0:
            locus = 'mt'
            num_locus_copies = 1
            # Using a 4-fold faster mutation rate for a mitochondrial locus
            # (all branches of species tree are multiplied by 4) so that the
            # difference in ploidy cancels out. E.g., 2Nu = 0.5N4u
            sp_tree_path = project_util.SP_TREE_MT_PATH
        tree_path = os.path.join(project_util.GENE_TREE_DIR,
                'gene-tree-locus-{0}.nex'.format(locus))
        seq_path = os.path.join(project_util.SEQ_DIR,
                'locus-{0}.nex'.format(locus))
        check_for_output(tree_path)
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
        seq_file_paths.append(seq_path)
    tree_workers = manager.Manager.run_workers(tree_workers,
            num_processors = multiprocessing.cpu_count())
    seq_workers = manager.Manager.run_workers(seq_workers,
            num_processors = multiprocessing.cpu_count())
    return seq_file_paths


def parse_alignments(paths,
        to_cull = {
                '1': {'mt': (0, 0, 0),
                      '1':  (0, 1, 11),
                      '2':  (1, 2, 0),
                      '3':  (2, 1, 26),
                      '4':  (1, 0, 5),
                      '5':  (1, 1, 33)},
                '2': {'mt': (1, 0, 51),
                      '1':  (2, 0, 0),
                      '3':  (0, 1, 0),
                      '4':  (1, 1, 0),
                      '5':  (0, 0, 0)},
                '3': {'mt': (0, 1, 13),
                      '1':  (0, 2, 33),
                      '3':  (1, 0, 9),
                      '4':  (2, 1, 17)},
                },
        ):
    temp_fs = tempfs.TempFileSystem(
            parent = project_util.BIN_DIR)
    if not os.path.exists(project_util.LIZARD_SEQ_DIR):
        functions.mkdr(project_util.LIZARD_SEQ_DIR)
    if not os.path.exists(project_util.LIZARD_CONFIG_DIR):
        functions.mkdr(project_util.LIZARD_CONFIG_DIR)
    config_path = os.path.join(project_util.LIZARD_CONFIG_DIR, 'sample-table.txt')
    check_for_output(config_path)
    config_lines = []
    pi_path = os.path.join(project_util.LIZARD_CONFIG_DIR, 'theta-estimates.txt')
    check_for_output(pi_path)
    pi_stream = open(pi_path, 'w')
    pi_stream.write('species\tpopulation\tlocus\tpi\n')
    pi_summary = stats.SampleSummarizer()
    id_pattern = re.compile(
            r'^species(?P<species>\d)_pop(?P<pop>\d)_tip(?P<tip>\d)$')
    data_sets = {}
    for p in paths:
        num_locus_copies = 2
        locus = os.path.splitext(p)[0].split('-')[-1]
        if locus.lower() == 'mt':
            num_locus_copies = 1
        sequences = dataio.get_seq_dict(p, format = 'nexus')
        for seq_id, seq in sequences.iteritems():
            m = id_pattern.match(seq_id)
            sp = m.group('species')
            pop = m.group('pop')
            if not sp in data_sets:
                data_sets[sp] = {}
            if not pop in data_sets[sp]:
                data_sets[sp][pop] = {}
            data_sets[sp][pop][seq_id] = seq
        for species in sorted(data_sets.keys()):
            try:
                cull_tup = to_cull[species][locus]
            except KeyError:
                continue
            pop_dict = data_sets[species]
            seq_list = []
            parameters = {}
            nsamples = []
            for population, seq_dict in pop_dict.iteritems():
                if cull_tup[2] > 0:
                    for s in seq_dict.itervalues():
                        s.seq = s.seq[:-cull_tup[2]]
                seq_ids = sorted(seq_dict.keys())
                cull_index = (cull_tup[int(population) - 1] *
                        num_locus_copies)
                if cull_index > 0:
                    seq_ids = seq_ids[:-cull_index]
                nsamples.append(len(seq_ids))
                seq_list.extend((seq_dict[k] for k in seq_ids))
                pi = seqstats.average_number_of_pairwise_differences(
                    seq_iter = (seq_dict[k] for k in seq_ids),
                    aligned = True,
                    per_site = True)
                pi_summary.add_sample(pi)
                pi_stream.write('\t{0}\t{1}\t{2}\t{3}\n'.format(species,
                        population, locus, pi))
            fasta_path = os.path.join(project_util.LIZARD_SEQ_DIR,
                    'species-{0}-locus-{1}.fasta'.format(species, locus))
            check_for_output(fasta_path)
            dataio.write_seqs(seq_list,
                    dest = fasta_path,
                    format = 'fasta')
            tmp_nex_data_path = temp_fs.get_file_path()
            tmp_nex_exe_path = temp_fs.get_file_path()
            nseqs = dataio.convert_format(in_file = fasta_path,
                    out_file = tmp_nex_data_path,
                    in_format = 'fasta',
                    out_format = 'nexus')
            write_paup_hky_file(tmp_nex_exe_path, tmp_nex_data_path)
            stdout_path = temp_fs.get_file_path()
            pw = PaupWorker(nex_path = tmp_nex_exe_path,
                    stdout_path = stdout_path)
            pw.start()
            parameters = parse_paup_log_file(stdout_path)
            parameters['species'] = 'species-{0}'.format(species)
            parameters['locus'] = 'locus-{0}'.format(locus)
            if locus.lower() == 'mt':
                parameters['ploidy_multiplier'] = 0.25
                parameters['rate_multiplier'] = 4.0
            else:
                parameters['ploidy_multiplier'] = 1.0
                parameters['rate_multiplier'] = 1.0
            parameters['nsamples1'] = nsamples[0]
            parameters['nsamples2'] = nsamples[1]
            assert parameters['nsamples1'] + parameters['nsamples2'] == nseqs
            parameters['nsites'] = len(pop_dict['1'].values()[0].seq)
            parameters['path'] = os.path.relpath(fasta_path,
                    os.path.dirname(config_path))
            if ((float(parameters['kappa']) < 1.0) or
                    (float(parameters['kappa']) > 1000.0)):
                parameters['kappa'] = 1.0
            config_lines.append('{species}\t{locus}\t{ploidy_multiplier}\t'
                    '{rate_multiplier}\t{nsamples1}\t{nsamples2}\t{kappa}\t'
                    '{nsites}\t{a}\t{c}\t{g}\t{path}\n'.format(**parameters))
    sys.stdout.write('Summary of pi estimates:\n{0}\n'.format(str(pi_summary)))
    config_lines.sort()
    with open(config_path, 'w') as out:
        for l in config_lines:
            out.write(l)
    pi_stream.close()
    temp_fs.purge()

def main():
    paths = simulate_alignments()
    parse_alignments(paths)


if __name__ == '__main__':
    main()

