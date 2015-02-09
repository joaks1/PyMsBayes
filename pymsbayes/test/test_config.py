#! /usr/bin/env python

import unittest
import os
from cStringIO import StringIO

from pymsbayes.config import MsBayesConfig, AlignmentConfig, SampleTable
from pymsbayes.utils.probability import *
from pymsbayes.utils import errors
from pymsbayes.test.support.pymsbayes_test_case import PyMsBayesTestCase
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class AlignmentConfigTestCase(unittest.TestCase):
    def test_empty_string(self):
        self.assertRaises(errors.SampleTableRowError, AlignmentConfig, '')

    def test_empty_list(self):
        self.assertRaises(errors.SampleTableRowError, AlignmentConfig, [])

    def test_column_error(self):
        row = "pair0\tlocus0\t1\t1\t10\t10\t15.6\t0.5\t1000\t0.35\t0.28\t0.13\tpair0.locus0.fasta"
        self.assertRaises(errors.SampleTableRowError, AlignmentConfig, row)
        row = "pair0\tlocus0\t1\t1\t10\t10\t1000\t0.35\t0.28\t0.13\tpair0.locus0.fasta"
        self.assertRaises(errors.SampleTableRowError, AlignmentConfig, row)
        row = ["pair0", "locus0", "1", "1", "10", "10", "15.6", "0.5", "1000",
                "0.35", "0.28", "0.13", "pair0.locus0.fasta"]
        self.assertRaises(errors.SampleTableRowError, AlignmentConfig, row)
        row = ["pair0", "locus0", "1", "1", "10", "10", "1000",
                "0.35", "0.28", "0.13", "pair0.locus0.fasta"]
        self.assertRaises(errors.SampleTableRowError, AlignmentConfig, row)

    def test_string_row(self):
        row = ("species-3\tlocus-mt\t0.25\t4.0\t5\t4\t11.423634\t587\t0.227487"
                "\t0.222071\t0.259081\t../sequences/species-3-locus-mt.fasta")
        ac = AlignmentConfig(row)
        self.assertEqual(ac.taxon_name, "species-3")
        self.assertEqual(ac.locus_name, "locus-mt")
        self.assertEqual(ac.ploidy_multiplier, 0.25)
        self.assertEqual(ac.mutation_rate_multiplier, 4.0)
        self.assertEqual(ac.number_of_gene_copies, (5, 4))
        self.assertEqual(ac.kappa, 11.423634)
        self.assertEqual(ac.length, 587)
        self.assertEqual(ac.base_frequencies, (0.227487, 0.222071, 0.259081))
        self.assertEqual(ac.path, "../sequences/species-3-locus-mt.fasta")
        self.assertEqual(str(ac), row)

    def test_tuple_row(self):
        row = ("species-3", "locus-mt", "0.25", "4.0", "5", "4", "11.423634",
                "587", "0.227487", "0.222071", "0.259081",
                "../sequences/species-3-locus-mt.fasta")
        ac = AlignmentConfig(row)
        self.assertEqual(ac.taxon_name, "species-3")
        self.assertEqual(ac.locus_name, "locus-mt")
        self.assertEqual(ac.ploidy_multiplier, 0.25)
        self.assertEqual(ac.mutation_rate_multiplier, 4.0)
        self.assertEqual(ac.number_of_gene_copies, (5, 4))
        self.assertEqual(ac.kappa, 11.423634)
        self.assertEqual(ac.length, 587)
        self.assertEqual(ac.base_frequencies, (0.227487, 0.222071, 0.259081))
        self.assertEqual(ac.path, "../sequences/species-3-locus-mt.fasta")
        self.assertEqual([str(x) for x in ac.get_sample_table_row_elements()],
                list(row))
        self.assertEqual(str(ac), '\t'.join(row))

    def test_generator_row(self):
        r = ("species-3\tlocus-mt\t0.25\t4.0\t5\t4\t11.423634\t587\t0.227487"
                "\t0.222071\t0.259081\t../sequences/species-3-locus-mt.fasta")
        row = (s for s in r.split('\t'))
        ac = AlignmentConfig(row)
        self.assertEqual(ac.taxon_name, "species-3")
        self.assertEqual(ac.locus_name, "locus-mt")
        self.assertEqual(ac.ploidy_multiplier, 0.25)
        self.assertEqual(ac.mutation_rate_multiplier, 4.0)
        self.assertEqual(ac.number_of_gene_copies, (5, 4))
        self.assertEqual(ac.kappa, 11.423634)
        self.assertEqual(ac.length, 587)
        self.assertEqual(ac.base_frequencies, (0.227487, 0.222071, 0.259081))
        self.assertEqual(ac.path, "../sequences/species-3-locus-mt.fasta")
        self.assertEqual(str(ac), r)

    def test_equality(self):
        row1 = ("species-3", "locus-mt", "0.25", "4.0", "5", "4", "11.423634",
                "587", "0.227487", "0.222071", "0.259081",
                "../sequences/species-3-locus-mt.fasta")
        ac1 = AlignmentConfig(row1)
        ac2 = AlignmentConfig(row1)
        self.assertFalse(ac1 == ac2)
        self.assertTrue(ac1.equals(ac2))
        row2 = ("species-3", "locus-mt", "0.25", "4.0", "5", "4", "11.423633",
                "587", "0.227487", "0.222071", "0.259082",
                "../sequences/species-3-locus-mt.fasta")
        ac3 = AlignmentConfig(row2)
        self.assertFalse(ac1.equals(ac3))

class SampleTableTestCase(unittest.TestCase):
    def setUp(self):
        self.rows = [
            ["species-1", "locus-1", "1.0", "1.0", "10", "8", "32.422050",
                    "389", "0.271215", "0.240174", "0.266343",
                    "../sequences/species-1-locus-1.fasta"],
            ["species-1", "locus-2", "1.0", "1.0", "8", "6", "5.507905", "500",
                    "0.254861", "0.225477", "0.246877",
                    "../sequences/species-1-locus-2.fasta"],
            ["species-1", "locus-3", "1.0", "1.0", "6", "8", "8.379708", "524",
                    "0.260506", "0.233269", "0.266142",
                    "../sequences/species-1-locus-3.fasta"],
            ["species-1", "locus-4", "1.0", "1.0", "8", "10", "5.204980",
                    "345", "0.251830", "0.232765", "0.249506",
                    "../sequences/species-1-locus-4.fasta"],
            ["species-1", "locus-5", "1.0", "1.0", "8", "8", "29.592792",
                    "417", "0.272341", "0.237232", "0.210548",
                    "../sequences/species-1-locus-5.fasta"],
            ["species-1", "locus-mt", "0.25", "4.0", "5", "5", "8.153262",
                    "600", "0.222976", "0.242721", "0.271977",
                    "../sequences/species-1-locus-mt.fasta"],
            ["species-2", "locus-1", "1.0", "1.0", "6", "10", "7.536519",
                    "400", "0.256404", "0.246540", "0.266092",
                    "../sequences/species-2-locus-1.fasta"],
            ["species-2", "locus-3", "1.0", "1.0", "10", "8", "11.148510",
                    "550", "0.270202", "0.229906", "0.249895",
                    "../sequences/species-2-locus-3.fasta"],
            ["species-2", "locus-4", "1.0", "1.0", "8", "8", "9.391906", "350",
                    "0.246659", "0.242283", "0.237685",
                    "../sequences/species-2-locus-4.fasta"],
            ["species-2", "locus-5", "1.0", "1.0", "10", "10", "13.327843",
                    "450", "0.264189", "0.240497", "0.227266",
                    "../sequences/species-2-locus-5.fasta"],
            ["species-2", "locus-mt", "0.25", "4.0", "4", "5", "7.595008",
                    "549", "0.233664", "0.264141", "0.234924",
                    "../sequences/species-2-locus-mt.fasta"],
            ["species-3", "locus-1", "1.0", "1.0", "10", "6", "17.035406",
                    "367", "0.258149", "0.231107", "0.276950",
                    "../sequences/species-3-locus-1.fasta"],
            ["species-3", "locus-3", "1.0", "1.0", "8", "10", "59.177467",
                    "541", "0.262631", "0.225555", "0.251191",
                    "../sequences/species-3-locus-3.fasta"],
            ["species-3", "locus-4", "1.0", "1.0", "6", "8", "6.901196", "333",
                    "0.287292", "0.230559", "0.215738",
                    "../sequences/species-3-locus-4.fasta"],
            ["species-3", "locus-mt", "0.25", "4.0", "5", "4", "11.423634",
                    "587", "0.227487", "0.211233", "0.259081",
                    "../sequences/species-3-locus-mt.fasta"],
        ]
        self.row_dict = None
        self.taxa = None
        self.loci = None
        self.stream = None
        self.update_stream()

    def update_stream(self, start = "BEGIN SAMPLE_TBL", end = "END SAMPLE_TBL"):
        self.row_dict = {}
        for r in self.rows:
            if not r[0] in self.row_dict:
                self.row_dict[r[0]] = {}
                self.row_dict[r[0]][r[1]] = r
                continue
            if not r[1] in self.row_dict[r[0]]:
                self.row_dict[r[0]][r[1]] = r
                continue

        self.taxa = []
        self.loci = []
        for r in self.rows:
            if not r[0] in self.taxa:
                self.taxa.append(r[0])
            if not r[1] in self.loci:
                self.loci.append(r[1])


        self.table_string = "\n".join([
                start,
                "\n".join(["\t".join(r) for r in self.rows]),
                end])
        self.stream = StringIO()
        self.stream.write("{0}\n".format(self.table_string))
        self.stream.seek(0)

    def test_empty_string(self):
        self.assertRaises(IOError, SampleTable, '')

    def test_empty_buffer(self):
        s = StringIO()
        self.assertRaises(errors.SampleTableError, SampleTable, s)

    def test_bad_row(self):
        self.rows[-1].pop(-1)
        self.update_stream()
        self.assertRaises(errors.SampleTableRowError, SampleTable, self.stream)

    def test_no_start(self):
        self.update_stream(start="")
        self.assertRaises(errors.SampleTableError, SampleTable, self.stream)

    def test_no_end(self):
        self.update_stream(end="")
        self.assertRaises(errors.SampleTableError, SampleTable, self.stream)

    def test_simple(self):
        st = SampleTable(self.stream)
        self.assertIsInstance(st, SampleTable)
        self.assertEqual(st.taxa, self.taxa)
        self.assertEqual(st.loci, self.loci)
        self.assertEqual(st.npairs, len(self.taxa))
        self.assertEqual(str(st), self.table_string)

    def test_equality(self):
        st1 = SampleTable(self.stream)
        self.update_stream()
        st2 = SampleTable(self.stream)
        self.assertFalse(st1 == st2)
        self.assertTrue(st1.equals(st2))
        self.rows[1][-2] = "0.25"
        self.update_stream()
        st3 = SampleTable(self.stream)
        self.assertFalse(st1.equals(st3))
        self.rows.pop(-1)
        self.update_stream()
        st4 = SampleTable(self.stream)
        self.assertFalse(st4.equals(st3))

    def test_duplicate_taxon_locus_pair(self):
        self.rows[-1][0] = 'species-1'
        self.rows[-1][1] = 'locus-1'
        self.update_stream()
        self.assertRaises(errors.SampleTableError, SampleTable, self.stream)

class MsBayesConfigTestCase(PyMsBayesTestCase):
    
    def setUp(self):
        self.cfg = StringIO()

    def _update_config(self, cfg, params, multi_locus=False, new_impl=False,
            time_in_subs_per_site = False):
        if not time_in_subs_per_site:
            if new_impl:
                cfg.write("""concentrationShape = {c_shape}
concentrationScale = {c_scale}
thetaShape = {theta_shape}
thetaScale = {theta_scale}
ancestralThetaShape = {a_theta_shape}
ancestralThetaScale = {a_theta_scale}
thetaParameters = {theta_params}
tauShape = {tau_shape}
tauScale = {tau_scale}
migrationShape = {migration_shape}
migrationScale = {migration_scale}
recombinationShape = {recombination_shape}
recombinationScale = {recombination_scale}
numTauClasses = {psi}
constrain = 0
subParamConstrain = 111111111
""".format(**params))
            else:
                cfg.write("""upperTheta = 0.1
lowerTheta = {ltheta}
upperTau = {utau}
numTauClasses = {psi}
upperMig = {umig}
upperRec = {urec}
upperAncPopSize = {atheta}
constrain = 0
subParamConstrain = 111111111
""".format(**params))
        else:
            if new_impl:
                cfg.write("""concentrationShape = {c_shape}
concentrationScale = {c_scale}
thetaShape = {theta_shape}
thetaScale = {theta_scale}
ancestralThetaShape = {a_theta_shape}
ancestralThetaScale = {a_theta_scale}
thetaParameters = {theta_params}
tauShape = {tau_shape}
tauScale = {tau_scale}
timeInSubsPerSite = 1
migrationShape = {migration_shape}
migrationScale = {migration_scale}
recombinationShape = {recombination_shape}
recombinationScale = {recombination_scale}
numTauClasses = {psi}
constrain = 0
subParamConstrain = 111111111
""".format(**params))
            else:
                cfg.write("""upperTheta = 0.1
lowerTheta = {ltheta}
upperTau = {utau}
timeInSubsPerSite = 1
numTauClasses = {psi}
upperMig = {umig}
upperRec = {urec}
upperAncPopSize = {atheta}
constrain = 0
subParamConstrain = 111111111
""".format(**params))
        if not multi_locus:
            cfg.write("""
BEGIN SAMPLE_TBL
pair0	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair0.locus0.fasta
pair1	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair2.locus0.fasta
pair2	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair3.locus0.fasta
pair3	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair3.locus0.fasta
END SAMPLE_TBL
""")
        else:
            cfg.write("""
BEGIN SAMPLE_TBL
pair0	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair0.locus0.fasta
pair0	locus1	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair0.locus0.fasta
pair1	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair2.locus0.fasta
pair1	locus1	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair2.locus0.fasta
pair1	locus2	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair2.locus0.fasta
pair1	locus3	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair2.locus0.fasta
pair2	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair3.locus0.fasta
pair3	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair3.locus0.fasta
pair3	locus1	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair3.locus0.fasta
END SAMPLE_TBL
""")
        cfg.seek(0)

    def test_is_config(self):
        p = {'ltheta': 0.0001,
             'utheta': 0.1,
             'utau': 10.0,
             'psi': 0,
             'umig': 0.0,
             'urec': 0.0,
             'atheta': 1.0,}
        self._update_config(self.cfg, p)
        self.assertTrue(MsBayesConfig.is_config(self.cfg))
        cfg = StringIO()
        cfg.write('not an\nmsbayes\nconfig file\nblah\n\n')
        self.assertFalse(MsBayesConfig.is_config(cfg))

    def test_single_locus(self):
        p = {'ltheta': 0.0001,
             'utheta': 0.1,
             'utau': 10.0,
             'psi': 0,
             'umig': 0.0,
             'urec': 0.0,
             'atheta': 1.0,}
        self._update_config(self.cfg, p)
        _LOG.debug('testing config:\n\n{0}\n'.format(self.cfg.getvalue()))
        c = MsBayesConfig(self.cfg)
        self.assertEqual(c.npairs, 4)
        self.assertIsInstance(c.psi, DiscreteUniformDistribution)
        self.assertIsInstance(c.tau, ContinuousUniformDistribution)
        self.assertIsInstance(c.recombination, ContinuousUniformDistribution)
        self.assertIsInstance(c.migration, ContinuousUniformDistribution)
        self.assertIsInstance(c.a_theta, ContinuousUniformDistribution)
        self.assertIsInstance(c.theta, ContinuousUniformDistribution)
        self.assertIsInstance(c.d_theta, BetaDistribution)
        self.assertEqual(c.div_model_prior, 'psi')
        self.assertEqual(c.dpp_concentration, None)
        self.assertEqual(c.theta_parameters, None)
        self.assertEqual(c.implementation, 'old')
        self.assertFalse(c.time_in_subs_per_site)

        self.assertSameDistributions(c.psi, DiscreteUniformDistribution(1,4))
        self.assertSameDistributions(c.tau, ContinuousUniformDistribution(0,10))
        self.assertSameDistributions(c.recombination, ContinuousUniformDistribution(0,0))
        self.assertSameDistributions(c.migration, ContinuousUniformDistribution(0,0))
        self.assertSameDistributions(c.a_theta, ContinuousUniformDistribution(0.0001,0.1))
        self.assertSameDistributions(c.theta, ContinuousUniformDistribution(0.0001,0.1))
        self.assertSameDistributions(c.d_theta, BetaDistribution(1,1,2))

    def test_single_locus_round_trip(self):
        p = {'ltheta': 0.0001,
             'utheta': 0.1,
             'utau': 10.0,
             'psi': 0,
             'umig': 0.0,
             'urec': 0.0,
             'atheta': 1.0,}
        self._update_config(self.cfg, p)
        _LOG.debug('testing config:\n\n{0}\n'.format(self.cfg.getvalue()))
        c = MsBayesConfig(self.cfg)
        s = StringIO()
        c.write(s)
        _LOG.debug('written config:\n\n{0}\n'.format(s.getvalue()))
        s.seek(0)
        c2 = MsBayesConfig(s)
        self.assertSameConfigs([c, c2])

    def test_single_locus_time_in_subs(self):
        p = {'ltheta': 0.0001,
             'utheta': 0.1,
             'utau': 10.0,
             'psi': 0,
             'umig': 0.0,
             'urec': 0.0,
             'atheta': 1.0,}
        self._update_config(self.cfg, p, time_in_subs_per_site = True)
        _LOG.debug('testing config:\n\n{0}\n'.format(self.cfg.getvalue()))
        c = MsBayesConfig(self.cfg)
        self.assertEqual(c.npairs, 4)
        self.assertIsInstance(c.psi, DiscreteUniformDistribution)
        self.assertIsInstance(c.tau, ContinuousUniformDistribution)
        self.assertIsInstance(c.recombination, ContinuousUniformDistribution)
        self.assertIsInstance(c.migration, ContinuousUniformDistribution)
        self.assertIsInstance(c.a_theta, ContinuousUniformDistribution)
        self.assertIsInstance(c.theta, ContinuousUniformDistribution)
        self.assertIsInstance(c.d_theta, BetaDistribution)
        self.assertEqual(c.div_model_prior, 'psi')
        self.assertEqual(c.dpp_concentration, None)
        self.assertEqual(c.theta_parameters, None)
        self.assertEqual(c.implementation, 'old')
        self.assertTrue(c.time_in_subs_per_site)

        self.assertSameDistributions(c.psi, DiscreteUniformDistribution(1,4))
        self.assertSameDistributions(c.tau, ContinuousUniformDistribution(0,10))
        self.assertSameDistributions(c.recombination, ContinuousUniformDistribution(0,0))
        self.assertSameDistributions(c.migration, ContinuousUniformDistribution(0,0))
        self.assertSameDistributions(c.a_theta, ContinuousUniformDistribution(0.0001,0.1))
        self.assertSameDistributions(c.theta, ContinuousUniformDistribution(0.0001,0.1))
        self.assertSameDistributions(c.d_theta, BetaDistribution(1,1,2))

    def test_single_locus_time_in_subs_round_trip(self):
        p = {'ltheta': 0.0001,
             'utheta': 0.1,
             'utau': 10.0,
             'psi': 0,
             'umig': 0.0,
             'urec': 0.0,
             'atheta': 1.0,}
        self._update_config(self.cfg, p, time_in_subs_per_site = True)
        _LOG.debug('testing config:\n\n{0}\n'.format(self.cfg.getvalue()))
        c = MsBayesConfig(self.cfg)
        s = StringIO()
        c.write(s)
        _LOG.debug('written config:\n\n{0}\n'.format(s.getvalue()))
        s.seek(0)
        c2 = MsBayesConfig(s)
        self.assertTrue(c.time_in_subs_per_site)
        self.assertTrue(c2.time_in_subs_per_site)
        self.assertSameConfigs([c, c2])

    def test_multi_locus(self):
        p = {'ltheta': 0.0001,
             'utheta': 0.1,
             'utau': 10.0,
             'psi': 0,
             'umig': 0.0,
             'urec': 0.0,
             'atheta': 1.0,}
        self._update_config(self.cfg, p, multi_locus=True)
        _LOG.debug('testing config:\n\n{0}\n'.format(self.cfg.getvalue()))
        c = MsBayesConfig(self.cfg)
        self.assertEqual(c.npairs, 4)
        self.assertIsInstance(c.psi, DiscreteUniformDistribution)
        self.assertIsInstance(c.tau, ContinuousUniformDistribution)
        self.assertIsInstance(c.recombination, ContinuousUniformDistribution)
        self.assertIsInstance(c.migration, ContinuousUniformDistribution)
        self.assertIsInstance(c.a_theta, ContinuousUniformDistribution)
        self.assertIsInstance(c.theta, ContinuousUniformDistribution)
        self.assertIsInstance(c.d_theta, BetaDistribution)
        self.assertEqual(c.div_model_prior, 'psi')
        self.assertEqual(c.dpp_concentration, None)
        self.assertEqual(c.theta_parameters, None)
        self.assertEqual(c.implementation, 'old')
        self.assertFalse(c.time_in_subs_per_site)

        self.assertSameDistributions(c.psi, DiscreteUniformDistribution(1,4))
        self.assertSameDistributions(c.tau, ContinuousUniformDistribution(0,10))
        self.assertSameDistributions(c.recombination, ContinuousUniformDistribution(0,0))
        self.assertSameDistributions(c.migration, ContinuousUniformDistribution(0,0))
        self.assertSameDistributions(c.a_theta, ContinuousUniformDistribution(0.0001,0.1))
        self.assertSameDistributions(c.theta, ContinuousUniformDistribution(0.0001,0.1))
        self.assertSameDistributions(c.d_theta, BetaDistribution(1,1,2))

    def test_multi_locus_round_trip(self):
        p = {'ltheta': 0.0001,
             'utheta': 0.1,
             'utau': 10.0,
             'psi': 0,
             'umig': 0.0,
             'urec': 0.0,
             'atheta': 1.0,}
        self._update_config(self.cfg, p, multi_locus=True)
        _LOG.debug('testing config:\n\n{0}\n'.format(self.cfg.getvalue()))
        c = MsBayesConfig(self.cfg)
        s = StringIO()
        c.write(s)
        _LOG.debug('written config:\n\n{0}\n'.format(s.getvalue()))
        s.seek(0)
        c2 = MsBayesConfig(s)
        self.assertSameConfigs([c, c2])

    def test_new_implementation(self):
        p = {'c_shape': 2,
             'c_scale': 5,
             'theta_shape': 2.0,
             'theta_scale': 0.002,
             'tau_shape': 1.0,
             'tau_scale': 5.0,
             'a_theta_shape': 1,
             'a_theta_scale': 0.01,
             'theta_params': '001',
             'migration_shape': 1.5,
             'migration_scale': 0.1,
             'recombination_shape': 0.0,
             'recombination_scale': 0.0,
             'psi': 0,}
        self._update_config(self.cfg, p, new_impl=True)
        _LOG.debug('testing config:\n\n{0}\n'.format(self.cfg.getvalue()))
        c = MsBayesConfig(self.cfg)
        self.assertEqual(c.npairs, 4)
        self.assertEqual(c.psi, None)
        self.assertIsInstance(c.tau, GammaDistribution)
        self.assertIsInstance(c.recombination, ContinuousUniformDistribution)
        self.assertIsInstance(c.migration, GammaDistribution)
        self.assertIsInstance(c.a_theta, GammaDistribution)
        self.assertEqual(c.theta, None)
        self.assertIsInstance(c.d_theta, GammaDistribution)
        self.assertEqual(c.div_model_prior, 'dpp')
        self.assertIsInstance(c.dpp_concentration, GammaDistribution)
        self.assertEqual(c.theta_parameters, [0, 0, 1])
        self.assertEqual(c.implementation, 'new')
        self.assertFalse(c.time_in_subs_per_site)

        self.assertSameDistributions(c.dpp_concentration, GammaDistribution(2, 5))
        self.assertSameDistributions(c.tau, GammaDistribution(1, 5))
        self.assertSameDistributions(c.recombination, ContinuousUniformDistribution(0,0))
        self.assertSameDistributions(c.migration, GammaDistribution(1.5, 0.1))
        self.assertSameDistributions(c.a_theta, GammaDistribution(1, 0.01))
        self.assertSameDistributions(c.d_theta, GammaDistribution(2, 0.002))

    def test_new_implementation_round_trip(self):
        p = {'c_shape': 2,
             'c_scale': 5,
             'theta_shape': 2.0,
             'theta_scale': 0.002,
             'tau_shape': 1.0,
             'tau_scale': 5.0,
             'a_theta_shape': 1,
             'a_theta_scale': 0.01,
             'theta_params': '001',
             'migration_shape': 1.5,
             'migration_scale': 0.1,
             'recombination_shape': 0.0,
             'recombination_scale': 0.0,
             'psi': 0,}
        self._update_config(self.cfg, p, new_impl=True)
        _LOG.debug('testing config:\n\n{0}\n'.format(self.cfg.getvalue()))
        c = MsBayesConfig(self.cfg)
        s = StringIO()
        c.write(s)
        _LOG.debug('written config:\n\n{0}\n'.format(s.getvalue()))
        s.seek(0)
        c2 = MsBayesConfig(s)
        self.assertSameConfigs([c, c2])

    def test_equal_sample_table(self):
        p1 = {'c_shape': 2,
             'c_scale': 5,
             'theta_shape': 2.0,
             'theta_scale': 0.002,
             'tau_shape': 1.0,
             'tau_scale': 5.0,
             'a_theta_shape': 1,
             'a_theta_scale': 0.01,
             'theta_params': '001',
             'migration_shape': 1.5,
             'migration_scale': 0.1,
             'recombination_shape': 0.0,
             'recombination_scale': 0.0,
             'psi': 0,}
        p2 = {'c_shape': 1000,
             'c_scale': 0.003,
             'theta_shape': 1.0,
             'theta_scale': 0.02,
             'tau_shape': 2.0,
             'tau_scale': 2.0,
             'a_theta_shape': 1,
             'a_theta_scale': 0.01,
             'theta_params': '001',
             'migration_shape': 1.5,
             'migration_scale': 0.1,
             'recombination_shape': 0.0,
             'recombination_scale': 0.0,
             'psi': 0,}
        cfg1 = StringIO()
        self._update_config(cfg1, p1, multi_locus=True, new_impl=True)
        c1 = MsBayesConfig(cfg1)
        cfg2 = StringIO()
        self._update_config(cfg2, p2, multi_locus=True, new_impl=True)
        c2 = MsBayesConfig(cfg2)
        cfg3 = StringIO()
        self._update_config(cfg3, p2, multi_locus=False, new_impl=True)
        c3 = MsBayesConfig(cfg3)
        self.assertTrue(c1.equal_sample_table(c2))
        self.assertFalse(c1.equal_sample_table(c3))


if __name__ == '__main__':
    unittest.main()

