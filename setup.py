from setuptools import setup, find_packages
import sys, os

from pymsbayes import __version__

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
SCRIPTS_DIR = os.path.relpath(os.path.join(BASE_DIR, 'scripts'), BASE_DIR)
SCRIPT_NAMES = [
        'dmc_dpp_summary.py',
        'dmc_plot_results.py',
        'dmc_posterior_probs.py',
        'dmc.py',
        'dmc_estimate_prior_probs.py',
        'dmc_plot_stats.py',
        # 'dmc_prob_shared_divergence.py',
        'dmc_saturation_plot.py',
        ]
SCRIPTS = [os.path.join(SCRIPTS_DIR, f) for f in SCRIPT_NAMES]

setup(name='PyMsBayes',
      version=__version__,
      description="Python msBayes wrapper",
      long_description="""\
Python msBayes wrapper""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords="msBayes ABC Approximate Bayesian Computation",
      author="Jamie Oaks",
      author_email="joaks1@gmail.com",
      url='',
      license="GPL",
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      scripts=SCRIPTS,
      include_package_data=True,
      zip_safe=False,
      test_suite="pymsbayes.test.get_unittest_suite",
      # test_suite="pymsbayes.test",
      # test_loader="unittest:TestLoader",
      install_requires=[
          'configobj'
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )

