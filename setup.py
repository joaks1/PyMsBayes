from setuptools import setup, find_packages
import sys, os

import pymsbayes.utils

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

version = '0.1'

def symlink_msbayes_requirements():
    for f in ['msDQH', 'sumstatsvector']:
        match = None
        for fname in os.listdir(os.path.join(pymsbayes.utils.BIN_DIR, 'new')):
            if fname.startswith(f):
                if match:
                    sys.stderr.write('ERROR: found mulitple {0!r}\n'.format(f))
                    sys.exit(1)
                match = fname
        if not match:
            sys.stderr.write('ERROR: could not find {0!r}\n'.format(f))
            sys.exit(1)
        src_path = os.path.join(pymsbayes.utils.BIN_DIR, 'new', match)
        dest_path = os.path.join(pymsbayes.utils.BIN_DIR, 'old', match)
        sys.stderr.write("\nCreating link: {0!r} => {1!r}\n".format(src_path,
                dest_path))
        if os.path.exists(dest_path):
            if os.path.islink(dest_path):
                real_dest = os.path.abspath(os.path.realpath(dest_path))
                if real_dest != os.path.abspath(os.path.realpath(src_path)):
                    sys.stderr.write('ERROR: Symbolic link {0!r} already '
                            'exists, but points to different source '
                            '{1!r}\n'.format(dest_path, real_dest))
                    sys.exit(1)
                else:
                    sys.stderr.write("Correct link already exists\n")
            else:
                sys.stderr.write('ERROR: Link target path {0!r} already '
                        'exists\n'.format(dest_path))
                sys.exit(1)
        else:
            os.symlink(src_path, dest_path)

setup(name='PyMsBayes',
      version=version,
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

