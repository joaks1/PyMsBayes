import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages
import sys, os
import platform
import subprocess

from pymsbayes import __version__

DEV_MODE = False
if 'develop' in sys.argv:
    DEV_MODE = True

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
SCRIPT_NAMES = [
        'dmc_dpp_summary.py',
        'dmc_plot_results.py',
        'dmc_posterior_probs.py',
        'dmc.py',
        'dmc_estimate_prior_probs.py',
        'dmc_plot_stats.py',
        # 'dmc_prob_shared_divergence.py',
        'dmc_saturation_plot.py',
        'dmc_sum_sims.py',
        ]
TOOL_NAMES = [
        'ABCestimator',
        'dpp-msbayes.pl',
        'dpp-msprior',
        'eureject',
        'msbayes.pl',
        'msDQH',
        'msprior',
        'msReject',
        'obsSumStats.pl',
        'regress_cli.r',
        'sumstatsvector',
        ]

def get_scripts_to_install():
    script_dir = os.path.relpath(os.path.join(BASE_DIR, 'scripts'), BASE_DIR)
    return [os.path.join(script_dir, f) for f in SCRIPT_NAMES]

def check_tool(path):
    try:
        p = subprocess.Popen([path],
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE)
        p.terminate()
    except subprocess.CalledProcessError:
        return path
    except OSError:
        return None
    return path

def get_tools_to_install():
    warning_msg = '''
**********************************************************************
****************************** WARNING *******************************
The bundled `dpp-msbayes` tools are not being installed, because some
of them are not executable on this system. The `PyMsBayes` package and
scripts will still be installed, however, you will need to build and
install `dpp-msbayes` (https://github.com/joaks1/dpp-msbayes) yourself
in order to use them.  Sorry for the inconvenience.
**********************************************************************
'''
    from pymsbayes.utils import ToolPathManager
    bin_dir = ToolPathManager._bin_dir
    if os.path.isdir(bin_dir):
        tools = [os.path.join(bin_dir, f) for f in TOOL_NAMES]
        ret = []
        for t in tools:
            if check_tool(t):
                tool_rel_path = os.path.relpath(t, BASE_DIR)
                ret.append(tool_rel_path)
            else:
                sys.stderr.write(warning_msg)
                return []
        return ret
    sys.stderr.write(warning_msg)
    return []

scripts = get_scripts_to_install()
if not DEV_MODE:
    scripts += get_tools_to_install()

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
      packages=find_packages(),
      scripts = scripts,
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

