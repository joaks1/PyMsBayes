from setuptools import setup, find_packages
import sys, os
import platform
import subprocess

from pymsbayes import __version__

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
    platform_name = platform.system().lower()
    # bin_dir = None
    # if platform_name == 'linux':
    #     bin_dir = os.path.join(BASE_DIR, "bin", "linux")
    # elif platform_name == 'darwin':
    bin_dir = os.path.join(BASE_DIR, "bin", "mac")
    # elif platform_name == 'windows':
    #     bin_dir = os.path.join(BASE_DIR, "bin", "win")
    # else:
    #     return []
    tools = [os.path.join(bin_dir, f) for f in TOOL_NAMES]
    ret = []
    for t in tools:
        if check_tool(t):
            ret.append(t)
        else:
            return []
    return ret


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
      scripts = get_scripts_to_install() + get_tools_to_install(),
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

