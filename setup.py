from setuptools import setup, find_packages
import sys, os

version = '0.1'

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
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
