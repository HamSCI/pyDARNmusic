#!/usr/bin/env python

from distutils.core import setup

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(name='pyDARNmusic',
      version='0.1',
      description='pyDARN-Compatible SuperDARN Multiple Signal Classification Algorithm for Detection of MSTIDs',
      author='Nathaniel A. Frissell',
      author_email='nathaniel.frissell@scranton.edu',
      url='https://hamsci.org',
      packages=['pyDARNmusic'],
      install_requires=requirements
     )
