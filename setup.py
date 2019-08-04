#!/usr/bin/env python
from setuptools import setup, find_packages

if __name__ == "__main__":
    setup(name='zchecker',
          version='2.4.3',
          description='ZTF moving target checker for short object lists',
          author="Michael S. P. Kelley",
          author_email="msk@astro.umd.edu",
          url="https://github.com/mkelley/zchecker",
          packages=find_packages(),
          requires=['numpy', 'astropy', 'requests'],
          scripts=['scripts/zchecker', 'scripts/zproject', 'scripts/zstack',
                   'scripts/zphot', 'scripts/ds9-ztf'],
          license='BSD',
          )
