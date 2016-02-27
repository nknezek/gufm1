metadata = dict( name= 'gufm1',
                 version = 0.1,
                 description='Functions to calculate Historial Magnetic Field for the Earth using GUFM1',
                 url='',
                 author='Nicholas Knezek',
                 author_email='nknezek@berkeley.edu',
                 license='MIT',
                 long_description='',
                 packages = ['gufm1'],
               )

#Try to use setuptools in order to check dependencies.
#if the system does not have setuptools, fall back on
#distutils.
try:
  from setuptools import setup, find_packages, Extension
  metadata['install_requires'] = ['numpy', 'scipy']
except ImportError:
  from distutils.core import setup, Extension

setup ( **metadata )