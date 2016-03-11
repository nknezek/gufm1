metadata = dict( name= 'gufm1',
                 version = '0.1',
                 description='Functions to calculate historial magnetic field for the Earth using GUFM1',
                 url='https://github.com/nknezek/gufm1',
                 author='Nicholas Knezek',
                 author_email='nknezek@gmail.com',
                 license='MIT',
                 long_description='',
                 packages = ['gufm1'],
                 package_dir = {'gufm1': 'gufm1'},
                 data_files = [('gufm1_data',['data/gufm1_data.txt'])]
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