metadata = dict( name= 'gufm1',
                 version = '0.1',
                 description='Functions to calculate historial magnetic field for the Earth using GUFM1',
                 url='https://github.com/nknezek/gufm1',
                 author='Nicholas Knezek',
                 author_email='nknezek@gmail.com',
                 license='MIT',
                 long_description='Functions to calculate historial magnetic field for the Earth using GUFM1',
                 packages = ['gufm1'],
                 package_dir = {'gufm1': 'gufm1'},
                 data_files = [('gufm1/data',['gufm1/data/gufm1_data.txt'])],
                 install_requires = ['numpy', 'scipy']
               )


from setuptools import setup, find_packages, Extension

setup ( **metadata )