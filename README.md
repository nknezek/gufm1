GUFM1
=====
This package provides Python code to calculate magnetic fields of Andrew Jackson's GUFM1 historical magnetic field model.

Details of GUFM1 can be found here: <http://jupiter.ethz.ch/~cfinlay/gufm1.html>

This code was largely based off of the Fortran code *eval_field.f* found on his website, but is completely re-written in native Python.

Usage
-----
To use, simply clone the git directory, then *import gufm1* to begin using the methods.

Users should be able to be calculate fields for any year from 1590 to 1990 from the main gufm1_data.txt file, but you can also use data files from individual years downloaded directly from his website if desired.
