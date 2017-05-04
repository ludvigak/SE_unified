Introduction
============

The package provides a Matlab implementation of the :doc:`Spectral Ewald <math>` method
(SE) for fast Ewald summation related to six different kernels:

- :mod:`SE`:          3-periodic Laplace (electrostatics)
- :mod:`SE2P`:        2-periodic Laplace (electrostatics)
- :mod:`SE1P`:        1-periodic Laplace (electrostatics)
- :mod:`SE_Stokes`:    3-periodic Stokeslet (viscous flow)
- :mod:`SE2P_Stokes`:  2-periodic Stokeslet (viscous flow)	
- :mod:`SE_Stresslet`: 3-periodic Stresslet (viscous flow)
- :mod:`SE_Rotlet`:    3-periodic Rotlet    (viscous flow)	

The aim of this package is to provide fast routines for computing the Fourier space sums
for the above kernels. Routines for the real space sums are provided for some kernels.

Most of the computational kernels are written in C, and quite some time has been spent
optimizing them using explicit SIMD instructions and OpenMP shared-memory parallelism. The
interfacing MEX layer is generally kept thin, so interfacing to another language than
Matlab should be quite straightforward.

Building
--------

Most of the C/MEX code is now built with CMake. To build it open a terminal and do

.. code-block:: shell

   cd build
   cmake ..
   make

To run the test suite:

.. code-block:: shell

   make test

Some of the C/MEX code is not yet integrated into CMake, and has to be built using Matlab
build script. The build script is then called ``make.m`` and is located in the base folder
for that kernel.

Getting started
---------------

Each directory contains the m-files related to that specific kernel. From any of these
directories do

.. code-block:: matlab

   init    % sets up paths 

then depending on directory you can do

.. code-block:: matlab

   demo    % Runs demo of spectral convergence + error estimates
           % for both real and Fourier space sums.

or

.. code-block:: matlab
		
   make             % build C/MEX code for this directory
   test_accuracy    % should display plot of spectral convergence

Example (in Matlab):

.. code-block:: matlab
   
   cd SE_Rotlet
   init
   demo

Most recent development has focused on 3P Stokes flow, so the directories related to that (SE_Stokes, SE_Rotlet, SE_Stresslet) are more developed.

Code examples for various are kernels can be found by looking at the tests, located in the
folder ``mfile_tests`` or ``tests`` in the folder for the respective kernel.

Testing
-------

Some directories contain a ``tests`` folder with unit tests. To run that test suite simply
execute

.. code-block:: matlab

    init
    run_unit_tests

To run the full suite, go to the root directory and run

.. code-block:: matlab

    run_unit_tests

This is recommended to do after building or before committing changes.

Additional files
----------------
The package also contains:

* ``util``: Common functions
* ``SE_fast_gridding``: C implementation of fast Gaussian gridding (below)
* ``SE_direct``: C-code for direct Ewald sums for Laplace 2P/3P
* ``SE_Stokes_direct``: C-code for direct Ewald sums for Stokes 2P/P3
* ``SE_leftovers``: Spectral Ewald, fast real-space and k=0 codes for Laplace and
  Stokes. These are unmaintained, but might prove useful.
