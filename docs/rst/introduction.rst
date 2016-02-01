Introduction
============

The package provides a Matlab implementation of the Spectral Ewald fast Ewald summation
method (SE) for six different kernels:

- ``SE``:          3-periodic Laplace (electrostatics)
- ``SE2P``:        2-periodic Laplace (electrostatics)
- ``SE_Stokes``:    3-periodic Stokeslet (viscous flow)
- ``SE2P_Stokes``:  2-periodic Stokeslet (viscous flow)	
- ``SE_Stresslet``: 3-periodic Stresslet (viscous flow)
- ``SE_Rotlet``:    3-periodic Rotlet    (viscous flow)	

Most of the computational kernels are written in optimized C for speed, and could very
well be interfaced to another language by mimicking the MEX interface provided.

Getting started
---------------
Each directory contains an m-file that implements the method, and a
basic accuracy/convergence script that shows how to use it. From any
one of these directories, do: 

.. code-block:: matlab

   init    % sets up paths 
   make    % builds all necessary mex C code 

then depending on directory you can do

.. code-block:: matlab

   demo    % Runs demo of spectral convergence + error estimates
           % for both real and Fourier space sums.

or

.. code-block:: matlab

   test_accuracy    % should display plot of spectral convergence

Most recent development has focused on 3P Stokes flow, so the directories related to that (SE_Stokes, SE_Rotlet, SE_Stresslet) are more developed.

Code examples for various are kernels can be found by looking at the tests, located in the
folder ``mfile_tests`` or ``tests`` in the folder for the respective kernel.

Testing
-------

Some directory contain a ``tests`` folder with unit tests. To run the full test suite in that folder simply execute

.. code-block:: matlab

    init
    run_unit_tests

To run a all full suite with all available tests, go to the root directory and run

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

The FGG and direct summation C implementations (which are substantial)
are better than any other versions scattered in the original
implementations of the four SE methods. Only bare-bones
implementations of SE methods included in this package. Associated
methods for e.g. fast real-space summation are not included. These can
be found in the original implementation directories, and may need to
be integrated before they work. Original implementations include
Matlab-implementations of FGG and direct summation, which are not
part of this package.

Notes on fast Gaussian gridding
-------------------------------

1. Common mex wrappers are included, compiled from application makefiles (above) or from
   the matlab verification tests found in 'Xp_matlab_impl'

2. Stand-alone C test code found in 'testing', useful for debugging and checking for memory leaks.

3. The FGG C-code contains comments that are meant to be helpful.

4. By default, manually implemented SSE2 implementations are ENABLED, see comments in  SE_fgg.c

5. There are remnants of OpenMP parallelization to be found throughout.  It is easy to add  appropriate work-sharing loop directives, particularly in the from-grid part. The most  mature OpenMP can be found in the original SE2P implementation (not in this package)
