
Welcome to Spectral Ewald's documentation!
==========================================

This is the documentation for the unified Spectral Ewald package, available on GitHub at
https://github.com/ludvigak/SE_unified . To obtain the lastest version simply execute

.. code-block:: none

    git clone https://github.com/ludvigak/SE_unified.git

The package provides a Matlab implementation of the Spectral Ewald (SE) fast Ewald
summation method for six different kernels:

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

   test_accuracy    % should display plot of spectral convergence

or

.. code-block:: matlab

   test_all    % run a larger test suite

Most recent development has focused on 3P Stokes flow, so the directories related to that are more developed.

Additional files
----------------
The package also contains:

* ``util``: Common functions
* ``SE_fast_gridding``: C implementation of fast Gaussian gridding (below)
* ``SE_direct``: C-code for direct Ewald sums for Laplace 2P/3P
* ``SE_Stokes_direct``: C-code for direct Ewald sums for Stokes 2P/P3

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



Contents
--------
.. toctree::
   :maxdepth: 2

   rst/cinterface.rst
