# The Spectral Ewald Unified package

This is a Matlab package that implements the Spectral Ewald (SE) fast 
Ewald summation method for six different kernels:

     - SE:          3-periodic Laplace (electrostatics)
     - SE2P:        2-periodic Laplace (electrostatics)
     - SE_Stokes    3-periodic Stokeslet (viscous flow)
     - SE2P_Stokes  2-periodic Stokeslet (viscous flow)	
     - SE_Stresslet 3-periodic Stresslet (viscous flow)
     - SE_Rotlet    3-periodic Rotlet    (viscous flow)	

## Getting started
Each directory contains an m-file that implements the method, and a
basic accuracy/convergence script that shows how to use it. From any
one of these directories, do: 

```matlab
init    % sets up paths 
make    % builds all necessary mex C code 
```
then depending on directory you can do
```matlab
test_accuracy    % should see plot of spectral convergence
```
or
```matlab
test_all    % run a larger test suite
```

Most recent development has focused on 3P Stokes flow, so the directories related to that are more developed.

## Additional files
The package also contains:

1. `util`: Common functions
2. `SE_fast_gridding`: C implementation of fast Gaussian gridding (below)
3. `SE_direct`: C-code for direct Ewald sums for Laplace 2P/3P
4. `SE_Stokes_direct`: C-code for direct Ewald sums for Stokes 2P/P3

The FGG and direct summation C implementations (which are substantial)
are better than any other versions scattered in the original
implementations of the four SE methods. Only bare-bones
implementations of SE methods included in this package. Associated
methods for e.g. fast real-space summation are not included. These can
be found in the original implementation directories, and may need to
be integrated before they work. Original implementations include
Matlab-implementations of FGG and direct summation, which are not
part of this package.

## Notes on fast Gaussian gridding

1. Common mex wrappers are included, compiled from application
  makefiles (above) or from the matlab verification tests found in
  'Xp_matlab_impl'

2. Stand-alone C test code found in 'testing', useful for debugging and
  checking for memory leaks.

3. The FGG C-code contains comments that are meant to be helpful.

4. By default, manually implemented SSE2 implementations are ENABLED,
  see comments in SE_fgg.c

5. There are remnants of OpenMP parallelization to be found throughout.
  It is easy to add appropriate work-sharing loop directives,
  particularly in the from-grid part. The most mature OpenMP can be 
  found in the original SE2P implementation (not in this package)

## License
This package is released under the MIT License, see the [LICENSE](./LICENSE) file for details.
