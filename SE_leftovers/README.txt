Spectral Ewald, fast RS and k0 codes for Laplace and Stokes

Notes:

1) Laplace 

a) 3P. Code for fast real-space (RS) summation using a combined Linked
   Cell List and Verlet method. Neighbour list generated explicitly,
   and has to be stored in memory (this is a limiting
   factor). Evaluation kernel is written in C, but it can be improved.

b) 2P. Code for fast evaluation of k0 sum using Chebyshev polynomials.

   No fast RS code. Adapting the 3P fast RS code is nearly trivial --
   just modify the part that does one-layer periodic extension from
   all three directions to just two.

2) Stokes

a) 3P. Fast RS code virtually indentical to 3P Laplace.

   The evalation kernel is over the Stokeslet, but neighbour list and
   pruning code is the same. Unifying all fast RS codes (four
   combinations Stokes/Laplace in 2P/3P) makes a lot of sense, but has
   not been done yet.

b) 2P. Code for fast evaluation of k0 sum using Chebyshev
   polynomials. Two versions: normal and matrix form. In the latter,
   all Stokeslet evauations are taken into a matrix, whose elements on
   particle-particle distances but not on the Stoekslet stengths,
   f. Thus, if particles are not moving, the k0 sum is obtained by
   normal matrix-vector multiplication.

   Also, code for a matrix-form of the real-space sum. Real-space
   interactions are assembled in to a large sparse matrix with 3x3
   block structure mirroring the decomposed Stokeslet tensor. After
   assembling this matrix, the real-space sum is computed as a single
   sparse matvec with the strengths f, which is _very_ fast. As with
   the matrix form of the k0 sum, the idea is to not have to
   reevaluate the costly Stokeslet kernel (erf and distance
   calculation) if particles have not moved.

   The assembly process is not "fast" (i.e. it's N^2). It could be
   fused with the LCL code, but this has not been done yet. Also, the
   assembled matrix takes a lot of memory. Here there's room for
   improvement. Firstly, the symmetry of the Stokeslet is not
   exploited. So only six of the nine blocks need to be computed and
   stored. Secondly, each block has the same non-zero structure. This
   means that the sparsity pattern (which is the same thing as the
   neighbour list) only needs to be stored once (instead of nine times
   as today). The reduction in memory on offer here is significant.

See also: SE_Unified for plain Ewald sum implementations, spectral
          k-space methods and Fast Gridding code.

Dag Lindbo, dag@kth.se, 2011
