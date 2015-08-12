Sparse Eigensolvers for Java (**SEJ**) builds on top of [MTJ](https://code.google.com/p/matrix-toolkits-java/) to provide eigenvalue and eigenvector solvers for sparse matrices.

## Status ##

SEJ is currently in ALPHA and the API is subject to change: we hope to create a common API for a large variety of solvers, making it as easy as possible for clients to change the solver backend.

The following sparse solvers have been implemented:

  * [LOBPCG](http://en.wikipedia.org/wiki/LOBPCG)

and the following are planned shortly:

  * [ARPACK](http://www.caam.rice.edu/software/ARPACK/)
  * Block Davidson


## Installation ##

SEJ is a Java library, not a standalone application. Add the SEJ (and MTJ) jars to your build path. If you are getting `java.lang.NoClassDefFoundError` exceptions, you probably forgot to do this step.

When the API becomes stable, SEJ will be distributed via Maven.

## Documentation ##

Clients are advised to link their IDEs to the source code of this project in order to see the Javadocs and to study the test cases for example use cases.

The [LobpcgJavaExamples](LobpcgJavaExamples.md) page shows some use case examples.

## References ##

  1. The Symmetric Eigenvalue Problem, Beresford N. Parlett, Society for  Industrial and Applied Mathematics (SIAM), Philadelphia, PA, 1998. Corrected reprint of the 1980 original.
  1. Toward The Optimal Preconditioned Eigensolver: Locally Optimal Block Preconditioned Conjugate Gradient Method, A. V. Knyazev, SIAM Journal on Scientic Computing 23 (2001), no. 2, pp. 517-541.
  1. Block Locally Optimal Preconditioned Eigenvalue Xolvers (BLOPEX) in hypre and PETSc, A. V. Knyazev, I. Lashuk, M. E. Argentati, and E. Ovchinnikov, SIAM Journal on Scientific Computing (SISC), Vol. 29 (2007),  No. 5, pp. 2224-2239.
  1. A Comparison of Eigensolvers for Large-scale 3D Modal Analysis Using AMG-Preconditioned Iterative Methods, (with Peter Arbenz, Ulrich Hetmaniuk, Ray Tuminaro), International Journal for Numerical Methods in Engineering, Volume 64, pp. 204-236, 2005.
  1. A Survey of Software for Sparse Eigenvalue Problems, V. Hernández, J. E. Román, A. Thomás and V. Vidal, SLEPc Technical Report STR-6, June 2007.